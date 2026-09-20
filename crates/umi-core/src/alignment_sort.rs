//! Stable coordinate sorting with bounded output buffers.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use rust_htslib::bam::{self, Read as _, Record};
use tempfile::TempPath;

use crate::alignment_io::coordinate_sort_key;

#[derive(Debug, thiserror::Error)]
pub enum SortError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Alignment(#[from] rust_htslib::errors::Error),
}

struct Run {
    path: TempPath,
    count: u64,
    order: u64,
}

/// Counts emitted reads even for TSV-only grouping, which has no alignment sink.
pub struct RecordOutput {
    writer: Option<bam::Writer>,
    sorter: Option<RecordSorter>,
    count: u64,
}

impl RecordOutput {
    pub fn new(writer: Option<bam::Writer>, sorted_header: Option<bam::Header>) -> Self {
        let sorter = if writer.is_some() {
            sorted_header.map(RecordSorter::new)
        } else {
            None
        };
        Self {
            writer,
            sorter,
            count: 0,
        }
    }

    pub fn push(&mut self, record: Record) -> Result<(), SortError> {
        if let Some(sorter) = &mut self.sorter {
            sorter.push(record)?;
        } else if let Some(writer) = &mut self.writer {
            writer.write(&record)?;
        }
        self.count += 1;
        Ok(())
    }

    pub fn finish(mut self) -> Result<u64, SortError> {
        if let (Some(sorter), Some(writer)) = (self.sorter.take(), self.writer.as_mut()) {
            sorter.finish(|record| Ok(writer.write(record)?))?;
        }
        Ok(self.count)
    }
}

pub struct RecordSorter {
    header: bam::Header,
    buffer: Vec<Record>,
    bytes: usize,
    levels: Vec<Vec<Run>>,
    byte_limit: usize,
    fan_in: usize,
    next_order: u64,
}

impl RecordSorter {
    pub const fn new(header: bam::Header) -> Self {
        Self {
            header,
            buffer: Vec::new(),
            bytes: 0,
            levels: Vec::new(),
            byte_limit: 64 * 1024 * 1024,
            fan_in: 64,
            next_order: 0,
        }
    }

    pub fn push(&mut self, record: Record) -> Result<(), SortError> {
        self.bytes += record.inner().m_data as usize + std::mem::size_of::<Record>();
        self.buffer.push(record);
        if self.bytes >= self.byte_limit {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> Result<(), SortError> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer.sort_by_key(coordinate_sort_key);
        let path = tempfile::NamedTempFile::new()?.into_temp_path();
        let mut writer = bam::Writer::from_path(&path, &self.header, bam::Format::Bam)?;
        for record in &self.buffer {
            writer.write(record)?;
        }
        drop(writer);
        let mut run = Run {
            path,
            count: self.buffer.len() as u64,
            order: self.next_order,
        };
        self.next_order += 1;
        self.buffer.clear();
        self.bytes = 0;
        let mut level = 0;
        loop {
            if level == self.levels.len() {
                self.levels.push(Vec::new());
            }
            self.levels[level].push(run);
            if self.levels[level].len() < self.fan_in {
                break;
            }
            let runs = std::mem::take(&mut self.levels[level]);
            run = self.merge_to_run(runs)?;
            level += 1;
        }
        Ok(())
    }

    fn merge_to_run(&self, runs: Vec<Run>) -> Result<Run, SortError> {
        let path = tempfile::NamedTempFile::new()?.into_temp_path();
        let mut writer = bam::Writer::from_path(&path, &self.header, bam::Format::Bam)?;
        let order = runs[0].order;
        let count = merge_runs(&runs, |record| Ok(writer.write(record)?))?;
        drop(writer);
        drop(runs);
        Ok(Run { path, count, order })
    }

    pub fn finish(
        mut self,
        mut emit: impl FnMut(&Record) -> Result<(), SortError>,
    ) -> Result<(), SortError> {
        if self.levels.is_empty() {
            self.buffer.sort_by_key(coordinate_sort_key);
            for record in &self.buffer {
                emit(record)?;
            }
            return Ok(());
        }
        self.spill()?;
        let mut runs: Vec<_> = std::mem::take(&mut self.levels)
            .into_iter()
            .flatten()
            .collect();
        runs.sort_by_key(|run| run.order);
        while runs.len() > self.fan_in {
            let mut next = Vec::new();
            let mut remaining = runs.into_iter();
            loop {
                let batch: Vec<_> = remaining.by_ref().take(self.fan_in).collect();
                if batch.is_empty() {
                    break;
                }
                next.push(self.merge_to_run(batch)?);
            }
            runs = next;
        }
        merge_runs(&runs, emit)?;
        Ok(())
    }
}

fn read_next(reader: &mut bam::Reader) -> Result<Option<Record>, SortError> {
    let mut record = Record::new();
    Ok(reader.read(&mut record).transpose()?.map(|()| record))
}

/// Run index breaks equal-coordinate ties, preserving original insertion order.
fn merge_runs(
    runs: &[Run],
    mut emit: impl FnMut(&Record) -> Result<(), SortError>,
) -> Result<u64, SortError> {
    let mut readers = Vec::new();
    let mut records = Vec::new();
    let mut heap = BinaryHeap::new();
    for (index, run) in runs.iter().enumerate() {
        let mut reader = bam::Reader::from_path(&run.path)?;
        let record = read_next(&mut reader)?;
        if let Some(record) = &record {
            heap.push(Reverse((coordinate_sort_key(record), index)));
        }
        records.push(record);
        readers.push(reader);
    }
    let mut count = 0;
    while let Some(Reverse((_, index))) = heap.pop() {
        let record = records[index].take().expect("heap entries have a record");
        emit(&record)?;
        count += 1;
        records[index] = read_next(&mut readers[index])?;
        if let Some(record) = &records[index] {
            heap.push(Reverse((coordinate_sort_key(record), index)));
        }
    }
    if count != runs.iter().map(|run| run.count).sum::<u64>() {
        return Err(std::io::Error::other("incomplete alignment sort spill file").into());
    }
    Ok(count)
}
#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (bam::Header, Vec<Record>) {
        let view = bam::HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:10000\n");
        let records = (0..80)
            .map(|i| {
                let (flag, contig, pos) = if i % 9 == 0 {
                    (4, "*", 0)
                } else {
                    (if i % 3 == 0 { 16 } else { 0 }, "chr1", (i % 7) * 10 + 1)
                };
                let cigar = if flag == 4 { "*" } else { "1M" };
                let line =
                    format!("r{i}\t{flag}\t{contig}\t{pos}\t60\t{cigar}\t*\t0\t0\tA\tI\tXX:Z:keep");
                Record::from_sam(&view, line.as_bytes()).unwrap()
            })
            .collect();
        (bam::Header::from_template(&view), records)
    }

    #[test]
    fn spilled_sort_matches_stable_in_memory_sort_and_cleans_up() {
        let (header, mut expected) = fixture();
        let mut sorter = RecordSorter::new(header);
        sorter.byte_limit = 1;
        sorter.fan_in = 3;
        for record in &expected {
            sorter.push(record.clone()).unwrap();
            assert!(
                sorter.buffer.is_empty(),
                "oversized records must be spilled immediately"
            );
        }
        let paths: Vec<_> = sorter
            .levels
            .iter()
            .flatten()
            .map(|run| run.path.to_path_buf())
            .collect();
        assert!(!paths.is_empty());
        expected.sort_by_key(coordinate_sort_key);
        let mut actual = Vec::new();
        sorter
            .finish(|record| {
                actual.push(record.clone());
                Ok(())
            })
            .unwrap();
        assert_eq!(actual, expected);
        assert!(paths.iter().all(|path| !path.exists()));
    }

    #[test]
    fn empty_and_in_memory_sorts_need_no_spill_files() {
        let (header, mut records) = fixture();
        let empty = RecordSorter::new(header.clone());
        empty
            .finish(|_| panic!("empty sort emitted a record"))
            .unwrap();
        let mut sorter = RecordSorter::new(header);
        for record in &records {
            sorter.push(record.clone()).unwrap();
        }
        assert!(sorter.levels.is_empty());
        records.sort_by_key(coordinate_sort_key);
        let mut actual = Vec::new();
        sorter
            .finish(|record| {
                actual.push(record.clone());
                Ok(())
            })
            .unwrap();
        assert_eq!(actual, records);
    }

    #[test]
    fn incomplete_spill_is_reported() {
        let (header, records) = fixture();
        let mut sorter = RecordSorter::new(header.clone());
        sorter.byte_limit = 1;
        sorter.push(records[0].clone()).unwrap();
        let path = sorter.levels[0][0].path.to_path_buf();
        drop(bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap());
        let error = sorter.finish(|_| Ok(())).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("incomplete alignment sort spill")
        );
        assert!(!path.exists());
    }

    #[test]
    fn destination_errors_propagate_and_remove_spills() {
        let (header, records) = fixture();
        let mut sorter = RecordSorter::new(header);
        sorter.byte_limit = 1;
        for record in records {
            sorter.push(record).unwrap();
        }
        let paths: Vec<_> = sorter
            .levels
            .iter()
            .flatten()
            .map(|run| run.path.to_path_buf())
            .collect();
        let result = sorter.finish(|_| Err(std::io::Error::other("injected write failure").into()));
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("injected write failure")
        );
        assert!(paths.iter().all(|path| !path.exists()));
    }
}
