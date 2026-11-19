use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};

use rust_htslib::bam::{self, Record};

use crate::alignment_io::{self, AlignmentFormat, AlignmentOutput};
use crate::alignment_sort::RecordOutput;
use crate::barcode::{Barcode, BarcodeError, BarcodeExtractor};
use crate::clustering::{MixedUmiLengths, cluster_umis};
use crate::dedup::{
    DedupMethod, GroupKey, PositionOptions, PythonRandom, five_prime_position, get_read_position,
    sampled_barcode,
};
use crate::gene::{Flush, GeneAssigner, GeneError, GeneOptions};
use crate::pairing::PairingOptions;

#[allow(clippy::struct_excessive_bools)]
pub struct GroupConfig {
    pub method: DedupMethod,
    pub ignore_umi: bool,
    pub barcode: BarcodeExtractor,
    pub umi_group_tag: Vec<u8>,
    pub random_seed: u64,
    pub output_path: Option<String>,
    pub output_format: AlignmentFormat,
    pub reference: Option<String>,
    pub output_bam: bool,
    pub no_sort_output: bool,
    /// Directory for the spills of the output sort.
    pub temp_dir: Option<String>,
    pub chrom: Option<String>,
    pub group_out: Option<String>,
    pub edit_distance_threshold: u32,
    pub position: PositionOptions,
    pub subset: Option<f64>,
    pub mapping_quality: u8,
    pub buffer_whole_contig: bool,
    pub gene: GeneOptions,
    pub pairing: PairingOptions,
    pub ignore_tlen: bool,
}

pub struct GroupStats {
    pub input_reads: u64,
    pub output_reads: u64,
}

struct GroupSlot {
    records: Vec<Record>,
    count: u32,
    insertion_order: u32,
}

type Drained<K> = BTreeMap<K, BTreeMap<GroupKey, HashMap<Vec<u8>, GroupSlot>>>;

struct GroupBuffer<K: Ord = i64> {
    groups: Drained<K>,
}

impl GroupBuffer<i64> {
    fn drain_up_to(&mut self, threshold: i64) -> Drained<i64> {
        let rest = self.groups.split_off(&(threshold + 1));
        std::mem::replace(&mut self.groups, rest)
    }
}

impl<K: Ord + Clone> GroupBuffer<K> {
    const fn new() -> Self {
        Self {
            groups: BTreeMap::new(),
        }
    }

    fn add(&mut self, record: Record, pos: K, key: GroupKey, umi: Vec<u8>) {
        let umi_map = self.groups.entry(pos).or_default().entry(key).or_default();

        if let Some(slot) = umi_map.get_mut(&umi) {
            slot.count += 1;
            slot.records.push(record);
            return;
        }

        let order = u32::try_from(umi_map.len()).unwrap_or(u32::MAX);
        umi_map.insert(
            umi,
            GroupSlot {
                records: vec![record],
                count: 1,
                insertion_order: order,
            },
        );
    }

    fn drain_key(&mut self, key: &K) -> Drained<K> {
        self.groups
            .remove(key)
            .map_or_default(|key_map| BTreeMap::from([(key.clone(), key_map)]))
    }

    fn drain_all(&mut self) -> Drained<K> {
        std::mem::take(&mut self.groups)
    }
}

fn assign_groups(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, GroupSlot>,
    edit_threshold: u32,
) -> Result<Vec<Vec<Vec<u8>>>, MixedUmiLengths> {
    Ok(cluster_umis(
        method,
        umi_map
            .iter()
            .map(|(umi, slot)| (umi.as_slice(), slot.count, slot.insertion_order)),
        edit_threshold,
    )?
    .into_iter()
    .map(|group| group.into_iter().map(<[u8]>::to_vec).collect())
    .collect())
}

/// Process drained position groups: assign UMI groups, annotate records, write TSV rows.
#[allow(clippy::cast_sign_loss)]
fn process_drained<K: Ord>(
    drained: Drained<K>,
    config: &GroupConfig,
    unique_id: &mut u32,
    tsv_writer: &mut Option<BufWriter<File>>,
    header_view: &bam::HeaderView,
    assigner: Option<&GeneAssigner>,
    output: &mut RecordOutput,
) -> Result<(), GroupError> {
    for key_map in drained.into_values() {
        for (_, mut umi_map) in key_map {
            let groups = assign_groups(config.method, &umi_map, config.edit_distance_threshold)?;

            for group in &groups {
                let top_umi = &group[0];
                let group_count: u32 = group.iter().map(|u| umi_map[u].count).sum();
                let top_umi_str = std::str::from_utf8(top_umi).unwrap_or("");

                for umi in group {
                    let slot = umi_map.remove(umi).expect("UMI must exist in umi_map");

                    for record in slot.records {
                        if let Some(w) = tsv_writer.as_mut() {
                            let read_name = std::str::from_utf8(record.qname()).unwrap_or("");
                            let contig_name =
                                std::str::from_utf8(header_view.tid2name(record.tid() as u32))
                                    .unwrap_or("");
                            let umi_str = std::str::from_utf8(umi).unwrap_or("");
                            let read_pos = five_prime_position(&record);
                            // umi_tools prints the contig, not the mapped gene, for --per-contig.
                            let gene_label = match assigner {
                                Some(assigner) if assigner.per_contig() => contig_name.to_string(),
                                Some(assigner) => assigner.gene(&record).map_or_default(|gene| {
                                    String::from_utf8_lossy(&gene).into_owned()
                                }),
                                None => "NA".to_string(),
                            };

                            writeln!(
                                w,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                read_name,
                                contig_name,
                                read_pos,
                                gene_label,
                                umi_str,
                                slot.count,
                                top_umi_str,
                                group_count,
                                *unique_id,
                            )
                            .map_err(|e| GroupError::TsvWrite(e.to_string()))?;
                        }

                        let mut tagged = record;
                        if config.output_bam {
                            alignment_io::set_aux(
                                &mut tagged,
                                b"UG",
                                rust_htslib::bam::record::Aux::U32(*unique_id),
                            )
                            .map_err(|e| GroupError::BamWrite(e.to_string()))?;
                            alignment_io::set_aux(
                                &mut tagged,
                                &config.umi_group_tag,
                                rust_htslib::bam::record::Aux::String(top_umi_str),
                            )
                            .map_err(|e| GroupError::BamWrite(e.to_string()))?;
                        }
                        output
                            .push(tagged)
                            .map_err(|e| GroupError::BamWrite(e.to_string()))?;
                    }
                }

                *unique_id += 1;
            }
        }
    }

    Ok(())
}

/// # Errors
///
/// Returns `GroupError` on BAM I/O failures or unknown chromosome filter.
#[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
pub fn run_group(config: &GroupConfig, input_path: &str) -> Result<GroupStats, GroupError> {
    let (mut source, assigner) = GeneAssigner::open(
        &config.gene,
        config.chrom.as_deref(),
        input_path,
        config.reference.as_deref(),
    )?;
    let mut flusher = assigner.as_ref().map(GeneAssigner::flusher);
    let header_view = source.header().clone();

    let writer = if config.output_bam {
        let header = if config.no_sort_output {
            bam::Header::from_template(&header_view)
        } else {
            alignment_io::coordinate_sorted_header(&header_view)
        };
        Some(
            alignment_io::open_writer(
                &header,
                AlignmentOutput {
                    path: config.output_path.as_deref(),
                    format: config.output_format,
                    reference: config.reference.as_deref(),
                },
            )
            .map_err(|e| GroupError::BamWrite(e.to_string()))?,
        )
    } else {
        None
    };

    // Optional chromosome filter
    let chrom_filter: Option<i32> = config
        .chrom
        .as_ref()
        .map(|c| {
            let tid = source
                .header()
                .tid(c.as_bytes())
                .ok_or_else(|| GroupError::UnknownChrom(c.clone()))?;
            #[allow(clippy::cast_possible_wrap)]
            Ok::<i32, GroupError>(tid as i32)
        })
        .transpose()?;

    // Open TSV writer
    let mut tsv_writer: Option<BufWriter<File>> = config
        .group_out
        .as_ref()
        .map(|path| {
            let file =
                File::create(path).map_err(|e| GroupError::TsvWrite(e.to_string()))?;
            let mut w = BufWriter::new(file);
            writeln!(
                w,
                "read_id\tcontig\tposition\tgene\tumi\tumi_count\tfinal_umi\tfinal_umi_count\tunique_id"
            )
            .map_err(|e| GroupError::TsvWrite(e.to_string()))?;
            Ok::<_, GroupError>(w)
        })
        .transpose()?;

    let mut buffer = GroupBuffer::<i64>::new();
    let mut gene_buffer = GroupBuffer::<Vec<u8>>::new();
    let mut stats = GroupStats {
        input_reads: 0,
        output_reads: 0,
    };

    #[allow(clippy::cast_possible_truncation)]
    let mut rng = PythonRandom::new(config.random_seed);

    let mut output = RecordOutput::new(
        writer,
        (!config.no_sort_output).then(|| alignment_io::coordinate_sorted_header(&header_view)),
        config.temp_dir.as_deref(),
    );
    let mut unique_id: u32 = 0;

    let mut last_start: i64 = 0;
    let mut last_chrom: i32 = -1;

    while let Some(record) = source
        .read_next()
        .map_err(|e| GroupError::BamRead(e.to_string()))?
    {
        let tid = record.tid();
        if chrom_filter.is_some_and(|filter_tid| tid != filter_tid) {
            continue;
        }

        let triage = config.pairing.triage(&record, true);
        for _ in 0..triage.copies {
            output
                .push(record.clone())
                .map_err(|e| GroupError::BamWrite(e.to_string()))?;
        }
        if triage.is_read2 {
            continue;
        }
        stats.input_reads += 1;
        if !triage.grouped {
            continue;
        }

        let Some(Barcode { umi, cell }) = sampled_barcode(
            &record,
            config.subset,
            &mut rng,
            config.mapping_quality,
            &config.barcode,
            config.ignore_umi,
        )?
        else {
            continue;
        };

        if let (Some(assigner), Some(flusher)) = (&assigner, flusher.as_mut()) {
            let Some(gene) = assigner.gene(&record) else {
                continue;
            };
            let done = match flusher.before_read(tid) {
                Flush::None => None,
                Flush::All => Some(gene_buffer.drain_all()),
                Flush::Gene(done) => Some(gene_buffer.drain_key(&done)),
            };
            if let Some(drained) = done {
                process_drained(
                    drained,
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    Some(assigner),
                    &mut output,
                )?;
            }
            flusher.after_read(tid, &gene);

            let key: GroupKey = GroupKey::for_gene(cell);
            gene_buffer.add(record, gene, key, umi);
        } else {
            // Standard coordinate mode
            let position = get_read_position(&record, config.position.soft_clip_threshold);
            let start = position.start;

            if tid != last_chrom {
                process_drained(
                    buffer.drain_all(),
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    assigner.as_ref(),
                    &mut output,
                )?;
            } else if !config.buffer_whole_contig && start > last_start + 1000 {
                let threshold = start - 1000;
                process_drained(
                    buffer.drain_up_to(threshold),
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    assigner.as_ref(),
                    &mut output,
                )?;
            }

            last_start = start;
            last_chrom = tid;

            let key = GroupKey::for_position(
                &record,
                &position,
                &config.position,
                config.pairing.paired && !config.ignore_tlen,
                cell,
            );

            buffer.add(record, position.pos, key, umi);
        }
    }

    process_drained(
        buffer.drain_all(),
        config,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
        assigner.as_ref(),
        &mut output,
    )?;
    process_drained(
        gene_buffer.drain_all(),
        config,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
        assigner.as_ref(),
        &mut output,
    )?;

    // Flush TSV
    if let Some(w) = tsv_writer.as_mut() {
        w.flush().map_err(|e| GroupError::TsvWrite(e.to_string()))?;
    }

    stats.output_reads = output
        .finish()
        .map_err(|e| GroupError::BamWrite(e.to_string()))?;

    Ok(stats)
}

#[derive(Debug, thiserror::Error)]
pub enum GroupError {
    #[error("failed to open BAM: {0}")]
    BamOpen(String),
    #[error("failed to read BAM record: {0}")]
    BamRead(String),
    #[error("failed to write BAM/SAM: {0}")]
    BamWrite(String),
    #[error("failed to write TSV: {0}")]
    TsvWrite(String),
    #[error("unknown chromosome: {0}")]
    UnknownChrom(String),
    #[error(transparent)]
    Barcode(#[from] BarcodeError),
    #[error(transparent)]
    Gene(#[from] GeneError),
    #[error(transparent)]
    UmiLengths(#[from] MixedUmiLengths),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(name: &str) -> Record {
        let view = bam::HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:10000\n");
        let line = format!("{name}\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI");
        Record::from_sam(&view, line.as_bytes()).unwrap()
    }

    fn buffer(reads: &[(&str, i64, &str, &str)]) -> GroupBuffer<i64> {
        let mut buffer = GroupBuffer::new();
        for &(name, position, cell, umi) in reads {
            let key = GroupKey::for_gene(cell.as_bytes().to_vec());
            buffer.add(record(name), position, key, umi.as_bytes().to_vec());
        }
        buffer
    }

    #[test]
    fn umis_are_numbered_in_first_seen_order_within_a_position_and_key() {
        let mut buffer = buffer(&[
            ("r1", 100, "", "AAAA"),
            ("r2", 100, "", "CCCC"),
            ("r3", 100, "", "AAAA"),
            ("r4", 100, "cell", "GGGG"),
            ("r5", 200, "", "TTTT"),
            ("r6", 100, "", "GGGG"),
        ]);
        let drained = buffer.drain_all();
        let slots = &drained[&100][&GroupKey::for_gene(Vec::new())];
        let summary = |umi: &[u8]| {
            let slot = &slots[umi];
            let names: Vec<&[u8]> = slot.records.iter().map(Record::qname).collect();
            (slot.insertion_order, slot.count, names.concat())
        };
        assert_eq!(summary(b"AAAA"), (0, 2, b"r1r3".to_vec()));
        assert_eq!(summary(b"CCCC"), (1, 1, b"r2".to_vec()));
        assert_eq!(summary(b"GGGG"), (2, 1, b"r6".to_vec()));
        let other_cell = &drained[&100][&GroupKey::for_gene(b"cell".to_vec())];
        assert_eq!(other_cell[&b"GGGG"[..]].insertion_order, 0);
        assert_eq!(
            drained[&200].values().next().unwrap()[&b"TTTT"[..]].insertion_order,
            0
        );
        assert!(buffer.drain_all().is_empty());
    }

    #[test]
    fn drain_up_to_includes_the_threshold_position() {
        let mut buffer = buffer(&[
            ("r1", 99, "", "AAAA"),
            ("r2", 100, "", "AAAA"),
            ("r3", 101, "", "AAAA"),
        ]);
        assert_eq!(
            buffer.drain_up_to(100).keys().copied().collect::<Vec<_>>(),
            [99, 100]
        );
        assert!(buffer.drain_up_to(100).is_empty());
        assert_eq!(
            buffer.drain_all().keys().copied().collect::<Vec<_>>(),
            [101]
        );
    }

    #[test]
    fn drain_key_takes_one_gene_and_leaves_the_rest() {
        let mut buffer = GroupBuffer::<Vec<u8>>::new();
        for (name, gene) in [("r1", "g1"), ("r2", "g2"), ("r3", "g1")] {
            let key = GroupKey::for_gene(Vec::new());
            buffer.add(
                record(name),
                gene.as_bytes().to_vec(),
                key,
                b"AAAA".to_vec(),
            );
        }
        let drained = buffer.drain_key(&b"g1".to_vec());
        assert_eq!(drained.len(), 1);
        assert_eq!(
            drained[&b"g1".to_vec()].values().next().unwrap()[&b"AAAA"[..]].count,
            2
        );
        assert!(buffer.drain_key(&b"g1".to_vec()).is_empty());
        assert_eq!(buffer.drain_all().len(), 1);
    }
}
