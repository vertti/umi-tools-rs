use std::collections::{BTreeMap, HashMap};
use std::io;

use rust_htslib::bam::{self, Read as BamRead, Record};

/// Trait for RNG used in reservoir-sampling tie-breaks.
///
/// Currently implemented by `PythonRandom` (MT19937 matching `CPython`) to get
/// identical output for compat tests. Can be swapped for any fast RNG once
/// exact-match testing is no longer needed.
trait TieBreakRng {
    /// Return a float in `[0, 1)`.
    fn random(&mut self) -> f64;
}

/// Mersenne Twister 19937 PRNG, matching `CPython`'s random module exactly.
///
/// Python `umi_tools` uses seeded `random.random()` for reservoir-sampling
/// tie-breaks in read selection. We replicate the identical float sequence.
struct PythonRandom {
    mt: [u32; 624],
    index: usize,
}

impl PythonRandom {
    const N: usize = 624;
    const M: usize = 397;
    const MATRIX_A: u32 = 0x9908_b0df;
    const UPPER_MASK: u32 = 0x8000_0000;
    const LOWER_MASK: u32 = 0x7fff_ffff;

    /// Seed the same way `CPython` `random.seed(int)` does:
    /// `init_genrand(19_650_218)` then `init_by_array(&[seed])`.
    fn new(seed: u32) -> Self {
        let mut rng = Self::init_genrand(19_650_218);
        rng.init_by_array(&[seed]);
        rng
    }

    #[allow(clippy::cast_possible_truncation)]
    fn init_genrand(seed: u32) -> Self {
        let mut mt = [0u32; Self::N];
        mt[0] = seed;
        for i in 1..Self::N {
            mt[i] = 1_812_433_253u32
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 30))
                .wrapping_add(i as u32); // i < 624, fits u32
        }
        Self { mt, index: Self::N }
    }

    #[allow(clippy::cast_possible_truncation)]
    fn init_by_array(&mut self, key: &[u32]) {
        let mut i: usize = 1;
        let mut j: usize = 0;
        let k = Self::N.max(key.len());
        for _ in 0..k {
            self.mt[i] = (self.mt[i]
                ^ ((self.mt[i - 1] ^ (self.mt[i - 1] >> 30)).wrapping_mul(1_664_525)))
            .wrapping_add(key[j])
            .wrapping_add(j as u32); // j < key.len(), fits u32
            i += 1;
            j += 1;
            if i >= Self::N {
                self.mt[0] = self.mt[Self::N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..Self::N - 1 {
            self.mt[i] = (self.mt[i]
                ^ ((self.mt[i - 1] ^ (self.mt[i - 1] >> 30)).wrapping_mul(1_566_083_941)))
            .wrapping_sub(i as u32); // i < 624, fits u32
            i += 1;
            if i >= Self::N {
                self.mt[0] = self.mt[Self::N - 1];
                i = 1;
            }
        }
        self.mt[0] = Self::UPPER_MASK;
    }

    fn generate(&mut self) {
        static MAG01: [u32; 2] = [0, PythonRandom::MATRIX_A];
        for kk in 0..Self::N - Self::M {
            let y = (self.mt[kk] & Self::UPPER_MASK) | (self.mt[kk + 1] & Self::LOWER_MASK);
            self.mt[kk] = self.mt[kk + Self::M] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
        }
        for kk in Self::N - Self::M..Self::N - 1 {
            let y = (self.mt[kk] & Self::UPPER_MASK) | (self.mt[kk + 1] & Self::LOWER_MASK);
            self.mt[kk] = self.mt[kk + Self::M - Self::N] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
        }
        let y = (self.mt[Self::N - 1] & Self::UPPER_MASK) | (self.mt[0] & Self::LOWER_MASK);
        self.mt[Self::N - 1] = self.mt[Self::M - 1] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
        self.index = 0;
    }

    fn next_u32(&mut self) -> u32 {
        if self.index >= Self::N {
            self.generate();
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
}

impl TieBreakRng for PythonRandom {
    /// `CPython` `genrand_res53`: 53-bit precision float in `[0, 1)`.
    fn random(&mut self) -> f64 {
        let a = self.next_u32() >> 5;
        let b = self.next_u32() >> 6;
        (f64::from(a) * 67_108_864.0 + f64::from(b)) * (1.0 / 9_007_199_254_740_992.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DedupMethod {
    Unique,
    Percentile,
    Cluster,
    Adjacency,
    Directional,
}

pub struct DedupConfig {
    pub method: DedupMethod,
    pub ignore_umi: bool,
    pub umi_separator: u8,
    pub random_seed: u64,
    pub out_sam: bool,
    pub chrom: Option<String>,
}

pub struct DedupStats {
    pub input_reads: u64,
    pub output_reads: u64,
    pub positions: u64,
}

/// Length of a trailing/leading soft-clip, or 0 if the CIGAR op isn't `S`.
fn soft_clip_len(op: Option<&rust_htslib::bam::record::Cigar>) -> i64 {
    match op {
        Some(c) if c.char() == 'S' => i64::from(c.len()),
        _ => 0,
    }
}

/// Returns `(start, pos)` for a read.
///
/// - `start`: leftmost aligned position (used for buffer-flush decisions)
/// - `pos`: 5′ coordinate accounting for soft-clipping (used for grouping)
///
/// Matches Python `get_read_position()` with default `soft_clip_threshold=4`.
fn get_read_position(record: &Record) -> (i64, i64) {
    let cigar = record.cigar();
    if record.is_reverse() {
        let start = record.pos();
        let pos = cigar.end_pos() + soft_clip_len(cigar.last());
        (start, pos)
    } else {
        let pos = record.pos() - soft_clip_len(cigar.first());
        (pos, pos)
    }
}

/// Sub-key within a position group: `(is_reverse, is_spliced, tlen, read_length)`.
/// With default options, this collapses to `(is_reverse, false, 0, 0)`.
type GroupKey = (bool, bool, i64, usize);

/// Holds per-UMI read selection state: best record + reservoir-sampling counter.
struct UmiSlot {
    record: Record,
    mapq: u8,
    tie_count: u32,
}

/// Buffered read collector that mirrors Python `umi_tools`' `reads_dict`.
///
/// Structure: `pos → key → umi → UmiSlot`
///
/// `pos` is the 5′ coordinate; `key` is `(is_reverse, …)`.
/// When flushing, positions are emitted in sorted order and keys within
/// each position are emitted in sorted order (matching Python's
/// `sorted(reads_dict[p].keys())`).
struct ReadBuffer<R: TieBreakRng> {
    groups: BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, UmiSlot>>>,
    rng: R,
}

impl<R: TieBreakRng> ReadBuffer<R> {
    const fn new(rng: R) -> Self {
        Self {
            groups: BTreeMap::new(),
            rng,
        }
    }

    /// Add a record to the buffer, performing reservoir-sampling read selection.
    fn add(&mut self, record: Record, pos: i64, key: GroupKey, umi: Vec<u8>) {
        let umi_map = self.groups.entry(pos).or_default().entry(key).or_default();

        let Some(slot) = umi_map.get_mut(&umi) else {
            let mapq = record.mapq();
            umi_map.insert(
                umi,
                UmiSlot {
                    record,
                    mapq,
                    tie_count: 0,
                },
            );
            return;
        };

        let record_mapq = record.mapq();
        match slot.mapq.cmp(&record_mapq) {
            std::cmp::Ordering::Greater => {}
            std::cmp::Ordering::Less => {
                slot.record = record;
                slot.mapq = record_mapq;
                slot.tie_count = 0;
            }
            std::cmp::Ordering::Equal => {
                slot.tie_count += 1;
                if self.rng.random() < 1.0 / f64::from(slot.tie_count) {
                    slot.record = record;
                }
            }
        }
    }

    /// Drain all position groups with `pos <= threshold`, yielding records in
    /// position-sorted, then key-sorted order.
    fn drain_up_to(&mut self, threshold: i64) -> Vec<Record> {
        let rest = self.groups.split_off(&(threshold + 1));
        let drained = std::mem::replace(&mut self.groups, rest);
        drained
            .into_values()
            .flat_map(BTreeMap::into_values)
            .flat_map(HashMap::into_values)
            .map(|slot| slot.record)
            .collect()
    }

    /// Drain all remaining position groups.
    fn drain_all(&mut self) -> Vec<Record> {
        std::mem::take(&mut self.groups)
            .into_values()
            .flat_map(BTreeMap::into_values)
            .flat_map(HashMap::into_values)
            .map(|slot| slot.record)
            .collect()
    }
}

/// # Errors
///
/// Returns `DedupError` on BAM I/O failures or unknown chromosome filter.
pub fn run_dedup(
    config: &DedupConfig,
    input_path: &str,
    output: &mut dyn io::Write,
) -> Result<DedupStats, DedupError> {
    let mut reader =
        bam::Reader::from_path(input_path).map_err(|e| DedupError::BamOpen(e.to_string()))?;
    let header = bam::Header::from_template(reader.header());

    let format = if config.out_sam {
        bam::Format::Sam
    } else {
        bam::Format::Bam
    };

    let mut writer = bam::Writer::from_stdout(&header, format)
        .map_err(|e| DedupError::BamWrite(e.to_string()))?;

    // Optional chromosome filter
    let chrom_filter: Option<i32> = config
        .chrom
        .as_ref()
        .map(|c| {
            let tid = reader
                .header()
                .tid(c.as_bytes())
                .ok_or_else(|| DedupError::UnknownChrom(c.clone()))?;
            #[allow(clippy::cast_possible_wrap)]
            Ok(tid as i32)
        })
        .transpose()?;

    #[allow(clippy::cast_possible_truncation)]
    let rng = PythonRandom::new(config.random_seed as u32);
    let mut buffer = ReadBuffer::new(rng);
    let mut stats = DedupStats {
        input_reads: 0,
        output_reads: 0,
        positions: 0,
    };

    // Collect all selected records, then sort by coordinate before writing.
    // Matches Python umi_tools which calls `pysam.sort()` after processing.
    let mut output_records: Vec<Record> = Vec::new();

    let mut last_start: i64 = 0;
    let mut last_chrom: i32 = -1;

    for result in reader.records() {
        let record = result.map_err(|e| DedupError::BamRead(e.to_string()))?;

        if record.is_unmapped() {
            continue;
        }

        let tid = record.tid();

        // Chromosome filter
        if chrom_filter.is_some_and(|filter_tid| tid != filter_tid) {
            continue;
        }

        stats.input_reads += 1;

        let (start, pos) = get_read_position(&record);

        // Flush buffer when moving far enough or changing chromosome.
        // Matches Python: `if start > last_pos + 1000 or current_chr != last_chr`
        if tid != last_chrom {
            output_records.extend(buffer.drain_all());
        } else if start > last_start + 1000 {
            let threshold = start - 1000;
            output_records.extend(buffer.drain_up_to(threshold));
        }

        last_start = start;
        last_chrom = tid;

        let key: GroupKey = (record.is_reverse(), false, 0, 0);

        let umi = if config.ignore_umi {
            Vec::new()
        } else {
            extract_umi_from_name(&record, config.umi_separator)
        };

        buffer.add(record, pos, key, umi);
    }

    output_records.extend(buffer.drain_all());

    // Sort by coordinate (tid, pos) to match `pysam.sort()` / `samtools sort`.
    output_records.sort_by(|a, b| a.tid().cmp(&b.tid()).then_with(|| a.pos().cmp(&b.pos())));

    stats.output_reads = output_records.len() as u64;
    for r in &output_records {
        writer
            .write(r)
            .map_err(|e| DedupError::BamWrite(e.to_string()))?;
    }

    // Drop writer to flush SAM/BAM output.
    // The output arg is unused for now (Writer writes to stdout directly).
    let _ = output;
    drop(writer);

    Ok(stats)
}

fn extract_umi_from_name(record: &Record, separator: u8) -> Vec<u8> {
    let name = record.qname();
    name.iter()
        .rposition(|&b| b == separator)
        .map_or_else(|| name.to_vec(), |pos| name[pos + 1..].to_vec())
}

#[derive(Debug, thiserror::Error)]
pub enum DedupError {
    #[error("failed to open BAM: {0}")]
    BamOpen(String),
    #[error("failed to read BAM record: {0}")]
    BamRead(String),
    #[error("failed to write BAM/SAM: {0}")]
    BamWrite(String),
    #[error("unknown chromosome: {0}")]
    UnknownChrom(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_random_matches() {
        let mut rng = PythonRandom::new(123_456_789);
        let expected = [
            0.641_400_616_185_872_6,
            0.542_189_268_096_949_5,
            0.993_175_066_283_272_1,
            0.843_252_136_686_916_6,
            0.811_733_928_337_940_6,
            0.397_173_710_078_000_4,
            0.937_095_107_912_042_5,
            0.689_102_653_165_816_2,
            0.397_110_488_525_983_74,
            0.351_025_192_423_044_75,
        ];
        for &exp in &expected {
            let got = rng.random();
            assert!(
                (got - exp).abs() < 1e-15,
                "mismatch: got {got:.20}, expected {exp:.20}"
            );
        }
    }
}
