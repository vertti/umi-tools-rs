use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::Write as IoWrite;

use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{Read as BamRead, Record};

use crate::alignment_io::{self, AlignmentFormat, AlignmentOutput};
use crate::barcode::{Barcode, BarcodeError, BarcodeExtractor};
use crate::clustering::{cluster_umis, hamming_distance, median};
use crate::gene::{Flush, GeneAssigner, GeneError, GeneOptions};
use crate::pairing::{PairingError, PairingOptions};

/// Trait for RNG used in reservoir-sampling tie-breaks.
///
/// Currently implemented by `PythonRandom` (MT19937 matching `CPython`) to get
/// identical output for compat tests. Can be swapped for any fast RNG once
/// exact-match testing is no longer needed.
pub(crate) trait TieBreakRng {
    /// Return a float in `[0, 1)`.
    fn random(&mut self) -> f64;
}

/// Mersenne Twister 19937 PRNG, matching `CPython`'s random module exactly.
///
/// Python `umi_tools` uses seeded `random.random()` for reservoir-sampling
/// tie-breaks in read selection. We replicate the identical float sequence.
pub(crate) struct PythonRandom {
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
    pub(crate) fn new(seed: u32) -> Self {
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

/// MT19937 PRNG matching `NumPy`'s `np.random.seed(int)` + `np.random.random()`.
///
/// `NumPy` seeds with `init_genrand(seed)` directly (unlike `CPython` which uses
/// `init_by_array`). Output generation (`genrand_res53`) is identical.
struct NumpyRandom {
    mt: [u32; 624],
    index: usize,
}

impl NumpyRandom {
    const N: usize = 624;

    fn new(seed: u32) -> Self {
        PythonRandom::init_genrand(seed).into()
    }

    fn random(&mut self) -> f64 {
        let a = self.next_u32() >> 5;
        let b = self.next_u32() >> 6;
        (f64::from(a) * 67_108_864.0 + f64::from(b)) * (1.0 / 9_007_199_254_740_992.0)
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

    fn generate(&mut self) {
        static MAG01: [u32; 2] = [0, PythonRandom::MATRIX_A];
        for kk in 0..PythonRandom::N - PythonRandom::M {
            let y = (self.mt[kk] & PythonRandom::UPPER_MASK)
                | (self.mt[kk + 1] & PythonRandom::LOWER_MASK);
            self.mt[kk] = self.mt[kk + PythonRandom::M] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
        }
        for kk in PythonRandom::N - PythonRandom::M..PythonRandom::N - 1 {
            let y = (self.mt[kk] & PythonRandom::UPPER_MASK)
                | (self.mt[kk + 1] & PythonRandom::LOWER_MASK);
            self.mt[kk] = self.mt[kk + PythonRandom::M - PythonRandom::N]
                ^ (y >> 1)
                ^ MAG01[(y & 1) as usize];
        }
        let y = (self.mt[PythonRandom::N - 1] & PythonRandom::UPPER_MASK)
            | (self.mt[0] & PythonRandom::LOWER_MASK);
        self.mt[PythonRandom::N - 1] =
            self.mt[PythonRandom::M - 1] ^ (y >> 1) ^ MAG01[(y & 1) as usize];
        self.index = 0;
    }
}

impl From<PythonRandom> for NumpyRandom {
    fn from(pr: PythonRandom) -> Self {
        Self {
            mt: pr.mt,
            index: pr.index,
        }
    }
}

pub use crate::clustering::DedupMethod;

#[allow(clippy::struct_excessive_bools)]
pub struct DedupConfig {
    pub method: DedupMethod,
    pub ignore_umi: bool,
    pub random_seed: u64,
    pub output_path: Option<String>,
    pub output_format: AlignmentFormat,
    pub reference: Option<String>,
    pub chrom: Option<String>,
    pub edit_distance_threshold: u32,
    pub position: PositionOptions,
    pub subset: Option<f32>,
    pub mapping_quality: u8,
    pub multimapping_detection: Option<MultimappingDetection>,
    pub buffer_whole_contig: bool,
    pub barcode: BarcodeExtractor,
    pub gene: GeneOptions,
    pub output_stats: Option<String>,
    pub pairing: PairingOptions,
    pub ignore_tlen: bool,
    pub umi_whitelist: Option<HashSet<Vec<u8>>>,
}

pub struct DedupStats {
    pub input_reads: u64,
    pub output_reads: u64,
    pub positions: u64,
}

/// Length of a trailing/leading soft-clip, or 0 if the CIGAR op isn't `S`.
pub(crate) fn soft_clip_len(op: Option<&Cigar>) -> i64 {
    match op {
        Some(c) if c.char() == 'S' => i64::from(c.len()),
        _ => 0,
    }
}

/// Grouping-key options shared by dedup and group.
#[derive(Debug, Clone, Copy)]
pub struct PositionOptions {
    /// Keep spliced and unspliced reads at the same position apart.
    pub spliced_is_unique: bool,
    /// A 5′ soft clip longer than this counts as splicing.
    pub soft_clip_threshold: f64,
    /// Add the read length to the grouping key.
    pub read_length: bool,
}

impl Default for PositionOptions {
    fn default() -> Self {
        Self {
            spliced_is_unique: false,
            soft_clip_threshold: 4.0,
            read_length: false,
        }
    }
}

impl PositionOptions {
    /// The `(splice, read_length)` parts of a `GroupKey` for a read.
    pub(crate) fn key_parts(self, position: &ReadPosition, record: &Record) -> (i64, usize) {
        let splice = if self.spliced_is_unique {
            position.splice_offset
        } else {
            0
        };
        let length = if self.read_length {
            record.seq_len()
        } else {
            0
        };
        (splice, length)
    }
}

/// A read's position as `umi_tools.sam_methods.get_read_position` computes it.
pub(crate) struct ReadPosition {
    /// Leftmost aligned position, used for buffer-flush decisions.
    pub(crate) start: i64,
    /// 5′ coordinate accounting for soft-clipping, used for grouping.
    pub(crate) pos: i64,
    /// Offset of the first splice from the 5′ end, or 0 when the read does not count as spliced.
    pub(crate) splice_offset: i64,
}

/// 5′ coordinate of a read accounting for soft-clipping.
pub(crate) fn five_prime_position(record: &Record) -> i64 {
    let cigar = record.cigar();
    if record.is_reverse() {
        cigar.end_pos() + soft_clip_len(cigar.last())
    } else {
        record.pos() - soft_clip_len(cigar.first())
    }
}

pub(crate) fn get_read_position(record: &Record, soft_clip_threshold: f64) -> ReadPosition {
    let cigar = record.cigar();
    let has_splice = cigar.iter().any(|op| op.char() == 'N');
    let pos = five_prime_position(record);
    if record.is_reverse() {
        #[allow(clippy::cast_precision_loss)]
        let clipped = soft_clip_len(cigar.first()) as f64 > soft_clip_threshold;
        let splice_offset = if has_splice || clipped {
            find_splice(cigar.iter().rev())
        } else {
            0
        };
        ReadPosition {
            start: record.pos(),
            pos,
            splice_offset,
        }
    } else {
        #[allow(clippy::cast_precision_loss)]
        let clipped = soft_clip_len(cigar.last()) as f64 > soft_clip_threshold;
        let splice_offset = if has_splice || clipped {
            find_splice(cigar.iter())
        } else {
            0
        };
        ReadPosition {
            start: pos,
            pos,
            splice_offset,
        }
    }
}

/// `umi_tools.sam_methods.find_splice`: reference offset of the first `N` or `S`
/// after a skipped leading soft clip, or 0 when there is none. Python returns
/// `False` for none, which compares and hashes equal to 0 in the grouping key.
fn find_splice<'a>(ops: impl Iterator<Item = &'a Cigar>) -> i64 {
    let mut ops = ops.peekable();
    let mut offset = ops
        .next_if(|op| op.char() == 'S')
        .map_or(0, |first| i64::from(first.len()));
    for op in ops {
        match op.char() {
            'N' | 'S' => return offset,
            'M' | 'D' | '=' | 'X' => offset += i64::from(op.len()),
            _ => {}
        }
    }
    0
}

/// Aligner tag consulted to break MAPQ ties between reads with the same position and UMI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultimappingDetection {
    Nh,
    X0,
    Xt,
}

impl MultimappingDetection {
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "NH" => Some(Self::Nh),
            "X0" => Some(Self::X0),
            "XT" => Some(Self::Xt),
            _ => None,
        }
    }

    const fn tag(self) -> &'static [u8; 2] {
        match self {
            Self::Nh => b"NH",
            Self::X0 => b"X0",
            Self::Xt => b"XT",
        }
    }

    /// `Some(true)` when `candidate` maps more uniquely than `current`, `Some(false)`
    /// when it maps less uniquely, `None` when the tag does not separate them.
    fn prefers_candidate(
        self,
        current: &Record,
        candidate: &Record,
    ) -> Result<Option<bool>, DedupError> {
        let tag = self.tag();
        let missing = |record: &Record| {
            DedupError::MissingTag(
                String::from_utf8_lossy(record.qname()).into_owned(),
                String::from_utf8_lossy(tag).into_owned(),
            )
        };
        match self {
            Self::Nh | Self::X0 => {
                let old = aux_int(current, tag).ok_or_else(|| missing(current))?;
                let new = aux_int(candidate, tag).ok_or_else(|| missing(candidate))?;
                Ok(match old.cmp(&new) {
                    std::cmp::Ordering::Less => Some(false),
                    std::cmp::Ordering::Greater => Some(true),
                    std::cmp::Ordering::Equal => None,
                })
            }
            Self::Xt => {
                let old = aux_char(current, tag).ok_or_else(|| missing(current))?;
                let new = aux_char(candidate, tag).ok_or_else(|| missing(candidate))?;
                Ok(if old == b'U' {
                    Some(false)
                } else if new == b'U' {
                    Some(true)
                } else {
                    None
                })
            }
        }
    }
}

fn aux_int(record: &Record, tag: &[u8]) -> Option<i64> {
    match record.aux(tag).ok()? {
        Aux::I8(v) => Some(i64::from(v)),
        Aux::U8(v) => Some(i64::from(v)),
        Aux::I16(v) => Some(i64::from(v)),
        Aux::U16(v) => Some(i64::from(v)),
        Aux::I32(v) => Some(i64::from(v)),
        Aux::U32(v) => Some(i64::from(v)),
        _ => None,
    }
}

fn aux_char(record: &Record, tag: &[u8]) -> Option<u8> {
    match record.aux(tag).ok()? {
        Aux::Char(c) => Some(c),
        Aux::String(s) => s.bytes().next(),
        _ => None,
    }
}

/// Sub-key within a position group: `(is_reverse, splice_offset, tlen, read_length, cell)`.
/// With default options, this collapses to `(is_reverse, 0, 0, 0, [])`.
pub(crate) type GroupKey = (bool, i64, i64, usize, Vec<u8>);

/// Holds per-UMI read selection state: best record + reservoir-sampling counter.
pub(crate) struct UmiSlot {
    pub(crate) record: Record,
    pub(crate) mapq: u8,
    pub(crate) tie_count: u32,
    pub(crate) count: u32,
    /// Insertion order within the (pos, key) group — used for deterministic
    /// tiebreaking to match Python dict insertion order.
    pub(crate) insertion_order: u32,
}

/// Buffered read collector that mirrors Python `umi_tools`' `reads_dict`.
///
/// Structure: `pos → key → umi → UmiSlot`
///
/// `pos` is the 5′ coordinate; `key` is `(is_reverse, …)`.
/// When flushing, positions are emitted in sorted order and keys within
/// each position are emitted in sorted order (matching Python's
/// `sorted(reads_dict[p].keys())`).
struct ReadBuffer<K: Ord = i64> {
    groups: BTreeMap<K, BTreeMap<GroupKey, HashMap<Vec<u8>, UmiSlot>>>,
    /// Per-(pos, key) insertion counters for deterministic ordering.
    insertion_counters: BTreeMap<K, BTreeMap<GroupKey, u32>>,
}

impl ReadBuffer<i64> {
    /// Drain all position groups with `pos <= threshold`, applying UMI dedup selection.
    fn drain_up_to(
        &mut self,
        threshold: i64,
        method: DedupMethod,
        edit_threshold: u32,
        stats_ctx: &mut Option<StatsContext>,
        umi_whitelist: Option<&HashSet<Vec<u8>>>,
    ) -> Vec<Record> {
        let rest = self.groups.split_off(&(threshold + 1));
        let drained = std::mem::replace(&mut self.groups, rest);
        // Clean up insertion counters for drained positions
        let rest_counters = self.insertion_counters.split_off(&(threshold + 1));
        let _ = std::mem::replace(&mut self.insertion_counters, rest_counters);
        Self::apply_selection(drained, method, edit_threshold, stats_ctx, umi_whitelist)
    }
}

impl<K: Ord + Clone> ReadBuffer<K> {
    const fn new() -> Self {
        Self {
            groups: BTreeMap::new(),
            insertion_counters: BTreeMap::new(),
        }
    }

    /// Add a record to the buffer, performing reservoir-sampling read selection.
    fn add(
        &mut self,
        record: Record,
        pos: K,
        key: GroupKey,
        umi: Vec<u8>,
        rng: &mut impl TieBreakRng,
        detection: Option<MultimappingDetection>,
    ) -> Result<(), DedupError> {
        let umi_map = self
            .groups
            .entry(pos.clone())
            .or_default()
            .entry(key.clone())
            .or_default();

        let Some(slot) = umi_map.get_mut(&umi) else {
            let counter = self
                .insertion_counters
                .entry(pos)
                .or_default()
                .entry(key)
                .or_default();
            let order = *counter;
            *counter += 1;
            let mapq = record.mapq();
            umi_map.insert(
                umi,
                UmiSlot {
                    record,
                    mapq,
                    tie_count: 0,
                    count: 1,
                    insertion_order: order,
                },
            );
            return Ok(());
        };

        slot.count += 1;

        let record_mapq = record.mapq();
        match slot.mapq.cmp(&record_mapq) {
            std::cmp::Ordering::Greater => {}
            std::cmp::Ordering::Less => {
                slot.record = record;
                slot.mapq = record_mapq;
                slot.tie_count = 0;
            }
            std::cmp::Ordering::Equal => {
                let verdict = detection
                    .map(|d| d.prefers_candidate(&slot.record, &record))
                    .transpose()?
                    .flatten();
                match verdict {
                    Some(false) => {}
                    Some(true) => {
                        // umi_tools resets the tie counter and then still draws once with probability 1.
                        slot.record = record;
                        slot.tie_count = 1;
                        rng.random();
                    }
                    None => {
                        slot.tie_count += 1;
                        if rng.random() < 1.0 / f64::from(slot.tie_count) {
                            slot.record = record;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Drain one group, applying UMI dedup selection.
    fn drain_key(
        &mut self,
        key: &K,
        method: DedupMethod,
        edit_threshold: u32,
        stats_ctx: &mut Option<StatsContext>,
        umi_whitelist: Option<&HashSet<Vec<u8>>>,
    ) -> Vec<Record> {
        let Some(key_map) = self.groups.remove(key) else {
            return Vec::new();
        };
        self.insertion_counters.remove(key);
        let drained = BTreeMap::from([(key.clone(), key_map)]);
        Self::apply_selection(drained, method, edit_threshold, stats_ctx, umi_whitelist)
    }

    /// Drain all remaining position groups, applying UMI dedup selection.
    fn drain_all(
        &mut self,
        method: DedupMethod,
        edit_threshold: u32,
        stats_ctx: &mut Option<StatsContext>,
        umi_whitelist: Option<&HashSet<Vec<u8>>>,
    ) -> Vec<Record> {
        let drained = std::mem::take(&mut self.groups);
        self.insertion_counters.clear();
        Self::apply_selection(drained, method, edit_threshold, stats_ctx, umi_whitelist)
    }

    /// Apply method-specific UMI selection to drained position groups.
    fn apply_selection(
        groups: BTreeMap<K, BTreeMap<GroupKey, HashMap<Vec<u8>, UmiSlot>>>,
        method: DedupMethod,
        edit_threshold: u32,
        stats_ctx: &mut Option<StatsContext>,
        umi_whitelist: Option<&HashSet<Vec<u8>>>,
    ) -> Vec<Record> {
        let mut records = Vec::new();
        for key_map in groups.into_values() {
            for umi_map in key_map.into_values() {
                if stats_ctx.is_some() {
                    let selected_with_counts =
                        select_umis_with_cluster_counts(method, &umi_map, edit_threshold);
                    let mut bundle_records: Vec<&Record> = Vec::new();
                    let mut selected_umis = Vec::new();
                    let mut cluster_counts = Vec::new();
                    for (umi, cluster_count) in &selected_with_counts {
                        if umi_whitelist.is_some_and(|wl| !wl.contains(umi)) {
                            continue;
                        }
                        if let Some(slot) = umi_map.get(umi) {
                            bundle_records.push(&slot.record);
                            selected_umis.push(umi.clone());
                            cluster_counts.push(*cluster_count);
                        }
                    }
                    if let Some(ctx) = stats_ctx.as_mut() {
                        ctx.collector.record_bundle(
                            &umi_map,
                            &selected_umis,
                            &cluster_counts,
                            &bundle_records,
                            &ctx.barcode,
                            &mut ctx.read_gen,
                        );
                    }
                    for r in bundle_records {
                        records.push(r.clone());
                    }
                } else {
                    let selected = select_umis(method, &umi_map, edit_threshold);
                    for umi in &selected {
                        if umi_whitelist.is_some_and(|wl| !wl.contains(umi)) {
                            continue;
                        }
                        if let Some(slot) = umi_map.get(umi) {
                            records.push(slot.record.clone());
                        }
                    }
                }
            }
        }
        records
    }
}

/// Bundles the stats collector + read generator for passing through drain calls.
struct StatsContext {
    collector: StatsCollector,
    read_gen: RandomReadGenerator,
    barcode: BarcodeExtractor,
}

/// Select the representative UMI of each group.
pub(crate) fn select_umis(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, UmiSlot>,
    edit_threshold: u32,
) -> Vec<Vec<u8>> {
    slot_groups(method, umi_map, edit_threshold)
        .into_iter()
        .map(|group| group[0].to_vec())
        .collect()
}

fn slot_groups(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, UmiSlot>,
    edit_threshold: u32,
) -> Vec<Vec<&[u8]>> {
    cluster_umis(
        method,
        umi_map
            .iter()
            .map(|(umi, slot)| (umi.as_slice(), slot.count, slot.insertion_order)),
        edit_threshold,
    )
}

/// Count deduplicated UMI groups from raw count/order maps.
///
/// # Panics
/// Panics if an input UMI has no corresponding insertion order.
#[allow(clippy::implicit_hasher)]
#[must_use]
pub fn count_umis(
    method: DedupMethod,
    counts: &HashMap<Vec<u8>, u32>,
    orders: &HashMap<Vec<u8>, u32>,
    edit_threshold: u32,
) -> usize {
    if method == DedupMethod::Unique {
        return counts.len();
    }
    cluster_umis(
        method,
        counts
            .iter()
            .map(|(umi, &count)| (umi.as_slice(), count, orders[umi])),
        edit_threshold,
    )
    .len()
}

fn select_umis_with_cluster_counts(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, UmiSlot>,
    edit_threshold: u32,
) -> Vec<(Vec<u8>, u32)> {
    slot_groups(method, umi_map, edit_threshold)
        .into_iter()
        .map(|group| {
            (
                group[0].to_vec(),
                group.iter().map(|umi| umi_map[*umi].count).sum(),
            )
        })
        .collect()
}

/// Mean pairwise Hamming distance between UMIs. Returns -1.0 for single UMI.
#[allow(clippy::cast_precision_loss)]
fn get_average_umi_distance(umis: &[&[u8]]) -> f64 {
    if umis.len() <= 1 {
        return -1.0;
    }
    let mut total: u64 = 0;
    let mut count: u64 = 0;
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            total += u64::from(hamming_distance(umis[i], umis[j]));
            count += 1;
        }
    }
    total as f64 / count as f64
}

/// Pre-scans BAM to build UMI frequency distribution for null model sampling.
struct RandomReadGenerator {
    keys: Vec<Vec<u8>>,
    cdf: Vec<f64>,
    rng: NumpyRandom,
    random_umis: Vec<Vec<u8>>,
    random_ix: usize,
    fill_size: usize,
}

impl RandomReadGenerator {
    fn new(
        bam_path: &str,
        reference: Option<&str>,
        barcode: &BarcodeExtractor,
        chrom: Option<&str>,
        seed: u32,
    ) -> Result<Self, DedupError> {
        let mut reader = alignment_io::open_reader(bam_path, reference)
            .map_err(|e| DedupError::BamOpen(e.to_string()))?;

        let chrom_tid: Option<i32> = chrom
            .map(|c| {
                let tid = reader
                    .header()
                    .tid(c.as_bytes())
                    .ok_or_else(|| DedupError::UnknownChrom(c.to_string()))?;
                #[allow(clippy::cast_possible_wrap)]
                Ok::<i32, DedupError>(tid as i32)
            })
            .transpose()?;

        // Count UMI frequencies, preserving insertion order (order of first appearance).
        let mut umi_order: Vec<Vec<u8>> = Vec::new();
        let mut umi_counts: HashMap<Vec<u8>, u64> = HashMap::new();

        for result in reader.records() {
            let record = result.map_err(|e| DedupError::BamRead(e.to_string()))?;
            if record.is_unmapped() {
                continue;
            }
            if record.is_last_in_template() {
                continue;
            }
            if let Some(filter_tid) = chrom_tid
                && record.tid() != filter_tid
            {
                continue;
            }
            let umi = match barcode.extract(&record) {
                Ok(barcode) => barcode.umi,
                Err(BarcodeError::MissingTag(_)) => continue,
                Err(e) => return Err(e.into()),
            };
            let entry = umi_counts.entry(umi.clone());
            if matches!(entry, std::collections::hash_map::Entry::Vacant(_)) {
                umi_order.push(umi);
            }
            *entry.or_insert(0) += 1;
        }

        // Build CDF from frequencies in insertion order.
        #[allow(clippy::cast_precision_loss)]
        let total: f64 = umi_counts.values().sum::<u64>() as f64;
        let mut cdf = Vec::with_capacity(umi_order.len());
        let mut cumsum = 0.0;
        for key in &umi_order {
            #[allow(clippy::cast_precision_loss)]
            {
                cumsum += umi_counts[key] as f64 / total;
            }
            cdf.push(cumsum);
        }

        let mut rng = Self {
            keys: umi_order,
            cdf,
            rng: NumpyRandom::new(seed),
            random_umis: Vec::new(),
            random_ix: 0,
            fill_size: 100_000,
        };
        rng.refill();
        Ok(rng)
    }

    fn refill(&mut self) {
        self.random_umis.clear();
        if self.keys.is_empty() {
            return;
        }
        self.random_umis.reserve(self.fill_size);
        for _ in 0..self.fill_size {
            let r = self.rng.random();
            let idx = self
                .cdf
                .partition_point(|&c| c <= r)
                .min(self.keys.len() - 1);
            self.random_umis.push(self.keys[idx].clone());
        }
        self.random_ix = 0;
    }

    fn get_umis(&mut self, n: usize) -> Vec<Vec<u8>> {
        if self.keys.is_empty() {
            return Vec::new();
        }
        if n >= self.fill_size - self.random_ix {
            if n > self.fill_size {
                self.fill_size = n * 2;
            }
            self.refill();
        }
        let result = self.random_umis[self.random_ix..self.random_ix + n].to_vec();
        self.random_ix += n;
        result
    }
}

/// Accumulates per-bundle stats during dedup for the 3 stats output files.
struct StatsCollector {
    // Per-UMI-per-position: (umi, count) tuples
    pre_umi_counts: Vec<(Vec<u8>, u32)>,
    post_umi_counts: Vec<(Vec<u8>, u32)>,
    // Edit distance stats per bundle
    pre_cluster_stats: Vec<f64>,
    post_cluster_stats: Vec<f64>,
    pre_cluster_stats_null: Vec<f64>,
    post_cluster_stats_null: Vec<f64>,
}

impl StatsCollector {
    const fn new() -> Self {
        Self {
            pre_umi_counts: Vec::new(),
            post_umi_counts: Vec::new(),
            pre_cluster_stats: Vec::new(),
            post_cluster_stats: Vec::new(),
            pre_cluster_stats_null: Vec::new(),
            post_cluster_stats_null: Vec::new(),
        }
    }

    fn record_bundle(
        &mut self,
        umi_map: &HashMap<Vec<u8>, UmiSlot>,
        selected_umis: &[Vec<u8>],
        cluster_counts: &[u32],
        selected_records: &[&Record],
        barcode: &BarcodeExtractor,
        read_gen: &mut RandomReadGenerator,
    ) {
        // Pre-dedup: all UMIs in the bundle
        let pre_umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
        for (umi, slot) in umi_map {
            self.pre_umi_counts.push((umi.clone(), slot.count));
        }
        let avg_dist = get_average_umi_distance(&pre_umis);
        self.pre_cluster_stats.push(avg_dist);

        let cluster_size = pre_umis.len();
        let random_umis = read_gen.get_umis(cluster_size);
        let random_refs: Vec<&[u8]> = random_umis.iter().map(Vec::as_slice).collect();
        let avg_null = get_average_umi_distance(&random_refs);
        self.pre_cluster_stats_null.push(avg_null);

        // Post-dedup: selected UMIs with cluster-aggregated counts
        for (umi, &count) in selected_umis.iter().zip(cluster_counts) {
            self.post_umi_counts.push((umi.clone(), count));
        }

        // Post-dedup edit distance from the actual output records' UMIs
        let post_umis: Vec<Vec<u8>> = selected_records
            .iter()
            .map(|r| barcode.extract(r).map(|b| b.umi).unwrap_or_default())
            .collect();
        let post_refs: Vec<&[u8]> = post_umis.iter().map(Vec::as_slice).collect();
        let avg_post = get_average_umi_distance(&post_refs);
        self.post_cluster_stats.push(avg_post);

        let post_size = post_umis.len();
        let random_umis_post = read_gen.get_umis(post_size);
        let random_post_refs: Vec<&[u8]> = random_umis_post.iter().map(Vec::as_slice).collect();
        let avg_null_post = get_average_umi_distance(&random_post_refs);
        self.post_cluster_stats_null.push(avg_null_post);
    }

    fn write_files(&self, prefix: &str, method_name: &str) -> Result<(), DedupError> {
        self.write_per_umi_per_position(prefix)?;
        self.write_per_umi(prefix)?;
        self.write_edit_distance(prefix, method_name)?;
        Ok(())
    }

    fn write_per_umi_per_position(&self, prefix: &str) -> Result<(), DedupError> {
        let mut pre_counts: HashMap<u32, u32> = HashMap::new();
        let mut post_counts: HashMap<u32, u32> = HashMap::new();
        for (_, count) in &self.pre_umi_counts {
            *pre_counts.entry(*count).or_default() += 1;
        }
        for (_, count) in &self.post_umi_counts {
            *post_counts.entry(*count).or_default() += 1;
        }

        let mut all_counts: Vec<u32> = pre_counts
            .keys()
            .chain(post_counts.keys())
            .copied()
            .collect::<HashSet<u32>>()
            .into_iter()
            .collect();
        all_counts.sort_unstable();

        let path = format!("{prefix}_per_umi_per_position.tsv");
        let mut f =
            File::create(&path).map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        writeln!(f, "counts\tinstances_pre\tinstances_post")
            .map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        for count in &all_counts {
            let pre = pre_counts.get(count).unwrap_or(&0);
            let post = post_counts.get(count).unwrap_or(&0);
            writeln!(f, "{count}\t{pre}\t{post}")
                .map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        }
        Ok(())
    }

    fn write_per_umi(&self, prefix: &str) -> Result<(), DedupError> {
        // Aggregate per UMI: median_counts, times_observed, total_counts
        let pre_agg = Self::aggregate_per_umi(&self.pre_umi_counts);
        let post_agg = Self::aggregate_per_umi(&self.post_umi_counts);

        // Sorted union of UMI keys
        let mut all_umis: Vec<Vec<u8>> = pre_agg
            .keys()
            .chain(post_agg.keys())
            .cloned()
            .collect::<HashSet<Vec<u8>>>()
            .into_iter()
            .collect();
        all_umis.sort();

        let path = format!("{prefix}_per_umi.tsv");
        let mut f =
            File::create(&path).map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        writeln!(
            f,
            "UMI\tmedian_counts_pre\ttimes_observed_pre\ttotal_counts_pre\t\
             median_counts_post\ttimes_observed_post\ttotal_counts_post"
        )
        .map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;

        for umi in &all_umis {
            let (med_pre, obs_pre, tot_pre) = pre_agg.get(umi).unwrap_or(&(0, 0, 0));
            let (med_post, obs_post, tot_post) = post_agg.get(umi).unwrap_or(&(0, 0, 0));
            let umi_str = std::str::from_utf8(umi).unwrap_or("?");
            writeln!(
                f,
                "{umi_str}\t{med_pre}\t{obs_pre}\t{tot_pre}\t{med_post}\t{obs_post}\t{tot_post}"
            )
            .map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        }
        Ok(())
    }

    /// Returns map: umi → (`median_counts`, `times_observed`, `total_counts`)
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    fn aggregate_per_umi(umi_counts: &[(Vec<u8>, u32)]) -> HashMap<Vec<u8>, (i64, i64, i64)> {
        let mut grouped: HashMap<Vec<u8>, Vec<u32>> = HashMap::new();
        for (umi, count) in umi_counts {
            grouped.entry(umi.clone()).or_default().push(*count);
        }
        grouped
            .into_iter()
            .map(|(umi, counts)| {
                let times_observed = counts.len() as i64;
                let total: i64 = counts.iter().map(|&c| i64::from(c)).sum();
                let med = median(&counts);
                // Python: .fillna(0).astype(int) truncates toward zero (same as floor for positive)
                let median_int = med as i64;
                (umi, (median_int, times_observed, total))
            })
            .collect()
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn write_edit_distance(&self, prefix: &str, method_name: &str) -> Result<(), DedupError> {
        // Find max edit distance across all stats
        let all_stats = self
            .pre_cluster_stats
            .iter()
            .chain(&self.post_cluster_stats)
            .chain(&self.pre_cluster_stats_null)
            .chain(&self.post_cluster_stats_null);
        let max_ed = all_stats.copied().fold(0.0_f64, f64::max) as i32;

        // bins = range(-1, max_ed + 2)  →  [-1, 0, 1, ..., max_ed+1]
        let bins: Vec<i32> = (-1..=max_ed + 1).collect();
        let nbins = bins.len();

        let digitize = |values: &[f64]| -> Vec<usize> {
            // np.digitize(values, bins, right=True): returns i such that
            // bins[i-1] < v <= bins[i]. Equivalent to searchsorted(side='left').
            values
                .iter()
                .map(|&v| bins.partition_point(|&b| f64::from(b) < v).min(nbins))
                .collect()
        };

        let bincount = |binned: &[usize], minlength: usize| -> Vec<u64> {
            let mut counts = vec![0u64; minlength];
            for &b in binned {
                if b < counts.len() {
                    counts[b] += 1;
                }
            }
            counts
        };

        let minlength = (max_ed + 3) as usize;

        let pre_binned = digitize(&self.pre_cluster_stats);
        let post_binned = digitize(&self.post_cluster_stats);
        let pre_null_binned = digitize(&self.pre_cluster_stats_null);
        let post_null_binned = digitize(&self.post_cluster_stats_null);

        let pre_counts = bincount(&pre_binned, minlength);
        let post_counts = bincount(&post_binned, minlength);
        let pre_null_counts = bincount(&pre_null_binned, minlength);
        let post_null_counts = bincount(&post_null_binned, minlength);

        let path = format!("{prefix}_edit_distance.tsv");
        let mut f =
            File::create(&path).map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        writeln!(
            f,
            "unique\tunique_null\t{method_name}\t{method_name}_null\tedit_distance"
        )
        .map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;

        for i in 0..minlength {
            let ed_label = if i == 0 {
                "Single_UMI".to_string()
            } else if i < bins.len() {
                bins[i].to_string()
            } else {
                (i - 1).to_string()
            };
            let pre = pre_counts.get(i).unwrap_or(&0);
            let post = post_counts.get(i).unwrap_or(&0);
            let pre_null = pre_null_counts.get(i).unwrap_or(&0);
            let post_null = post_null_counts.get(i).unwrap_or(&0);
            writeln!(f, "{pre}\t{pre_null}\t{post}\t{post_null}\t{ed_label}")
                .map_err(|e| DedupError::StatsWrite(path.clone(), e.to_string()))?;
        }
        Ok(())
    }
}

/// # Errors
///
/// Returns `DedupError` on BAM I/O failures or unknown chromosome filter.
#[allow(clippy::too_many_lines)]
pub fn run_dedup(config: &DedupConfig, input_path: &str) -> Result<DedupStats, DedupError> {
    config.pairing.validate(false)?;
    let mut source = alignment_io::RecordSource::whole(input_path, config.reference.as_deref())
        .map_err(|e| DedupError::BamOpen(e.to_string()))?;
    let assigner = GeneAssigner::new(&config.gene, source.header())?;
    if let Some(map) = assigner.as_ref().and_then(GeneAssigner::transcript_map) {
        source = alignment_io::RecordSource::by_contig(
            input_path,
            config.reference.as_deref(),
            map.genes.clone(),
        )
        .map_err(|e| DedupError::BamOpen(e.to_string()))?;
    }
    let mut flusher = assigner.as_ref().map(GeneAssigner::flusher);
    let header = alignment_io::coordinate_sorted_header(source.header());

    let mut writer = alignment_io::open_writer(
        &header,
        AlignmentOutput {
            path: config.output_path.as_deref(),
            format: config.output_format,
            reference: config.reference.as_deref(),
        },
    )
    .map_err(|e| DedupError::BamWrite(e.to_string()))?;

    // Optional chromosome filter
    let chrom_filter: Option<i32> = config
        .chrom
        .as_ref()
        .map(|c| {
            let tid = source
                .header()
                .tid(c.as_bytes())
                .ok_or_else(|| DedupError::UnknownChrom(c.clone()))?;
            #[allow(clippy::cast_possible_wrap)]
            Ok::<i32, DedupError>(tid as i32)
        })
        .transpose()?;

    #[allow(clippy::cast_possible_truncation)]
    let mut rng = PythonRandom::new(config.random_seed as u32);
    let mut buffer = ReadBuffer::<i64>::new();
    let mut gene_buffer = ReadBuffer::<Vec<u8>>::new();
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

    // Stats collection (optional, only when --output-stats is set)
    #[allow(clippy::cast_possible_truncation)]
    let mut stats_ctx: Option<StatsContext> = config
        .output_stats
        .as_ref()
        .map(|_| {
            let read_gen = RandomReadGenerator::new(
                input_path,
                config.reference.as_deref(),
                &config.barcode,
                config.chrom.as_deref(),
                config.random_seed as u32,
            )?;
            Ok::<_, DedupError>(StatsContext {
                collector: StatsCollector::new(),
                read_gen,
                barcode: config.barcode.clone(),
            })
        })
        .transpose()?;

    let wl_ref = config.umi_whitelist.as_ref();

    while let Some(record) = source
        .read_next()
        .map_err(|e| DedupError::BamRead(e.to_string()))?
    {
        let tid = record.tid();
        if chrom_filter.is_some_and(|filter_tid| tid != filter_tid) {
            continue;
        }

        let triage = config.pairing.triage(&record, false);
        if triage.is_read2 {
            continue;
        }
        stats.input_reads += 1;
        if !triage.grouped {
            continue;
        }

        // Subset check consumes one RNG call per mapped read (before buffer.add)
        if config.subset.is_some_and(|s| rng.random() >= f64::from(s)) {
            continue;
        }

        if record.mapq() < config.mapping_quality {
            continue;
        }

        let Some(Barcode { umi, cell }) =
            config.barcode.for_grouping(&record, config.ignore_umi)?
        else {
            continue;
        };

        if let (Some(assigner), Some(flusher)) = (&assigner, flusher.as_mut()) {
            let Some(gene) = assigner.gene(&record) else {
                continue;
            };
            match flusher.before_read(tid) {
                Flush::None => {}
                Flush::All => output_records.extend(gene_buffer.drain_all(
                    config.method,
                    config.edit_distance_threshold,
                    &mut stats_ctx,
                    wl_ref,
                )),
                Flush::Gene(done) => output_records.extend(gene_buffer.drain_key(
                    &done,
                    config.method,
                    config.edit_distance_threshold,
                    &mut stats_ctx,
                    wl_ref,
                )),
            }
            flusher.after_read(tid, &gene);
            gene_buffer.add(
                record,
                gene,
                (false, 0, 0, 0, cell),
                umi,
                &mut rng,
                config.multimapping_detection,
            )?;
        } else {
            let position = get_read_position(&record, config.position.soft_clip_threshold);
            let start = position.start;

            // Flush buffer when moving far enough or changing chromosome.
            if tid != last_chrom {
                output_records.extend(buffer.drain_all(
                    config.method,
                    config.edit_distance_threshold,
                    &mut stats_ctx,
                    wl_ref,
                ));
            } else if !config.buffer_whole_contig && start > last_start + 1000 {
                let threshold = start - 1000;
                output_records.extend(buffer.drain_up_to(
                    threshold,
                    config.method,
                    config.edit_distance_threshold,
                    &mut stats_ctx,
                    wl_ref,
                ));
            }

            last_start = start;
            last_chrom = tid;

            let tlen = if config.pairing.paired && !config.ignore_tlen {
                record.insert_size()
            } else {
                0
            };
            let (splice, length) = config.position.key_parts(&position, &record);
            let key: GroupKey = (record.is_reverse(), splice, tlen, length, cell);
            buffer.add(
                record,
                position.pos,
                key,
                umi,
                &mut rng,
                config.multimapping_detection,
            )?;
        }
    }

    output_records.extend(buffer.drain_all(
        config.method,
        config.edit_distance_threshold,
        &mut stats_ctx,
        wl_ref,
    ));
    output_records.extend(gene_buffer.drain_all(
        config.method,
        config.edit_distance_threshold,
        &mut stats_ctx,
        wl_ref,
    ));

    // Paired mode: second pass to find R2 mates of surviving R1 reads.
    if config.pairing.paired {
        let mut mate_set: HashSet<(Vec<u8>, i32, i64)> = HashSet::new();
        for r1 in &output_records {
            mate_set.insert((r1.qname().to_vec(), r1.mtid(), r1.mpos()));
        }
        let mut reader2 = alignment_io::open_reader(input_path, config.reference.as_deref())
            .map_err(|e| DedupError::BamOpen(e.to_string()))?;
        for result in reader2.records() {
            let record = result.map_err(|e| DedupError::BamRead(e.to_string()))?;
            if record.is_unmapped() || record.is_mate_unmapped() {
                continue;
            }
            if record.is_first_in_template() {
                continue;
            }
            let key = (record.qname().to_vec(), record.tid(), record.pos());
            if mate_set.remove(&key) {
                output_records.push(record);
            }
        }
    }

    // Sort by coordinate (tid, pos) to match `pysam.sort()` / `samtools sort`.
    output_records.sort_by_key(alignment_io::coordinate_sort_key);

    stats.output_reads = output_records.len() as u64;
    for r in &output_records {
        writer
            .write(r)
            .map_err(|e| DedupError::BamWrite(e.to_string()))?;
    }

    drop(writer);

    // Write stats files if requested
    if let (Some(prefix), Some(ctx)) = (&config.output_stats, &stats_ctx) {
        let method_name = match config.method {
            DedupMethod::Unique => "unique",
            DedupMethod::Percentile => "percentile",
            DedupMethod::Cluster => "cluster",
            DedupMethod::Adjacency => "adjacency",
            DedupMethod::Directional => "directional",
        };
        ctx.collector.write_files(prefix, method_name)?;
    }

    Ok(stats)
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
    #[error("invalid regex: {0}")]
    InvalidRegex(String),
    #[error("failed to write stats file {0}: {1}")]
    StatsWrite(String, String),
    #[error("read {0} has no {1} tag for --multimapping-detection-method")]
    MissingTag(String, String),
    #[error(transparent)]
    Barcode(#[from] BarcodeError),
    #[error(transparent)]
    Gene(#[from] GeneError),
    #[error(transparent)]
    Pairing(#[from] PairingError),
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::HeaderView;

    use super::*;

    fn record(flag: u16, cigar: &str) -> Record {
        let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:100000\n");
        let line = format!("r\t{flag}\tchr1\t101\t60\t{cigar}\t*\t0\t0\t*\t*");
        Record::from_sam(&header, line.as_bytes()).unwrap()
    }

    #[test]
    fn read_position_matches_umi_tools() {
        // Expected values from umi_tools.sam_methods.get_read_position via pysam (0-based).
        let cases: [(u16, &str, f64, i64, i64, i64); 13] = [
            (0, "10M", 4.0, 100, 100, 0),
            (0, "5S10M", 4.0, 95, 95, 0),
            (0, "10M5S", 4.0, 100, 100, 10),
            (0, "10M5S", 5.0, 100, 100, 0),
            (0, "3S10M5S", 4.0, 97, 97, 13),
            (0, "10M100N10M", 4.0, 100, 100, 10),
            (0, "2S10M100N10M", 4.0, 98, 98, 12),
            (16, "10M5S", 4.0, 100, 115, 0),
            (16, "5S10M", 4.0, 100, 110, 10),
            (16, "10M100N10M5S", 4.0, 100, 225, 15),
            (16, "10M100N10M", 4.0, 100, 220, 10),
            (0, "4M2I6M50N10M", 4.0, 100, 100, 10),
            (0, "4M2D6M50N10M", 4.0, 100, 100, 12),
        ];
        for (flag, cigar, threshold, start, pos, splice) in cases {
            let p = get_read_position(&record(flag, cigar), threshold);
            assert_eq!(
                (p.start, p.pos, p.splice_offset),
                (start, pos, splice),
                "flag={flag} cigar={cigar} threshold={threshold}"
            );
        }
    }

    struct FixedRng {
        value: f64,
        draws: u32,
    }

    impl TieBreakRng for FixedRng {
        fn random(&mut self) -> f64 {
            self.draws += 1;
            self.value
        }
    }

    fn tagged(name: &str, tags: &str) -> Record {
        let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:100000\n");
        let line = format!("{name}\t0\tchr1\t101\t60\t10M\t*\t0\t0\t*\t*\t{tags}");
        Record::from_sam(&header, line.as_bytes()).unwrap()
    }

    #[test]
    fn multimapping_tag_breaks_mapq_ties_like_umi_tools() {
        let key: GroupKey = (false, 0, 0, 0, Vec::new());
        let umi = b"ACGT".to_vec();
        let nh = Some(MultimappingDetection::Nh);
        let mut rng = FixedRng {
            value: 0.5,
            draws: 0,
        };
        let mut buffer = ReadBuffer::new();
        let selected = |buffer: &ReadBuffer| {
            let slot = &buffer.groups[&100][&key][&umi];
            (slot.record.qname().to_vec(), slot.tie_count)
        };

        buffer
            .add(
                tagged("a", "NH:i:3"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                nh,
            )
            .unwrap();
        buffer
            .add(
                tagged("b", "NH:i:1"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                nh,
            )
            .unwrap();
        assert_eq!(selected(&buffer), (b"b".to_vec(), 1), "fewer hits wins");
        assert_eq!(rng.draws, 1, "the replacement still consumes one draw");

        buffer
            .add(
                tagged("c", "NH:i:5"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                nh,
            )
            .unwrap();
        assert_eq!(selected(&buffer), (b"b".to_vec(), 1), "more hits loses");
        assert_eq!(rng.draws, 1, "losing consumes no draw");

        buffer
            .add(
                tagged("d", "NH:i:1"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                nh,
            )
            .unwrap();
        assert_eq!(selected(&buffer), (b"b".to_vec(), 2), "equal hits sample");
        assert_eq!(rng.draws, 2);

        let err = buffer
            .add(
                tagged("e", "AS:i:1"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                nh,
            )
            .unwrap_err();
        assert!(matches!(err, DedupError::MissingTag(_, _)));
    }

    #[test]
    fn xt_unique_beats_repeat() {
        let key: GroupKey = (false, 0, 0, 0, Vec::new());
        let umi = b"ACGT".to_vec();
        let xt = Some(MultimappingDetection::Xt);
        let mut rng = FixedRng {
            value: 0.5,
            draws: 0,
        };
        let mut buffer = ReadBuffer::new();
        buffer
            .add(
                tagged("r", "XT:A:R"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                xt,
            )
            .unwrap();
        buffer
            .add(
                tagged("u", "XT:A:U"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                xt,
            )
            .unwrap();
        buffer
            .add(
                tagged("r2", "XT:A:R"),
                100,
                key.clone(),
                umi.clone(),
                &mut rng,
                xt,
            )
            .unwrap();
        let slot = &buffer.groups[&100][&key][&umi];
        assert_eq!(slot.record.qname(), b"u");
        assert_eq!(rng.draws, 1);
    }

    #[test]
    fn key_parts_follow_options() {
        // read_length is the SEQ length as stored, soft clips included.
        let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:100000\n");
        let line = b"r\t0\tchr1\t101\t60\t10M5S\t*\t0\t0\tACGTACGTACGTACG\t*";
        let read = Record::from_sam(&header, line).unwrap();
        let position = get_read_position(&read, 4.0);
        assert_eq!(
            PositionOptions::default().key_parts(&position, &read),
            (0, 0)
        );
        let all = PositionOptions {
            spliced_is_unique: true,
            soft_clip_threshold: 4.0,
            read_length: true,
        };
        assert_eq!(all.key_parts(&position, &read), (10, 15));
    }

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
