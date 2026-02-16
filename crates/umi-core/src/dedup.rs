use std::collections::{BTreeMap, HashMap, HashSet};
use std::io;

use rust_htslib::bam::{self, Read as BamRead, Record};

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
    pub edit_distance_threshold: u32,
    pub subset: Option<f32>,
    pub extract_umi_method: String,
    pub umi_tag: Option<String>,
}

pub struct DedupStats {
    pub input_reads: u64,
    pub output_reads: u64,
    pub positions: u64,
}

/// Length of a trailing/leading soft-clip, or 0 if the CIGAR op isn't `S`.
pub(crate) fn soft_clip_len(op: Option<&rust_htslib::bam::record::Cigar>) -> i64 {
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
pub(crate) fn get_read_position(record: &Record) -> (i64, i64) {
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
pub(crate) type GroupKey = (bool, bool, i64, usize);

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
struct ReadBuffer {
    groups: BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, UmiSlot>>>,
    /// Per-(pos, key) insertion counters for deterministic ordering.
    insertion_counters: BTreeMap<i64, BTreeMap<GroupKey, u32>>,
}

impl ReadBuffer {
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
        pos: i64,
        key: GroupKey,
        umi: Vec<u8>,
        rng: &mut impl TieBreakRng,
    ) {
        let umi_map = self.groups.entry(pos).or_default().entry(key).or_default();

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
            return;
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
                slot.tie_count += 1;
                if rng.random() < 1.0 / f64::from(slot.tie_count) {
                    slot.record = record;
                }
            }
        }
    }

    /// Drain all position groups with `pos <= threshold`, applying UMI dedup selection.
    fn drain_up_to(
        &mut self,
        threshold: i64,
        method: DedupMethod,
        edit_threshold: u32,
    ) -> Vec<Record> {
        let rest = self.groups.split_off(&(threshold + 1));
        let drained = std::mem::replace(&mut self.groups, rest);
        // Clean up insertion counters for drained positions
        let rest_counters = self.insertion_counters.split_off(&(threshold + 1));
        let _ = std::mem::replace(&mut self.insertion_counters, rest_counters);
        Self::apply_selection(drained, method, edit_threshold)
    }

    /// Drain all remaining position groups, applying UMI dedup selection.
    fn drain_all(&mut self, method: DedupMethod, edit_threshold: u32) -> Vec<Record> {
        let drained = std::mem::take(&mut self.groups);
        self.insertion_counters.clear();
        Self::apply_selection(drained, method, edit_threshold)
    }

    /// Apply method-specific UMI selection to drained position groups.
    /// Emits records in the order returned by `select_umis`, preserving
    /// the deterministic ordering from clustering (count-desc, lex tiebreak).
    fn apply_selection(
        groups: BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, UmiSlot>>>,
        method: DedupMethod,
        edit_threshold: u32,
    ) -> Vec<Record> {
        let mut records = Vec::new();
        for key_map in groups.into_values() {
            for umi_map in key_map.into_values() {
                let selected = select_umis(method, &umi_map, edit_threshold);
                for umi in &selected {
                    if let Some(slot) = umi_map.get(umi) {
                        records.push(slot.record.clone());
                    }
                }
            }
        }
        records
    }
}

/// Hamming distance between two byte slices of equal length.
/// Returns `u32::MAX` if lengths differ (matching Python's `np.inf` return).
#[allow(clippy::cast_possible_truncation)]
pub(crate) fn edit_distance(a: &[u8], b: &[u8]) -> u32 {
    if a.len() != b.len() {
        return u32::MAX;
    }
    // UMIs are 5-12bp; count always fits u32
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count() as u32
}

/// Build undirected adjacency list (for cluster + adjacency methods).
/// Edge between A and B iff `edit_distance(A, B) <= threshold`.
pub(crate) fn build_adjacency_list(
    umis: &[&[u8]],
    threshold: u32,
) -> HashMap<Vec<u8>, Vec<Vec<u8>>> {
    let mut adj: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();
    for umi in umis {
        adj.entry(umi.to_vec()).or_default();
    }
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            if edit_distance(umis[i], umis[j]) <= threshold {
                adj.get_mut(umis[i]).unwrap().push(umis[j].to_vec());
                adj.get_mut(umis[j]).unwrap().push(umis[i].to_vec());
            }
        }
    }
    adj
}

/// Build directed adjacency list (for directional method).
/// Edge A→B iff `edit_distance(A,B) <= threshold AND counts[A] >= 2*counts[B] - 1`.
pub(crate) fn build_directional_adjacency_list(
    umis: &[&[u8]],
    counts: &HashMap<&[u8], u32>,
    threshold: u32,
) -> HashMap<Vec<u8>, Vec<Vec<u8>>> {
    let mut adj: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();
    for umi in umis {
        adj.entry(umi.to_vec()).or_default();
    }
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            if edit_distance(umis[i], umis[j]) <= threshold {
                let ca = counts[umis[i]];
                let cb = counts[umis[j]];
                if ca >= (2 * cb).saturating_sub(1) {
                    adj.get_mut(umis[i]).unwrap().push(umis[j].to_vec());
                }
                if cb >= (2 * ca).saturating_sub(1) {
                    adj.get_mut(umis[j]).unwrap().push(umis[i].to_vec());
                }
            }
        }
    }
    adj
}

/// BFS from `start`, following edges in `adj_list`. Returns the connected component.
pub(crate) fn bfs(start: &[u8], adj_list: &HashMap<Vec<u8>, Vec<Vec<u8>>>) -> Vec<Vec<u8>> {
    let mut searched: HashSet<Vec<u8>> = HashSet::new();
    let mut queue: Vec<Vec<u8>> = Vec::new();
    searched.insert(start.to_vec());
    queue.push(start.to_vec());
    while let Some(node) = queue.pop() {
        if let Some(neighbors) = adj_list.get(&node) {
            for next_node in neighbors {
                if searched.insert(next_node.clone()) {
                    queue.push(next_node.clone());
                }
            }
        }
    }
    let mut result: Vec<Vec<u8>> = searched.into_iter().collect();
    result.sort();
    result
}

/// Find connected components by iterating UMIs in count-descending order,
/// running BFS from each unvisited node. Matches Python `_get_connected_components_adjacency`.
pub(crate) fn connected_components(
    umis: &[&[u8]],
    counts: &HashMap<&[u8], u32>,
    orders: &HashMap<&[u8], u32>,
    adj_list: &HashMap<Vec<u8>, Vec<Vec<u8>>>,
) -> Vec<Vec<Vec<u8>>> {
    // Sort UMIs by count descending, then insertion order ascending for ties
    let mut sorted_umis: Vec<&[u8]> = umis.to_vec();
    sorted_umis.sort_by(|a, b| {
        counts[b]
            .cmp(&counts[a])
            .then_with(|| orders[a].cmp(&orders[b]))
    });

    let mut found: HashSet<Vec<u8>> = HashSet::new();
    let mut components: Vec<Vec<Vec<u8>>> = Vec::new();
    for umi in &sorted_umis {
        if !found.contains(*umi) {
            let component = bfs(umi, adj_list);
            for node in &component {
                found.insert(node.clone());
            }
            components.push(component);
        }
    }
    components
}

/// Greedy min-set-cover: select fewest UMIs (by descending count) to "cover"
/// all UMIs in the cluster via adjacency. Matches Python `_get_best_min_account`.
pub(crate) fn min_set_cover(
    cluster: &[Vec<u8>],
    adj_list: &HashMap<Vec<u8>, Vec<Vec<u8>>>,
    counts: &HashMap<&[u8], u32>,
) -> Vec<Vec<u8>> {
    if cluster.len() == 1 {
        return cluster.to_vec();
    }
    let mut sorted_nodes: Vec<&Vec<u8>> = cluster.iter().collect();
    // Sort by count desc, lex asc (BFS output is lex-sorted; Python's stable sort preserves that)
    sorted_nodes.sort_by(|a, b| {
        counts[b.as_slice()]
            .cmp(&counts[a.as_slice()])
            .then_with(|| a.cmp(b))
    });
    for i in 0..sorted_nodes.len() - 1 {
        let selected: Vec<&[u8]> = sorted_nodes[..=i].iter().map(|v| v.as_slice()).collect();
        // Compute covered nodes: selected nodes + their neighbors
        let mut covered: HashSet<&[u8]> = HashSet::new();
        for s in &selected {
            covered.insert(s);
            if let Some(neighbors) = adj_list.get(*s) {
                for n in neighbors {
                    covered.insert(n.as_slice());
                }
            }
        }
        // Check if all cluster nodes are covered
        let remaining: usize = cluster
            .iter()
            .filter(|n| !covered.contains(n.as_slice()))
            .count();
        if remaining == 0 {
            return selected.into_iter().map(<[u8]>::to_vec).collect();
        }
    }
    // Fallback: all nodes (shouldn't reach here for valid inputs)
    sorted_nodes.into_iter().cloned().collect()
}

/// Select UMIs to keep for one (pos, key) group. Returns UMIs whose records to emit.
#[allow(clippy::too_many_lines)]
pub(crate) fn select_umis(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, UmiSlot>,
    edit_threshold: u32,
) -> Vec<Vec<u8>> {
    // Build count and insertion-order maps for sorting (matches Python dict insertion order)
    let counts: HashMap<&[u8], u32> = umi_map
        .iter()
        .map(|(k, v)| (k.as_slice(), v.count))
        .collect();
    let orders: HashMap<&[u8], u32> = umi_map
        .iter()
        .map(|(k, v)| (k.as_slice(), v.insertion_order))
        .collect();
    // Sort key for within-component representative selection: count desc, lex asc.
    // BFS produces lex-sorted components; Python's stable sort preserves that.
    let lex_sort = |a: &[u8], b: &[u8]| -> std::cmp::Ordering {
        counts[b].cmp(&counts[a]).then_with(|| a.cmp(b))
    };

    match method {
        DedupMethod::Unique => {
            // Python returns UMIs in dict insertion order (no count sorting)
            let mut umis: Vec<Vec<u8>> = umi_map.keys().cloned().collect();
            umis.sort_by(|a, b| orders[a.as_slice()].cmp(&orders[b.as_slice()]));
            umis
        }

        DedupMethod::Percentile => {
            if counts.len() <= 1 {
                return umi_map.keys().cloned().collect();
            }
            let all_counts: Vec<u32> = counts.values().copied().collect();
            let threshold = median(&all_counts) / 100.0;
            // Python filters then preserves dict insertion order
            let mut umis: Vec<Vec<u8>> = umi_map
                .iter()
                .filter(|(_, slot)| f64::from(slot.count) > threshold)
                .map(|(umi, _)| umi.clone())
                .collect();
            umis.sort_by(|a, b| orders[a.as_slice()].cmp(&orders[b.as_slice()]));
            umis
        }

        DedupMethod::Cluster => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            // Representative per component: highest count, lex tiebreak
            components
                .into_iter()
                .map(|mut comp| {
                    comp.sort_by(|a, b| lex_sort(a, b));
                    comp.into_iter().next().unwrap()
                })
                .collect()
        }

        DedupMethod::Adjacency => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            let mut result = Vec::new();
            for component in components {
                if component.len() == 1 {
                    result.push(component.into_iter().next().unwrap());
                } else {
                    let lead_umis = min_set_cover(&component, &adj_list, &counts);
                    result.extend(lead_umis);
                }
            }
            result
        }

        DedupMethod::Directional => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_directional_adjacency_list(&umis, &counts, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            let mut observed: HashSet<Vec<u8>> = HashSet::new();
            let mut result = Vec::new();
            for component in components {
                if component.len() == 1 {
                    let umi = component.into_iter().next().unwrap();
                    observed.insert(umi.clone());
                    result.push(umi);
                } else {
                    // Sort by count desc, lex asc (BFS output is lex-sorted,
                    // Python's stable sort preserves that for equal counts)
                    let mut sorted_comp = component;
                    sorted_comp.sort_by(|a, b| lex_sort(a, b));
                    let mut group_lead = None;
                    for node in sorted_comp {
                        if observed.insert(node.clone()) && group_lead.is_none() {
                            group_lead = Some(node);
                        }
                    }
                    if let Some(lead) = group_lead {
                        result.push(lead);
                    }
                }
            }
            result
        }
    }
}

/// Compute the median of a slice of u32 values, returned as f64.
pub(crate) fn median(values: &[u32]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    if n.is_multiple_of(2) {
        f64::midpoint(f64::from(sorted[n / 2 - 1]), f64::from(sorted[n / 2]))
    } else {
        f64::from(sorted[n / 2])
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
    let mut rng = PythonRandom::new(config.random_seed as u32);
    let mut buffer = ReadBuffer::new();
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

        // Subset check consumes one RNG call per mapped read (before buffer.add)
        if config.subset.is_some_and(|s| rng.random() >= f64::from(s)) {
            continue;
        }

        let (start, pos) = get_read_position(&record);

        // Flush buffer when moving far enough or changing chromosome.
        // Matches Python: `if start > last_pos + 1000 or current_chr != last_chr`
        if tid != last_chrom {
            output_records.extend(buffer.drain_all(config.method, config.edit_distance_threshold));
        } else if start > last_start + 1000 {
            let threshold = start - 1000;
            output_records.extend(buffer.drain_up_to(
                threshold,
                config.method,
                config.edit_distance_threshold,
            ));
        }

        last_start = start;
        last_chrom = tid;

        let key: GroupKey = (record.is_reverse(), false, 0, 0);

        let umi = if config.ignore_umi {
            Vec::new()
        } else if config.extract_umi_method == "tag" {
            match extract_umi_from_tag(&record, config.umi_tag.as_deref().unwrap_or("RX")) {
                Some(u) => u,
                None => continue,
            }
        } else {
            extract_umi_from_name(&record, config.umi_separator)
        };

        buffer.add(record, pos, key, umi, &mut rng);
    }

    output_records.extend(buffer.drain_all(config.method, config.edit_distance_threshold));

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

pub(crate) fn extract_umi_from_tag(record: &Record, tag: &str) -> Option<Vec<u8>> {
    match record.aux(tag.as_bytes()) {
        Ok(rust_htslib::bam::record::Aux::String(s)) => Some(s.as_bytes().to_vec()),
        _ => None,
    }
}

pub(crate) fn extract_umi_from_name(record: &Record, separator: u8) -> Vec<u8> {
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
