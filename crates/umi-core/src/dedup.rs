use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, Write as IoWrite};

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DedupMethod {
    Unique,
    Percentile,
    Cluster,
    Adjacency,
    Directional,
}

#[allow(clippy::struct_excessive_bools)]
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
    pub per_gene: bool,
    pub gene_tag: Option<String>,
    pub skip_tags_regex: Option<String>,
    pub output_stats: Option<String>,
    pub paired: bool,
    pub ignore_tlen: bool,
    pub umi_whitelist: Option<HashSet<Vec<u8>>>,
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
        groups: BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, UmiSlot>>>,
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
                            ctx.umi_separator,
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
    umi_separator: u8,
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
pub(crate) fn build_adjacency_list<'a>(
    umis: &[&'a [u8]],
    threshold: u32,
) -> HashMap<&'a [u8], Vec<&'a [u8]>> {
    let mut adj: HashMap<&'a [u8], Vec<&'a [u8]>> = HashMap::new();
    for umi in umis {
        adj.entry(umi).or_default();
    }
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            if edit_distance(umis[i], umis[j]) <= threshold {
                adj.get_mut(umis[i])
                    .expect("UMI pre-inserted")
                    .push(umis[j]);
                adj.get_mut(umis[j])
                    .expect("UMI pre-inserted")
                    .push(umis[i]);
            }
        }
    }
    adj
}

/// Build directed adjacency list (for directional method).
/// Edge A→B iff `edit_distance(A,B) <= threshold AND counts[A] >= 2*counts[B] - 1`.
pub(crate) fn build_directional_adjacency_list<'a>(
    umis: &[&'a [u8]],
    counts: &HashMap<&[u8], u32>,
    threshold: u32,
) -> HashMap<&'a [u8], Vec<&'a [u8]>> {
    let mut adj: HashMap<&'a [u8], Vec<&'a [u8]>> = HashMap::new();
    for umi in umis {
        adj.entry(umi).or_default();
    }
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            if edit_distance(umis[i], umis[j]) <= threshold {
                let ca = counts[umis[i]];
                let cb = counts[umis[j]];
                if ca >= (2 * cb).saturating_sub(1) {
                    adj.get_mut(umis[i])
                        .expect("UMI pre-inserted")
                        .push(umis[j]);
                }
                if cb >= (2 * ca).saturating_sub(1) {
                    adj.get_mut(umis[j])
                        .expect("UMI pre-inserted")
                        .push(umis[i]);
                }
            }
        }
    }
    adj
}

/// BFS from `start`, following edges in `adj_list`. Returns the connected component.
pub(crate) fn bfs<'a>(
    start: &'a [u8],
    adj_list: &HashMap<&'a [u8], Vec<&'a [u8]>>,
) -> Vec<&'a [u8]> {
    let mut searched: HashSet<&'a [u8]> = HashSet::new();
    let mut queue: Vec<&'a [u8]> = Vec::new();
    searched.insert(start);
    queue.push(start);
    while let Some(node) = queue.pop() {
        if let Some(neighbors) = adj_list.get(node) {
            for &next_node in neighbors {
                if searched.insert(next_node) {
                    queue.push(next_node);
                }
            }
        }
    }
    let mut result: Vec<&'a [u8]> = searched.into_iter().collect();
    result.sort();
    result
}

/// Find connected components by iterating UMIs in count-descending order,
/// running BFS from each unvisited node. Matches Python `_get_connected_components_adjacency`.
pub(crate) fn connected_components<'a>(
    umis: &[&'a [u8]],
    counts: &HashMap<&[u8], u32>,
    orders: &HashMap<&[u8], u32>,
    adj_list: &HashMap<&'a [u8], Vec<&'a [u8]>>,
) -> Vec<Vec<&'a [u8]>> {
    // Sort UMIs by count descending, then insertion order ascending for ties
    let mut sorted_umis: Vec<&[u8]> = umis.to_vec();
    sorted_umis.sort_by(|a, b| {
        counts[b]
            .cmp(&counts[a])
            .then_with(|| orders[a].cmp(&orders[b]))
    });

    let mut found: HashSet<&[u8]> = HashSet::new();
    let mut components: Vec<Vec<&'a [u8]>> = Vec::new();
    for umi in &sorted_umis {
        if !found.contains(*umi) {
            let component = bfs(umi, adj_list);
            for &node in &component {
                found.insert(node);
            }
            components.push(component);
        }
    }
    components
}

/// Greedy min-set-cover: select fewest UMIs (by descending count) to "cover"
/// all UMIs in the cluster via adjacency. Matches Python `_get_best_min_account`.
pub(crate) fn min_set_cover<'a>(
    cluster: &[&'a [u8]],
    adj_list: &HashMap<&'a [u8], Vec<&'a [u8]>>,
    counts: &HashMap<&[u8], u32>,
) -> Vec<&'a [u8]> {
    if cluster.len() == 1 {
        return cluster.to_vec();
    }
    let mut sorted_nodes: Vec<&'a [u8]> = cluster.to_vec();
    // Sort by count desc, lex asc (BFS output is lex-sorted; Python's stable sort preserves that)
    sorted_nodes.sort_by(|a, b| counts[*b].cmp(&counts[*a]).then_with(|| a.cmp(b)));
    for i in 0..sorted_nodes.len() - 1 {
        let selected = &sorted_nodes[..=i];
        // Compute covered nodes: selected nodes + their neighbors
        let mut covered: HashSet<&[u8]> = HashSet::new();
        for &s in selected {
            covered.insert(s);
            if let Some(neighbors) = adj_list.get(s) {
                for &n in neighbors {
                    covered.insert(n);
                }
            }
        }
        // Check if all cluster nodes are covered
        let remaining: usize = cluster.iter().filter(|n| !covered.contains(*n)).count();
        if remaining == 0 {
            return selected.to_vec();
        }
    }
    // Fallback: all nodes (shouldn't reach here for valid inputs)
    sorted_nodes
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
                    comp.into_iter()
                        .next()
                        .expect("component is non-empty")
                        .to_vec()
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
                    result.push(component[0].to_vec());
                } else {
                    let lead_umis = min_set_cover(&component, &adj_list, &counts);
                    result.extend(lead_umis.into_iter().map(<[u8]>::to_vec));
                }
            }
            result
        }

        DedupMethod::Directional => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_directional_adjacency_list(&umis, &counts, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            let mut observed: HashSet<&[u8]> = HashSet::new();
            let mut result = Vec::new();
            for component in components {
                if component.len() == 1 {
                    let umi = component[0];
                    observed.insert(umi);
                    result.push(umi.to_vec());
                } else {
                    // Sort by count desc, lex asc (BFS output is lex-sorted,
                    // Python's stable sort preserves that for equal counts)
                    let mut sorted_comp = component;
                    sorted_comp.sort_by(|a, b| lex_sort(a, b));
                    let mut group_lead = None;
                    for node in sorted_comp {
                        if observed.insert(node) && group_lead.is_none() {
                            group_lead = Some(node);
                        }
                    }
                    if let Some(lead) = group_lead {
                        result.push(lead.to_vec());
                    }
                }
            }
            result
        }
    }
}

/// Count deduplicated UMI groups from raw count/order maps.
///
/// Same logic as `select_umis` but takes `HashMap<Vec<u8>, u32>` instead of
/// `UmiSlot`, and returns only the count of surviving UMI groups.
#[allow(clippy::implicit_hasher)]
#[must_use]
pub fn count_umis(
    method: DedupMethod,
    counts: &HashMap<Vec<u8>, u32>,
    orders: &HashMap<Vec<u8>, u32>,
    edit_threshold: u32,
) -> usize {
    let count_refs: HashMap<&[u8], u32> = counts.iter().map(|(k, v)| (k.as_slice(), *v)).collect();
    let order_refs: HashMap<&[u8], u32> = orders.iter().map(|(k, v)| (k.as_slice(), *v)).collect();
    let lex_sort = |a: &[u8], b: &[u8]| -> std::cmp::Ordering {
        count_refs[b].cmp(&count_refs[a]).then_with(|| a.cmp(b))
    };

    match method {
        DedupMethod::Unique => counts.len(),

        DedupMethod::Percentile => {
            if counts.len() <= 1 {
                return counts.len();
            }
            let all_counts: Vec<u32> = counts.values().copied().collect();
            let threshold = median(&all_counts) / 100.0;
            counts
                .values()
                .filter(|&&c| f64::from(c) > threshold)
                .count()
        }

        DedupMethod::Cluster => {
            let umis: Vec<&[u8]> = counts.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &count_refs, &order_refs, &adj_list);
            components.len()
        }

        DedupMethod::Adjacency => {
            let umis: Vec<&[u8]> = counts.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &count_refs, &order_refs, &adj_list);
            let mut total = 0;
            for component in components {
                if component.len() == 1 {
                    total += 1;
                } else {
                    total += min_set_cover(&component, &adj_list, &count_refs).len();
                }
            }
            total
        }

        DedupMethod::Directional => {
            let umis: Vec<&[u8]> = counts.keys().map(Vec::as_slice).collect();
            let adj_list = build_directional_adjacency_list(&umis, &count_refs, edit_threshold);
            let components = connected_components(&umis, &count_refs, &order_refs, &adj_list);
            let mut observed: HashSet<&[u8]> = HashSet::new();
            let mut total = 0;
            for component in components {
                if component.len() == 1 {
                    let umi = component[0];
                    observed.insert(umi);
                    total += 1;
                } else {
                    let mut sorted_comp = component;
                    sorted_comp.sort_by(|a, b| lex_sort(a, b));
                    let mut found_lead = false;
                    for node in sorted_comp {
                        if observed.insert(node) && !found_lead {
                            found_lead = true;
                            total += 1;
                        }
                    }
                }
            }
            total
        }
    }
}

/// Extract UMI and optional cell barcode from a read name using the `umis` method.
///
/// Splits `qname` by `:` and looks for `UMI_<seq>` and `CELL_<barcode>` prefixed
/// fields. Returns `(umi, Option<cell>)`.
#[must_use]
pub fn extract_umi_umis(qname: &[u8]) -> (Vec<u8>, Option<Vec<u8>>) {
    let mut umi = None;
    let mut cell = None;
    for part in qname.split(|&b| b == b':') {
        if part.starts_with(b"UMI_") {
            umi = Some(part[4..].to_vec());
        } else if part.starts_with(b"CELL_") {
            cell = Some(part[5..].to_vec());
        }
    }
    (umi.unwrap_or_default(), cell)
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

/// Like `select_umis`, but also returns the total count for each cluster
/// (sum of all UMI counts in the cluster, not just the representative).
/// Returns `(selected_umi, cluster_total_count)` pairs.
#[allow(clippy::too_many_lines)]
fn select_umis_with_cluster_counts(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, UmiSlot>,
    edit_threshold: u32,
) -> Vec<(Vec<u8>, u32)> {
    let counts: HashMap<&[u8], u32> = umi_map
        .iter()
        .map(|(k, v)| (k.as_slice(), v.count))
        .collect();
    let orders: HashMap<&[u8], u32> = umi_map
        .iter()
        .map(|(k, v)| (k.as_slice(), v.insertion_order))
        .collect();
    let lex_sort = |a: &[u8], b: &[u8]| -> std::cmp::Ordering {
        counts[b].cmp(&counts[a]).then_with(|| a.cmp(b))
    };

    match method {
        DedupMethod::Unique => {
            let mut umis: Vec<Vec<u8>> = umi_map.keys().cloned().collect();
            umis.sort_by(|a, b| orders[a.as_slice()].cmp(&orders[b.as_slice()]));
            umis.into_iter()
                .map(|u| {
                    let c = counts[u.as_slice()];
                    (u, c)
                })
                .collect()
        }

        DedupMethod::Percentile => {
            if counts.len() <= 1 {
                return umi_map.iter().map(|(u, s)| (u.clone(), s.count)).collect();
            }
            let all_counts: Vec<u32> = counts.values().copied().collect();
            let threshold = median(&all_counts) / 100.0;
            let mut umis: Vec<Vec<u8>> = umi_map
                .iter()
                .filter(|(_, slot)| f64::from(slot.count) > threshold)
                .map(|(umi, _)| umi.clone())
                .collect();
            umis.sort_by(|a, b| orders[a.as_slice()].cmp(&orders[b.as_slice()]));
            umis.into_iter()
                .map(|u| {
                    let c = counts[u.as_slice()];
                    (u, c)
                })
                .collect()
        }

        DedupMethod::Cluster => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            components
                .into_iter()
                .map(|mut comp| {
                    let cluster_count: u32 = comp.iter().map(|u| counts[*u]).sum();
                    comp.sort_by(|a, b| lex_sort(a, b));
                    (
                        comp.into_iter()
                            .next()
                            .expect("component is non-empty")
                            .to_vec(),
                        cluster_count,
                    )
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
                    let c = counts[component[0]];
                    result.push((component[0].to_vec(), c));
                } else {
                    let lead_umis = min_set_cover(&component, &adj_list, &counts);
                    // Each lead UMI's cluster: itself + its unobserved neighbors
                    let mut observed: HashSet<&[u8]> = HashSet::new();
                    for &lead in &lead_umis {
                        let mut cluster_count = counts[lead];
                        observed.insert(lead);
                        if let Some(neighbors) = adj_list.get(lead) {
                            for &n in neighbors {
                                if observed.insert(n) {
                                    cluster_count += counts[n];
                                }
                            }
                        }
                        result.push((lead.to_vec(), cluster_count));
                    }
                }
            }
            result
        }

        DedupMethod::Directional => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_directional_adjacency_list(&umis, &counts, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            let mut observed: HashSet<&[u8]> = HashSet::new();
            let mut result = Vec::new();
            for component in components {
                if component.len() == 1 {
                    let umi = component[0];
                    let c = counts[umi];
                    observed.insert(umi);
                    result.push((umi.to_vec(), c));
                } else {
                    let mut sorted_comp = component;
                    sorted_comp.sort_by(|a, b| lex_sort(a, b));
                    let mut group_lead = None;
                    let mut cluster_count: u32 = 0;
                    for node in sorted_comp {
                        if observed.insert(node) {
                            cluster_count += counts[node];
                            if group_lead.is_none() {
                                group_lead = Some(node);
                            }
                        }
                    }
                    if let Some(lead) = group_lead {
                        result.push((lead.to_vec(), cluster_count));
                    }
                }
            }
            result
        }
    }
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
            total += u64::from(edit_distance(umis[i], umis[j]));
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
        umi_separator: u8,
        extract_method: &str,
        umi_tag: Option<&str>,
        chrom: Option<&str>,
        seed: u32,
    ) -> Result<Self, DedupError> {
        let mut reader =
            bam::Reader::from_path(bam_path).map_err(|e| DedupError::BamOpen(e.to_string()))?;

        let chrom_tid: Option<i32> = chrom
            .map(|c| {
                let tid = reader
                    .header()
                    .tid(c.as_bytes())
                    .ok_or_else(|| DedupError::UnknownChrom(c.to_string()))?;
                #[allow(clippy::cast_possible_wrap)]
                Ok(tid as i32)
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
            let umi = if extract_method == "tag" {
                match extract_umi_from_tag(&record, umi_tag.unwrap_or("RX")) {
                    Some(u) => u,
                    None => continue,
                }
            } else {
                extract_umi_from_name(&record, umi_separator)
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
        umi_separator: u8,
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
            .map(|r| extract_umi_from_name(r, umi_separator))
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

    // Per-gene mode: gene tag value → sequential i64 ID used as "position".
    let skip_regex = config
        .skip_tags_regex
        .as_ref()
        .map(|s| regex::Regex::new(s).map_err(|e| DedupError::InvalidRegex(e.to_string())))
        .transpose()?;
    let mut gene_ids: HashMap<Vec<u8>, i64> = HashMap::new();
    let mut next_gene_id: i64 = 0;

    // Stats collection (optional, only when --output-stats is set)
    #[allow(clippy::cast_possible_truncation)]
    let mut stats_ctx: Option<StatsContext> = config
        .output_stats
        .as_ref()
        .map(|_| {
            let read_gen = RandomReadGenerator::new(
                input_path,
                config.umi_separator,
                &config.extract_umi_method,
                config.umi_tag.as_deref(),
                config.chrom.as_deref(),
                config.random_seed as u32,
            )?;
            Ok(StatsContext {
                collector: StatsCollector::new(),
                read_gen,
                umi_separator: config.umi_separator,
            })
        })
        .transpose()?;

    let wl_ref = config.umi_whitelist.as_ref();

    for result in reader.records() {
        let record = result.map_err(|e| DedupError::BamRead(e.to_string()))?;

        if record.is_unmapped() {
            continue;
        }

        // Paired mode: skip R2 reads and R1s with unmapped mates.
        if config.paired {
            if record.is_last_in_template() {
                continue;
            }
            if record.is_mate_unmapped() {
                continue;
            }
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

        if config.per_gene {
            // Per-gene mode: group by gene tag value instead of position.
            let gene_tag_name = config.gene_tag.as_deref().unwrap_or("XF");
            let Some(gene) = extract_umi_from_tag(&record, gene_tag_name) else {
                continue;
            };
            if skip_regex
                .as_ref()
                .is_some_and(|re| re.is_match(std::str::from_utf8(&gene).unwrap_or("")))
            {
                continue;
            }
            let gene_id = *gene_ids.entry(gene).or_insert_with(|| {
                let id = next_gene_id;
                next_gene_id += 1;
                id
            });
            buffer.add(record, gene_id, (false, false, 0, 0), umi, &mut rng);
        } else {
            let (start, pos) = get_read_position(&record);

            // Flush buffer when moving far enough or changing chromosome.
            if tid != last_chrom {
                output_records.extend(buffer.drain_all(
                    config.method,
                    config.edit_distance_threshold,
                    &mut stats_ctx,
                    wl_ref,
                ));
            } else if start > last_start + 1000 {
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

            let tlen = if config.paired && !config.ignore_tlen {
                record.insert_size()
            } else {
                0
            };
            let key: GroupKey = (record.is_reverse(), false, tlen, 0);
            buffer.add(record, pos, key, umi, &mut rng);
        }
    }

    output_records.extend(buffer.drain_all(
        config.method,
        config.edit_distance_threshold,
        &mut stats_ctx,
        wl_ref,
    ));

    // Paired mode: second pass to find R2 mates of surviving R1 reads.
    if config.paired {
        let mut mate_set: HashSet<(Vec<u8>, i32, i64)> = HashSet::new();
        for r1 in &output_records {
            mate_set.insert((r1.qname().to_vec(), r1.mtid(), r1.mpos()));
        }
        let mut reader2 =
            bam::Reader::from_path(input_path).map_err(|e| DedupError::BamOpen(e.to_string()))?;
        for result in reader2.records() {
            let record = result.map_err(|e| DedupError::BamRead(e.to_string()))?;
            if record.is_unmapped() || record.is_mate_unmapped() {
                continue;
            }
            if !record.is_last_in_template() {
                continue;
            }
            let key = (record.qname().to_vec(), record.tid(), record.pos());
            if mate_set.remove(&key) {
                output_records.push(record);
            }
        }
    }

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
    #[error("invalid regex: {0}")]
    InvalidRegex(String),
    #[error("failed to write stats file {0}: {1}")]
    StatsWrite(String, String),
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
