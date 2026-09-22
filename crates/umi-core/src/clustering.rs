//! UMI clustering shared by deduplication, grouping, and counting.

use std::collections::{HashMap, HashSet};

use crate::neighbors::NeighborIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DedupMethod {
    Unique,
    Percentile,
    Cluster,
    Adjacency,
    Directional,
}

impl DedupMethod {
    /// Parse a method name shared by deduplication, grouping, and counting.
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "unique" => Some(Self::Unique),
            "percentile" => Some(Self::Percentile),
            "cluster" => Some(Self::Cluster),
            "adjacency" => Some(Self::Adjacency),
            "directional" => Some(Self::Directional),
            _ => None,
        }
    }
}

/// Hamming distance between two byte slices of equal length.
/// Returns `u32::MAX` if lengths differ (matching Python's `np.inf` return).
#[allow(clippy::cast_possible_truncation)]
pub fn hamming_distance(a: &[u8], b: &[u8]) -> u32 {
    if a.len() != b.len() {
        return u32::MAX;
    }
    // UMIs are 5-12bp; count always fits u32
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count() as u32
}

/// Build undirected adjacency list (for cluster + adjacency methods).
/// Edge between A and B iff `hamming_distance(A, B) <= threshold`.
fn build_adjacency_list<'a>(umis: &[&'a [u8]], threshold: u32) -> HashMap<&'a [u8], Vec<&'a [u8]>> {
    indexed_adjacency(umis, threshold, |_, _| true)
}

/// Build directed adjacency list (for directional method).
/// Edge A→B iff `hamming_distance(A,B) <= threshold AND counts[A] >= 2*counts[B] - 1`.
fn build_directional_adjacency_list<'a>(
    umis: &[&'a [u8]],
    counts: &HashMap<&[u8], u32>,
    threshold: u32,
) -> HashMap<&'a [u8], Vec<&'a [u8]>> {
    indexed_adjacency(umis, threshold, |left, right| {
        counts[left] >= (2 * counts[right]).saturating_sub(1)
    })
}

fn indexed_adjacency<'a>(
    umis: &[&'a [u8]],
    threshold: u32,
    accepts: impl Fn(&[u8], &[u8]) -> bool,
) -> HashMap<&'a [u8], Vec<&'a [u8]>> {
    // Tiny bundles need at most 56 comparisons; avoid allocating an index and
    // per-query candidate sets for them. Keep the same input-order neighbors.
    if umis.len() <= 8 {
        return umis
            .iter()
            .enumerate()
            .map(|(i, &umi)| {
                let neighbors = umis
                    .iter()
                    .enumerate()
                    .filter(|&(j, &other)| {
                        i != j
                            && umi.len() == other.len()
                            && hamming_distance(umi, other) <= threshold
                            && accepts(umi, other)
                    })
                    .map(|(_, &other)| other)
                    .collect();
                (umi, neighbors)
            })
            .collect();
    }
    let index = NeighborIndex::new(umis.to_vec(), threshold as usize);
    umis.iter()
        .enumerate()
        .map(|(i, &umi)| {
            let mut neighbors = index.matches(umi, usize::MAX);
            // Preserve input order, including adjacency-method assignment ties.
            neighbors.sort_unstable();
            let neighbors = neighbors
                .into_iter()
                .filter(|&j| j != i && accepts(umi, umis[j]))
                .map(|j| umis[j])
                .collect();
            (umi, neighbors)
        })
        .collect()
}

/// Nodes reachable from `start`, following edges in `adj_list`. Returns the connected component.
fn reachable_nodes<'a>(
    start: &'a [u8],
    adj_list: &HashMap<&'a [u8], Vec<&'a [u8]>>,
) -> Vec<&'a [u8]> {
    let mut searched: HashSet<&'a [u8]> = HashSet::new();
    let mut stack: Vec<&'a [u8]> = Vec::new();
    searched.insert(start);
    stack.push(start);
    while let Some(node) = stack.pop() {
        if let Some(neighbors) = adj_list.get(node) {
            for &next_node in neighbors {
                if searched.insert(next_node) {
                    stack.push(next_node);
                }
            }
        }
    }
    let mut result: Vec<&'a [u8]> = searched.into_iter().collect();
    result.sort();
    result
}

/// Find connected components by iterating UMIs in count-descending order,
/// traversing from each unvisited node. Matches Python `_get_connected_components_adjacency`.
fn connected_components<'a>(
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
            let component = reachable_nodes(umi, adj_list);
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
fn min_set_cover<'a>(
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

pub fn median(values: &[u32]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let n = sorted.len();
    if n.is_multiple_of(2) {
        f64::midpoint(f64::from(sorted[n / 2 - 1]), f64::from(sorted[n / 2]))
    } else {
        f64::from(sorted[n / 2])
    }
}

/// Ordered UMI groups; the first member is the representative.
/// Entries carry counts and first-seen order, independent of alignment storage.
pub fn cluster_umis<'a>(
    method: DedupMethod,
    entries: impl IntoIterator<Item = (&'a [u8], u32, u32)>,
    threshold: u32,
) -> Vec<Vec<&'a [u8]>> {
    let mut entries = entries.into_iter();
    let Some(first) = entries.next() else {
        return Vec::new();
    };
    let Some(second) = entries.next() else {
        return if method == DedupMethod::Percentile && first.1 == 0 {
            Vec::new()
        } else {
            vec![vec![first.0]]
        };
    };
    let entries = [first, second].into_iter().chain(entries);
    let (counts, orders): (HashMap<_, _>, HashMap<_, _>) = entries
        .map(|(umi, count, order)| ((umi, count), (umi, order)))
        .unzip();
    let mut umis: Vec<&[u8]> = counts.keys().copied().collect();
    umis.sort_by_key(|umi| orders[umi]);
    if matches!(method, DedupMethod::Unique | DedupMethod::Percentile) {
        let cutoff = if method == DedupMethod::Percentile && counts.len() > 1 {
            median(&counts.values().copied().collect::<Vec<_>>()) / 100.0
        } else {
            0.0
        };
        return umis
            .into_iter()
            .filter(|umi| method == DedupMethod::Unique || f64::from(counts[umi]) > cutoff)
            .map(|umi| vec![umi])
            .collect();
    }
    let adjacency = if method == DedupMethod::Directional {
        build_directional_adjacency_list(&umis, &counts, threshold)
    } else {
        build_adjacency_list(&umis, threshold)
    };
    let components = connected_components(&umis, &counts, &orders, &adjacency);
    let mut groups = Vec::new();
    let mut observed = HashSet::new();
    for mut component in components {
        component.sort_by(|a, b| counts[b].cmp(&counts[a]).then_with(|| a.cmp(b)));
        match method {
            DedupMethod::Cluster => groups.push(component),
            DedupMethod::Directional => {
                component.retain(|umi| observed.insert(*umi));
                if !component.is_empty() {
                    groups.push(component);
                }
            }
            DedupMethod::Adjacency => {
                let leads = min_set_cover(&component, &adjacency, &counts);
                observed.extend(leads.iter().copied());
                for lead in leads {
                    let mut group = vec![lead];
                    for &neighbor in &adjacency[lead] {
                        if observed.insert(neighbor) {
                            group.push(neighbor);
                        }
                    }
                    groups.push(group);
                }
            }
            DedupMethod::Unique | DedupMethod::Percentile => unreachable!("handled above"),
        }
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_and_singleton_groups_preserve_each_methods_count_rules() {
        let umi = b"ACGT".as_slice();
        for method in [
            DedupMethod::Unique,
            DedupMethod::Percentile,
            DedupMethod::Cluster,
            DedupMethod::Adjacency,
            DedupMethod::Directional,
        ] {
            assert!(cluster_umis(method, [], 1).is_empty());
            for count in [0, 1, 100, u32::MAX] {
                let expected = if method == DedupMethod::Percentile && count == 0 {
                    vec![]
                } else {
                    vec![vec![umi]]
                };
                assert_eq!(cluster_umis(method, [(umi, count, 7)], 1), expected);
            }
            // Repeated entries keep the last count, as the original maps do.
            let expected = if method == DedupMethod::Percentile {
                vec![]
            } else {
                vec![vec![umi]]
            };
            assert_eq!(
                cluster_umis(method, [(umi, 3, 0), (umi, 0, 1)], 1),
                expected
            );
        }
    }

    #[test]
    fn small_adjacency_lists_preserve_input_order_and_directed_edges() {
        let umis: Vec<&[u8]> = vec![
            b"TTT", b"AAC", b"", b"AAA", b"AA", b"ANT", b"NAA", b"TTA", b"ACA",
        ];
        let counts: HashMap<_, _> = umis
            .iter()
            .enumerate()
            .map(|(index, &umi)| (umi, u32::try_from(index % 3 + 1).unwrap()))
            .collect();
        for size in 0..=umis.len() {
            let umis = &umis[..size];
            for threshold in [0, 1, 2, 3, u32::MAX] {
                let undirected = build_adjacency_list(umis, threshold);
                let directed = build_directional_adjacency_list(umis, &counts, threshold);
                for &umi in umis {
                    let expected: Vec<_> = umis
                        .iter()
                        .copied()
                        .filter(|&other| {
                            other != umi
                                && other.len() == umi.len()
                                && hamming_distance(umi, other) <= threshold
                        })
                        .collect();
                    assert_eq!(undirected[umi], expected);
                    let expected: Vec<_> = expected
                        .into_iter()
                        .filter(|other| counts[umi] >= 2 * counts[other] - 1)
                        .collect();
                    assert_eq!(directed[umi], expected);
                }
            }
        }
    }

    #[test]
    fn adjacency_matches_exhaustive_search_for_mixed_barcodes() {
        let mut sequences = std::collections::BTreeSet::from([Vec::new(), b"AAAAA".to_vec()]);
        for mut code in 0..125 {
            let mut sequence = Vec::new();
            for _ in 0..3 {
                sequence.push(b"ACGTN"[code % 5]);
                code /= 5;
            }
            sequences.insert(sequence);
        }
        let umis: Vec<&[u8]> = sequences.iter().map(Vec::as_slice).collect();
        let counts: HashMap<&[u8], u32> = umis
            .iter()
            .enumerate()
            .map(|(index, &umi)| (umi, u32::try_from(index % 7 + 1).unwrap()))
            .collect();
        for threshold in 0..=5 {
            let undirected = build_adjacency_list(&umis, threshold);
            let directed = build_directional_adjacency_list(&umis, &counts, threshold);
            for &umi in &umis {
                let expected: Vec<_> = umis
                    .iter()
                    .copied()
                    .filter(|&other| {
                        other != umi
                            && other.len() == umi.len()
                            && hamming_distance(umi, other) <= threshold
                    })
                    .collect();
                let mut actual = undirected[umi].clone();
                actual.sort_unstable();
                assert_eq!(actual, expected);
                let expected: Vec<_> = expected
                    .into_iter()
                    .filter(|other| counts[umi] >= 2 * counts[other] - 1)
                    .collect();
                let mut actual = directed[umi].clone();
                actual.sort_unstable();
                assert_eq!(actual, expected);
            }
        }
    }
}
