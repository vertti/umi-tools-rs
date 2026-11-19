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

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("not all umis are the same length(!):  {min} - {max}")]
pub struct MixedUmiLengths {
    min: usize,
    max: usize,
}

/// # Errors
/// Returns the shortest and longest length when the UMIs differ in length.
pub fn require_equal_lengths<'a>(
    umis: impl IntoIterator<Item = &'a [u8]>,
) -> Result<(), MixedUmiLengths> {
    let mut lengths = umis.into_iter().map(<[u8]>::len);
    let Some(first) = lengths.next() else {
        return Ok(());
    };
    let (min, max) = lengths.fold((first, first), |(min, max), len| {
        (min.min(len), max.max(len))
    });
    if min == max {
        Ok(())
    } else {
        Err(MixedUmiLengths { min, max })
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
    let mut uncovered: HashSet<&[u8]> = cluster.iter().copied().collect();
    for (i, &node) in sorted_nodes[..sorted_nodes.len() - 1].iter().enumerate() {
        uncovered.remove(node);
        for neighbor in adj_list.get(node).into_iter().flatten() {
            uncovered.remove(neighbor);
        }
        if uncovered.is_empty() {
            sorted_nodes.truncate(i + 1);
            break;
        }
    }
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
///
/// # Errors
/// Returns an error when the UMIs differ in length.
pub fn cluster_umis<'a>(
    method: DedupMethod,
    entries: impl IntoIterator<Item = (&'a [u8], u32, u32)>,
    threshold: u32,
) -> Result<Vec<Vec<&'a [u8]>>, MixedUmiLengths> {
    let (counts, orders): (HashMap<_, _>, HashMap<_, _>) = entries
        .into_iter()
        .map(|(umi, count, order)| ((umi, count), (umi, order)))
        .unzip();
    let mut umis: Vec<&[u8]> = counts.keys().copied().collect();
    umis.sort_by_key(|umi| orders[umi]);
    require_equal_lengths(umis.iter().copied())?;
    if matches!(method, DedupMethod::Unique | DedupMethod::Percentile) {
        let cutoff = if method == DedupMethod::Percentile && counts.len() > 1 {
            median(&counts.values().copied().collect::<Vec<_>>()) / 100.0
        } else {
            0.0
        };
        return Ok(umis
            .into_iter()
            .filter(|umi| method == DedupMethod::Unique || f64::from(counts[umi]) > cutoff)
            .map(|umi| vec![umi])
            .collect());
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
    Ok(groups)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rescanning_min_set_cover<'a>(
        cluster: &[&'a [u8]],
        adj_list: &HashMap<&'a [u8], Vec<&'a [u8]>>,
        counts: &HashMap<&[u8], u32>,
    ) -> Vec<&'a [u8]> {
        let mut sorted_nodes = cluster.to_vec();
        sorted_nodes.sort_by(|a, b| counts[*b].cmp(&counts[*a]).then_with(|| a.cmp(b)));
        for i in 0..sorted_nodes.len() - 1 {
            let selected = &sorted_nodes[..=i];
            let covered: HashSet<&[u8]> = selected
                .iter()
                .flat_map(|&node| std::iter::once(node).chain(adj_list[node].iter().copied()))
                .collect();
            if cluster.iter().all(|node| covered.contains(node)) {
                return selected.to_vec();
            }
        }
        sorted_nodes
    }

    #[test]
    fn every_method_rejects_mixed_umi_lengths() {
        for method in [
            DedupMethod::Unique,
            DedupMethod::Percentile,
            DedupMethod::Cluster,
            DedupMethod::Adjacency,
            DedupMethod::Directional,
        ] {
            let entries: [(&[u8], u32, u32); 3] =
                [(b"ACGT", 1, 0), (b"ACGTAA", 1, 1), (b"ACG", 1, 2)];
            let error = cluster_umis(method, entries, 1).unwrap_err();
            assert_eq!(
                error.to_string(),
                "not all umis are the same length(!):  3 - 6"
            );
        }
    }

    /// UMIs with their counts, in first-seen order.
    const BUNDLES: &[&[(&str, u32)]] = &[
        &[
            ("ATAT", 456),
            ("GTAT", 2),
            ("CCAT", 90),
            ("ATAG", 72),
            ("CCAG", 2),
            ("AAAT", 1),
            ("TTTT", 5),
            ("TTTA", 5),
        ],
        &[
            ("AAAA", 10),
            ("AAAT", 9),
            ("AATT", 1),
            ("TTTT", 4),
            ("TTTA", 4),
            ("GGGG", 1),
        ],
        &[("ACGT", 3)],
        &[("ACGT", 1), ("ACGA", 1), ("ACAA", 1), ("AAAA", 1)],
    ];

    /// `UMIClusterer(method)(counts, threshold)` of `umi_tools` on each bundle.
    const EXPECTED_GROUPS: &[(usize, DedupMethod, u32, &[&[&str]])] = &[
        (
            0,
            DedupMethod::Unique,
            1,
            &[
                &["ATAT"],
                &["GTAT"],
                &["CCAT"],
                &["ATAG"],
                &["CCAG"],
                &["AAAT"],
                &["TTTT"],
                &["TTTA"],
            ],
        ),
        (
            0,
            DedupMethod::Percentile,
            1,
            &[
                &["ATAT"],
                &["GTAT"],
                &["CCAT"],
                &["ATAG"],
                &["CCAG"],
                &["AAAT"],
                &["TTTT"],
                &["TTTA"],
            ],
        ),
        (
            0,
            DedupMethod::Cluster,
            1,
            &[
                &["ATAT", "ATAG", "GTAT", "AAAT"],
                &["CCAT", "CCAG"],
                &["TTTA", "TTTT"],
            ],
        ),
        (
            0,
            DedupMethod::Cluster,
            2,
            &[&[
                "ATAT", "CCAT", "ATAG", "TTTA", "TTTT", "CCAG", "GTAT", "AAAT",
            ]],
        ),
        (
            0,
            DedupMethod::Adjacency,
            1,
            &[
                &["ATAT", "ATAG", "GTAT", "AAAT"],
                &["CCAT", "CCAG"],
                &["TTTA", "TTTT"],
            ],
        ),
        (
            0,
            DedupMethod::Adjacency,
            2,
            &[
                &["ATAT", "GTAT", "TTTT", "AAAT"],
                &["CCAT", "CCAG"],
                &["ATAG"],
                &["TTTA"],
            ],
        ),
        (
            0,
            DedupMethod::Directional,
            1,
            &[
                &["ATAT", "ATAG", "GTAT", "AAAT"],
                &["CCAT", "CCAG"],
                &["TTTT"],
                &["TTTA"],
            ],
        ),
        (
            0,
            DedupMethod::Directional,
            2,
            &[
                &["ATAT", "CCAT", "ATAG", "TTTT", "CCAG", "GTAT", "AAAT"],
                &["TTTA"],
            ],
        ),
        (
            1,
            DedupMethod::Unique,
            1,
            &[
                &["AAAA"],
                &["AAAT"],
                &["AATT"],
                &["TTTT"],
                &["TTTA"],
                &["GGGG"],
            ],
        ),
        (
            1,
            DedupMethod::Percentile,
            1,
            &[
                &["AAAA"],
                &["AAAT"],
                &["AATT"],
                &["TTTT"],
                &["TTTA"],
                &["GGGG"],
            ],
        ),
        (
            1,
            DedupMethod::Cluster,
            1,
            &[&["AAAA", "AAAT", "AATT"], &["TTTA", "TTTT"], &["GGGG"]],
        ),
        (
            1,
            DedupMethod::Cluster,
            2,
            &[&["AAAA", "AAAT", "TTTA", "TTTT", "AATT"], &["GGGG"]],
        ),
        (
            1,
            DedupMethod::Adjacency,
            1,
            &[&["AAAA"], &["AAAT", "AATT"], &["TTTA", "TTTT"], &["GGGG"]],
        ),
        (
            1,
            DedupMethod::Adjacency,
            2,
            &[&["AAAA", "AATT"], &["AAAT"], &["TTTA", "TTTT"], &["GGGG"]],
        ),
        (
            1,
            DedupMethod::Directional,
            1,
            &[
                &["AAAA"],
                &["AAAT", "AATT"],
                &["TTTT"],
                &["TTTA"],
                &["GGGG"],
            ],
        ),
        (
            1,
            DedupMethod::Directional,
            2,
            &[
                &["AAAA", "AATT"],
                &["AAAT"],
                &["TTTT"],
                &["TTTA"],
                &["GGGG"],
            ],
        ),
        (2, DedupMethod::Unique, 1, &[&["ACGT"]]),
        (2, DedupMethod::Percentile, 1, &[&["ACGT"]]),
        (2, DedupMethod::Cluster, 1, &[&["ACGT"]]),
        (2, DedupMethod::Cluster, 2, &[&["ACGT"]]),
        (2, DedupMethod::Adjacency, 1, &[&["ACGT"]]),
        (2, DedupMethod::Adjacency, 2, &[&["ACGT"]]),
        (2, DedupMethod::Directional, 1, &[&["ACGT"]]),
        (2, DedupMethod::Directional, 2, &[&["ACGT"]]),
        (
            3,
            DedupMethod::Unique,
            1,
            &[&["ACGT"], &["ACGA"], &["ACAA"], &["AAAA"]],
        ),
        (
            3,
            DedupMethod::Percentile,
            1,
            &[&["ACGT"], &["ACGA"], &["ACAA"], &["AAAA"]],
        ),
        (
            3,
            DedupMethod::Cluster,
            1,
            &[&["AAAA", "ACAA", "ACGA", "ACGT"]],
        ),
        (
            3,
            DedupMethod::Cluster,
            2,
            &[&["AAAA", "ACAA", "ACGA", "ACGT"]],
        ),
        (
            3,
            DedupMethod::Adjacency,
            1,
            &[&["AAAA"], &["ACAA"], &["ACGA", "ACGT"]],
        ),
        (
            3,
            DedupMethod::Adjacency,
            2,
            &[&["AAAA", "ACGA"], &["ACAA", "ACGT"]],
        ),
        (
            3,
            DedupMethod::Directional,
            1,
            &[&["AAAA", "ACAA", "ACGA", "ACGT"]],
        ),
        (
            3,
            DedupMethod::Directional,
            2,
            &[&["AAAA", "ACAA", "ACGA", "ACGT"]],
        ),
    ];

    #[test]
    fn groups_match_umi_tools_for_every_method() {
        for &(bundle, method, threshold, groups) in EXPECTED_GROUPS {
            let entries = BUNDLES[bundle]
                .iter()
                .zip(0..)
                .map(|(&(umi, count), order)| (umi.as_bytes(), count, order));
            let actual: Vec<Vec<&str>> = cluster_umis(method, entries, threshold)
                .unwrap()
                .into_iter()
                .map(|group| {
                    group
                        .into_iter()
                        .map(|umi| std::str::from_utf8(umi).unwrap())
                        .collect()
                })
                .collect();
            // Past the lead UMI, umi_tools lists an adjacency group in set order, which
            // varies with the Python hash seed.
            let normalised = |groups: Vec<Vec<&str>>| -> Vec<Vec<String>> {
                groups
                    .into_iter()
                    .map(|mut group| {
                        if method == DedupMethod::Adjacency {
                            group[1..].sort_unstable();
                        }
                        group.into_iter().map(str::to_string).collect()
                    })
                    .collect()
            };
            assert_eq!(
                normalised(actual),
                normalised(groups.iter().map(|group| group.to_vec()).collect()),
                "bundle {bundle}, {method:?}, threshold {threshold}"
            );
        }
    }

    #[test]
    fn median_of_an_even_number_of_counts_is_the_midpoint() {
        assert!((median(&[1, 2, 3, 10]) - 2.5).abs() < f64::EPSILON);
        assert!((median(&[7, 1, 3]) - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn min_set_cover_matches_rescanning_reference() {
        for size in [1_usize, 2, 3, 50, 500] {
            let names: Vec<Vec<u8>> = (0..size).map(|i| format!("{i:04}").into_bytes()).collect();
            let nodes: Vec<&[u8]> = names.iter().map(Vec::as_slice).collect();
            let neighbours_by_shape: [&dyn Fn(usize) -> Vec<usize>; 4] = [
                &|i| {
                    [i.wrapping_sub(1), i + 1]
                        .into_iter()
                        .filter(|&j| j < size)
                        .collect()
                },
                &|i| if i == 0 { (1..size).collect() } else { vec![0] },
                &|i| (0..size).filter(|&j| j != i).collect(),
                &|_| Vec::new(),
            ];
            for neighbours in neighbours_by_shape {
                let adj_list: HashMap<&[u8], Vec<&[u8]>> = (0..size)
                    .map(|i| {
                        (
                            nodes[i],
                            neighbours(i).into_iter().map(|j| nodes[j]).collect(),
                        )
                    })
                    .collect();
                for count_of in [
                    |_: usize| 1_u32,
                    |i: usize| u32::try_from(i % 7 + 1).unwrap(),
                ] {
                    let counts: HashMap<&[u8], u32> =
                        (0..size).map(|i| (nodes[i], count_of(i))).collect();
                    assert_eq!(
                        min_set_cover(&nodes, &adj_list, &counts),
                        rescanning_min_set_cover(&nodes, &adj_list, &counts)
                    );
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
