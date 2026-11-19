//! Exact Hamming neighbors, indexed by disjoint sequence blocks.

use std::collections::{HashMap, HashSet};

pub struct NeighborIndex<'a> {
    sequences: Vec<&'a [u8]>,
    threshold: usize,
    buckets: HashMap<(usize, usize, &'a [u8]), Vec<usize>>,
}

impl<'a> NeighborIndex<'a> {
    pub fn new(sequences: Vec<&'a [u8]>, threshold: usize) -> Self {
        let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
        for (index, sequence) in sequences.iter().enumerate() {
            for (part, block) in blocks(sequence, threshold).enumerate() {
                buckets
                    .entry((sequence.len(), part, block))
                    .or_default()
                    .push(index);
            }
        }
        Self {
            sequences,
            threshold,
            buckets,
        }
    }

    /// Return matching indices, stopping after `limit` matches when only
    /// uniqueness matters. Ordering follows index buckets, not sequence order.
    pub fn matches(&self, query: &[u8], limit: usize) -> Vec<usize> {
        let mut seen = HashSet::new();
        let mut matches = Vec::new();
        if limit == 0 {
            return matches;
        }
        for (part, block) in blocks(query, self.threshold).enumerate() {
            if let Some(candidates) = self.buckets.get(&(query.len(), part, block)) {
                for &index in candidates {
                    if seen.insert(index)
                        && within_distance(query, self.sequences[index], self.threshold)
                    {
                        matches.push(index);
                        if matches.len() == limit {
                            return matches;
                        }
                    }
                }
            }
        }
        matches
    }
}

fn blocks(sequence: &[u8], threshold: usize) -> impl Iterator<Item = &[u8]> {
    // With at most k mismatches, at least one of k+1 disjoint blocks is exact.
    // If k covers the whole length, an empty block includes every equal-length
    // sequence. This also handles empty barcodes and avoids k+1 overflow.
    let parts = if threshold >= sequence.len() {
        1
    } else {
        threshold + 1
    };
    (0..parts).map(move |part| {
        if threshold >= sequence.len() {
            return &sequence[..0];
        }
        let width = sequence.len() / parts;
        let start = part * width;
        let end = if part + 1 == parts {
            sequence.len()
        } else {
            start + width
        };
        &sequence[start..end]
    })
}

fn within_distance(left: &[u8], right: &[u8], threshold: usize) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .filter(|(a, b)| a != b)
            .take(threshold.saturating_add(1))
            .count()
            <= threshold
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sequences() -> Vec<Vec<u8>> {
        let mut sequences: Vec<Vec<u8>> = (0..256_usize)
            .map(|code| (0..4).map(|i| b"ACGT"[(code >> (2 * i)) & 3]).collect())
            .collect();
        sequences.extend([
            b"".to_vec(),
            b"ACGTA".to_vec(),
            b"ACGTT".to_vec(),
            b"TTTTT".to_vec(),
        ]);
        sequences
    }

    #[test]
    fn matches_are_exactly_the_sequences_within_the_threshold() {
        let owned = sequences();
        let all: Vec<&[u8]> = owned.iter().map(Vec::as_slice).collect();
        for threshold in 0..=6 {
            let index = NeighborIndex::new(all.clone(), threshold);
            for query in all.iter().step_by(7).chain(&all[256..]) {
                let expected: Vec<usize> = (0..all.len())
                    .filter(|&i| within_distance(query, all[i], threshold))
                    .collect();
                let mut actual = index.matches(query, usize::MAX);
                actual.sort_unstable();
                assert_eq!(actual, expected, "threshold {threshold}");
            }
        }
    }

    #[test]
    fn limit_stops_the_search_early() {
        let owned = sequences();
        let all: Vec<&[u8]> = owned.iter().map(Vec::as_slice).collect();
        let index = NeighborIndex::new(all.clone(), 1);
        let neighbours = index.matches(b"ACGT", usize::MAX).len();
        assert_eq!(neighbours, 13);
        for limit in 0..=3 {
            let matches = index.matches(b"ACGT", limit);
            assert_eq!(matches.len(), limit);
            assert!(matches.iter().all(|&i| within_distance(b"ACGT", all[i], 1)));
        }
        assert_eq!(index.matches(b"ACGTACGT", usize::MAX), Vec::<usize>::new());
    }

    #[test]
    fn blocks_are_disjoint_and_cover_the_sequence() {
        for (sequence, threshold, expected) in [
            (&b"ACGTAC"[..], 0, vec![&b"ACGTAC"[..]]),
            (b"ACGTAC", 1, vec![b"ACG", b"TAC"]),
            (b"ACGTACG", 2, vec![&b"AC"[..], b"GT", b"ACG"]),
            (b"ACG", 3, vec![b""]),
            (b"", 0, vec![b""]),
        ] {
            assert_eq!(blocks(sequence, threshold).collect::<Vec<_>>(), expected);
        }
    }
}
