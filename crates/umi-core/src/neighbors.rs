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
