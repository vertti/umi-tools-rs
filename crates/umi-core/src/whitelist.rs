use std::collections::HashMap;
use std::io::{BufWriter, Write};

use needletail::parser::{FastqReader, FastxReader};

use crate::error::ExtractError;
use crate::pattern::BarcodePattern;

/// Method for detecting the knee point in the barcode frequency distribution.
#[derive(Debug, Clone, Copy, Default)]
pub enum KneeMethod {
    #[default]
    Distance,
    Density,
}

/// How to handle whitelist barcodes whose edit distance to a higher-count
/// whitelist barcode is within the error-correct threshold.
#[derive(Debug, Clone, Copy)]
pub enum EdAboveThreshold {
    Discard,
    Correct,
}

/// Configuration for the whitelist subcommand.
pub struct WhitelistConfig {
    pub pattern: BarcodePattern,
    pub knee_method: KneeMethod,
    pub cell_number: Option<usize>,
    pub expect_cells: Option<usize>,
    pub error_correct_threshold: usize,
    pub ed_above_threshold: Option<EdAboveThreshold>,
    pub subset_reads: usize,
}

/// A single whitelisted barcode with its count and error-correction mappings.
pub struct WhitelistEntry {
    pub barcode: String,
    pub count: u64,
    pub corrections: Vec<(String, u64)>,
}

/// Statistics from a whitelist run.
pub struct WhitelistStats {
    pub input_reads: u64,
    pub no_match: u64,
}

/// Run the whitelist pipeline: count barcodes, find knee, build error correction, write TSV.
///
/// # Errors
/// Returns error on I/O or pattern-matching failures.
pub fn run_whitelist<R: std::io::Read + Send, W: Write>(
    config: &WhitelistConfig,
    input: R,
    output: W,
) -> Result<WhitelistStats, ExtractError> {
    let (all_counts, stats) = count_barcodes(&config.pattern, input, config.subset_reads)?;

    let whitelist = determine_whitelist(&all_counts, config.knee_method, config.cell_number);

    let corrections =
        build_error_correction_map(&all_counts, &whitelist, config.error_correct_threshold);

    let mut entries: Vec<WhitelistEntry> = whitelist
        .into_iter()
        .map(|bc| {
            let count = all_counts.get(&bc).copied().unwrap_or(0);
            let corr = corrections.get(&bc).cloned().unwrap_or_default();
            WhitelistEntry {
                barcode: bc,
                count,
                corrections: corr,
            }
        })
        .collect();

    entries.sort_by(|a, b| a.barcode.cmp(&b.barcode));

    let mut writer = BufWriter::new(output);
    write_whitelist_tsv(&entries, &mut writer)?;
    writer.flush()?;

    Ok(stats)
}

/// Read FASTQ, extract cell barcodes, count frequencies.
fn count_barcodes<R: std::io::Read + Send>(
    pattern: &BarcodePattern,
    input: R,
    subset_reads: usize,
) -> Result<(HashMap<String, u64>, WhitelistStats), ExtractError> {
    let mut counts: HashMap<String, u64> = HashMap::new();
    let mut stats = WhitelistStats {
        input_reads: 0,
        no_match: 0,
    };

    let mut reader = FastqReader::new(input);

    while let Some(result) = reader.next() {
        let record = result.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
        stats.input_reads += 1;

        if stats.input_reads > subset_reads as u64 {
            break;
        }

        let seq = record.seq();
        let qual = record
            .qual()
            .ok_or_else(|| ExtractError::FastqParse("missing quality scores".into()))?;

        match pattern.extract(&seq, qual) {
            Ok(extraction) => {
                let cell = String::from_utf8_lossy(&extraction.cell_barcode).into_owned();
                if !cell.is_empty() {
                    *counts.entry(cell).or_insert(0) += 1;
                }
            }
            Err(ExtractError::ReadTooShort { .. } | ExtractError::RegexNoMatch) => {
                stats.no_match += 1;
            }
            Err(e) => return Err(e),
        }
    }

    Ok((counts, stats))
}

/// Determine which barcodes to whitelist based on knee detection or explicit cell number.
fn determine_whitelist(
    all_counts: &HashMap<String, u64>,
    knee_method: KneeMethod,
    cell_number: Option<usize>,
) -> Vec<String> {
    let mut sorted_barcodes: Vec<(&String, &u64)> = all_counts.iter().collect();
    sorted_barcodes.sort_by(|a, b| b.1.cmp(a.1));

    if let Some(n) = cell_number {
        if n == 0 || sorted_barcodes.is_empty() {
            return Vec::new();
        }
        let threshold_idx = n.min(sorted_barcodes.len()) - 1;
        let threshold = *sorted_barcodes[threshold_idx].1;
        sorted_barcodes
            .iter()
            .filter(|(_, count)| **count > threshold)
            .map(|(bc, _)| (*bc).clone())
            .collect()
    } else {
        match knee_method {
            KneeMethod::Distance => {
                let counts: Vec<u64> = sorted_barcodes.iter().map(|(_, c)| **c).collect();
                if counts.is_empty() {
                    return Vec::new();
                }
                let knee = knee_distance(&counts);
                sorted_barcodes[..=knee]
                    .iter()
                    .map(|(bc, _)| (*bc).clone())
                    .collect()
            }
            KneeMethod::Density => {
                // Will be implemented in PR 6b
                Vec::new()
            }
        }
    }
}

/// Iterative distance-to-diagonal knee detection on cumulative counts.
fn knee_distance(sorted_desc_counts: &[u64]) -> usize {
    let values = cumulative_sum(sorted_desc_counts);
    let mut prev = 0;
    let mut knee = get_max_distance_index(&values);
    for _ in 0..100 {
        if knee == prev {
            break;
        }
        prev = knee;
        let end = (knee * 3).min(values.len());
        knee = get_max_distance_index(&values[..end]);
    }
    knee
}

/// Find the index with maximum perpendicular distance to the line from first to last point.
#[allow(clippy::cast_precision_loss)]
fn get_max_distance_index(values: &[f64]) -> usize {
    let n = values.len();
    if n <= 1 {
        return 0;
    }

    let first = (0.0_f64, values[0]);
    let last = ((n - 1) as f64, values[n - 1]);
    let line_vec = (last.0 - first.0, last.1 - first.1);
    let line_len = line_vec.0.hypot(line_vec.1);

    if line_len == 0.0 {
        return 0;
    }

    let line_norm = (line_vec.0 / line_len, line_vec.1 / line_len);

    let mut best_dist = 0.0_f64;
    let mut best_idx = 0;
    for (i, &val) in values.iter().enumerate() {
        let v = (i as f64 - first.0, val - first.1);
        let scalar_proj = v.0.mul_add(line_norm.0, v.1 * line_norm.1);
        let parallel = (scalar_proj * line_norm.0, scalar_proj * line_norm.1);
        let perp = (v.0 - parallel.0, v.1 - parallel.1);
        let dist = perp.0.hypot(perp.1);
        if dist > best_dist {
            best_dist = dist;
            best_idx = i;
        }
    }
    best_idx
}

/// Compute cumulative sum as f64 values.
#[allow(clippy::cast_precision_loss)]
fn cumulative_sum(counts: &[u64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(counts.len());
    let mut sum = 0.0_f64;
    for &c in counts {
        sum += c as f64;
        result.push(sum);
    }
    result
}

/// Build error correction map: for each non-whitelist barcode, find if it maps
/// uniquely to exactly one whitelist barcode within the Hamming distance threshold.
fn build_error_correction_map(
    all_counts: &HashMap<String, u64>,
    whitelist: &[String],
    threshold: usize,
) -> HashMap<String, Vec<(String, u64)>> {
    let mut corrections: HashMap<String, Vec<(String, u64)>> = HashMap::new();

    for (barcode, &count) in all_counts {
        if whitelist.contains(barcode) {
            continue;
        }

        let mut matches: Vec<&String> = Vec::new();
        for wl_bc in whitelist {
            if hamming_distance(barcode.as_bytes(), wl_bc.as_bytes()) <= threshold {
                matches.push(wl_bc);
            }
        }

        if matches.len() == 1 {
            corrections
                .entry(matches[0].clone())
                .or_default()
                .push((barcode.clone(), count));
        }
    }

    for corr_list in corrections.values_mut() {
        corr_list.sort_by(|a, b| a.0.cmp(&b.0));
    }

    corrections
}

/// Hamming distance between two byte strings. Returns `usize::MAX` if lengths differ.
fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    if a.len() != b.len() {
        return usize::MAX;
    }
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Write whitelist entries as 4-column TSV.
fn write_whitelist_tsv<W: Write>(
    entries: &[WhitelistEntry],
    writer: &mut W,
) -> Result<(), ExtractError> {
    for entry in entries {
        let error_barcodes: String = entry
            .corrections
            .iter()
            .map(|(bc, _)| bc.as_str())
            .collect::<Vec<_>>()
            .join(",");
        let error_counts: String = entry
            .corrections
            .iter()
            .map(|(_, count)| count.to_string())
            .collect::<Vec<_>>()
            .join(",");

        writer.write_all(entry.barcode.as_bytes())?;
        writer.write_all(b"\t")?;
        writer.write_all(error_barcodes.as_bytes())?;
        writer.write_all(b"\t")?;
        writer.write_all(entry.count.to_string().as_bytes())?;
        writer.write_all(b"\t")?;
        writer.write_all(error_counts.as_bytes())?;
        writer.write_all(b"\n")?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_distance_same() {
        assert_eq!(hamming_distance(b"ACGT", b"ACGT"), 0);
    }

    #[test]
    fn test_hamming_distance_one() {
        assert_eq!(hamming_distance(b"ACGT", b"ACGA"), 1);
    }

    #[test]
    fn test_hamming_distance_different_length() {
        assert_eq!(hamming_distance(b"ACGT", b"ACG"), usize::MAX);
    }

    #[test]
    fn test_cumulative_sum() {
        let counts = vec![10, 5, 3, 1];
        let result = cumulative_sum(&counts);
        assert_eq!(result, vec![10.0, 15.0, 18.0, 19.0]);
    }

    #[test]
    fn test_get_max_distance_index() {
        let values = vec![10.0, 15.0, 18.0, 19.0, 20.0];
        let idx = get_max_distance_index(&values);
        assert!(idx > 0 && idx < values.len() - 1);
    }
}
