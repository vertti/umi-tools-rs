use std::collections::HashMap;
use std::f64::consts::PI;
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
pub fn run_whitelist<R: std::io::Read + Send, W: Write, FW: Write>(
    config: &WhitelistConfig,
    input: R,
    output: W,
    filtered_out: Option<FW>,
) -> Result<WhitelistStats, ExtractError> {
    let (all_counts, stats) =
        count_barcodes(&config.pattern, input, config.subset_reads, filtered_out)?;

    let whitelist = determine_whitelist(
        &all_counts,
        config.knee_method,
        config.cell_number,
        config.expect_cells,
    );

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
/// Optionally writes non-matching reads to `filtered_out`.
fn count_barcodes<R: std::io::Read + Send, FW: Write>(
    pattern: &BarcodePattern,
    input: R,
    subset_reads: usize,
    filtered_out: Option<FW>,
) -> Result<(HashMap<String, u64>, WhitelistStats), ExtractError> {
    let mut counts: HashMap<String, u64> = HashMap::new();
    let mut stats = WhitelistStats {
        input_reads: 0,
        no_match: 0,
    };
    let mut filt_writer = filtered_out.map(BufWriter::new);

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
                if let Some(fw) = filt_writer.as_mut() {
                    write_fastq_record(fw, record.id(), &seq, qual)?;
                }
            }
            Err(e) => return Err(e),
        }
    }

    if let Some(fw) = filt_writer.as_mut() {
        fw.flush()?;
    }

    Ok((counts, stats))
}

/// Write a FASTQ record (used for filtered-out output).
fn write_fastq_record<W: Write>(
    writer: &mut W,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> Result<(), ExtractError> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// Determine which barcodes to whitelist based on knee detection or explicit cell number.
fn determine_whitelist(
    all_counts: &HashMap<String, u64>,
    knee_method: KneeMethod,
    cell_number: Option<usize>,
    expect_cells: Option<usize>,
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
            KneeMethod::Density => knee_density(&sorted_barcodes, expect_cells),
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

/// Density-based knee detection using Gaussian KDE on log10-transformed counts.
/// Matches scipy's `gaussian_kde(data, bw_method=0.1)` behavior.
#[allow(clippy::cast_precision_loss)]
fn knee_density(sorted_barcodes: &[(&String, &u64)], expect_cells: Option<usize>) -> Vec<String> {
    if sorted_barcodes.is_empty() {
        return Vec::new();
    }

    let max_count = *sorted_barcodes[0].1 as f64;
    let abundance_threshold = max_count * 0.001;

    // Log10-transform counts above abundance threshold
    let log_counts: Vec<f64> = sorted_barcodes
        .iter()
        .map(|(_, c)| **c as f64)
        .filter(|&c| c > abundance_threshold)
        .map(f64::log10)
        .collect();

    if log_counts.is_empty() {
        return Vec::new();
    }

    let bw = sample_std(&log_counts) * 0.1;
    if bw <= 0.0 {
        return Vec::new();
    }

    let log_min = log_counts.iter().copied().fold(f64::INFINITY, f64::min);
    let log_max = log_counts.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let num_points: usize = 10_000;
    let xx: Vec<f64> = (0..num_points)
        .map(|i| (log_max - log_min).mul_add(i as f64 / (num_points - 1) as f64, log_min))
        .collect();

    let density = gaussian_kde(&log_counts, bw, &xx);

    // Find local minima: density[i] < density[i-1] && density[i] < density[i+1]
    let local_mins: Vec<usize> = (1..density.len() - 1)
        .filter(|&i| density[i] < density[i - 1] && density[i] < density[i + 1])
        .collect();

    if local_mins.is_empty() {
        return Vec::new();
    }

    // Select the appropriate local minimum by iterating in reverse
    let mut selected_min: Option<usize> = None;
    for &min_idx in local_mins.iter().rev() {
        let threshold = 10.0_f64.powf(xx[min_idx]);
        let passing_count = sorted_barcodes
            .iter()
            .filter(|(_, c)| **c as f64 > threshold)
            .count();

        if let Some(expected) = expect_cells {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let lo = (expected as f64 * 0.1) as usize;
            if passing_count > lo && passing_count <= expected {
                selected_min = Some(min_idx);
                break;
            }
        } else {
            let xx_values = xx.len();
            let at_least_20pct = min_idx as f64 >= 0.2 * xx_values as f64;
            let far_from_max = log_max - xx[min_idx] > 0.5;
            let below_half_max = xx[min_idx] < log_max / 2.0;

            if at_least_20pct && (far_from_max || below_half_max) {
                selected_min = Some(min_idx);
                break;
            }
        }
    }

    let Some(min_idx) = selected_min else {
        return Vec::new();
    };

    let threshold = 10.0_f64.powf(xx[min_idx]);
    sorted_barcodes
        .iter()
        .filter(|(_, c)| **c as f64 > threshold)
        .map(|(bc, _)| (*bc).clone())
        .collect()
}

/// Gaussian KDE evaluation matching scipy's `gaussian_kde` behavior.
#[allow(clippy::cast_precision_loss)]
fn gaussian_kde(data: &[f64], bw: f64, points: &[f64]) -> Vec<f64> {
    let n = data.len() as f64;
    let coeff = 1.0 / (n * bw * (2.0 * PI).sqrt());
    points
        .iter()
        .map(|&x| {
            coeff
                * data
                    .iter()
                    .map(|&d| {
                        let z = (x - d) / bw;
                        (-0.5 * z * z).exp()
                    })
                    .sum::<f64>()
        })
        .collect()
}

/// Sample standard deviation (ddof=1, Bessel's correction).
#[allow(clippy::cast_precision_loss)]
fn sample_std(data: &[f64]) -> f64 {
    let n = data.len() as f64;
    if n <= 1.0 {
        return 0.0;
    }
    let mean = data.iter().sum::<f64>() / n;
    let var = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
    var.sqrt()
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

    #[test]
    fn test_sample_std() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let s = sample_std(&data);
        // Expected: sqrt(32/7) ≈ 2.138
        assert!((s - 2.138).abs() < 0.01);
    }

    #[test]
    fn test_gaussian_kde_single_point() {
        let data = vec![0.0];
        let bw = 1.0;
        let points = vec![0.0];
        let result = gaussian_kde(&data, bw, &points);
        // At data point with bw=1: 1/(1*1*sqrt(2*pi)) * exp(0) ≈ 0.3989
        assert!((result[0] - 0.3989).abs() < 0.001);
    }
}
