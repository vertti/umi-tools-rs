use std::ops::Range;

use crate::error::ExtractError;

/// Result of extracting barcodes from a single read's sequence and quality.
#[derive(Debug, Clone)]
pub struct ExtractionResult {
    pub umi: Vec<u8>,
    pub umi_quality: Vec<u8>,
    pub cell_barcode: Vec<u8>,
    pub trimmed_sequence: Vec<u8>,
    pub trimmed_quality: Vec<u8>,
}

/// Which end of the read to extract the barcode from (string method only).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimeEnd {
    Five,
    Three,
}

/// A parsed barcode pattern that knows how to extract UMI/cell/sample bases from a read.
#[derive(Debug, Clone)]
pub enum BarcodePattern {
    String(StringPattern),
    Regex(RegexPattern),
}

impl BarcodePattern {
    /// # Errors
    /// Returns error if the read is too short (string method) or doesn't match (regex method).
    pub fn extract(
        &self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<ExtractionResult, ExtractError> {
        match self {
            Self::String(p) => p.extract(sequence, quality),
            Self::Regex(p) => p.extract(sequence, quality),
        }
    }
}

/// String-method pattern using fixed-position characters.
///
/// Pattern characters:
/// - `N` — UMI base (extracted to read name)
/// - `C` — Cell barcode base (extracted to read name)
/// - `X` — Sample/discard base (stays in output sequence, removed from barcode region)
#[derive(Debug, Clone)]
pub struct StringPattern {
    umi_positions: Vec<usize>,
    cell_positions: Vec<usize>,
    sample_positions: Vec<usize>,
    umi_range: Option<Range<usize>>,
    cell_range: Option<Range<usize>>,
    sample_range: Option<Range<usize>>,
    pattern_length: usize,
    prime_end: PrimeEnd,
}

impl StringPattern {
    /// Parse a string-method pattern like `NNNXXXXNN`.
    ///
    /// # Errors
    /// Returns error if pattern is empty or contains characters other than N, X, C.
    pub fn parse(pattern_str: &str, prime_end: PrimeEnd) -> Result<Self, ExtractError> {
        if pattern_str.is_empty() {
            return Err(ExtractError::InvalidPattern(
                "pattern must not be empty".into(),
            ));
        }

        let mut umi_positions = Vec::new();
        let mut cell_positions = Vec::new();
        let mut sample_positions = Vec::new();

        for (i, ch) in pattern_str.chars().enumerate() {
            match ch {
                'N' => umi_positions.push(i),
                'C' => cell_positions.push(i),
                'X' => sample_positions.push(i),
                other => {
                    return Err(ExtractError::InvalidPattern(format!(
                        "pattern contains invalid character '{other}' at position {i}; \
                         only N, X, C are allowed"
                    )));
                }
            }
        }

        let umi_range = as_contiguous_range(&umi_positions);
        let cell_range = as_contiguous_range(&cell_positions);
        let sample_range = as_contiguous_range(&sample_positions);

        Ok(Self {
            umi_positions,
            cell_positions,
            sample_positions,
            umi_range,
            cell_range,
            sample_range,
            pattern_length: pattern_str.len(),
            prime_end,
        })
    }

    /// Extract barcodes from a read's sequence and quality strings.
    ///
    /// # Errors
    /// Returns error if the read is shorter than the pattern.
    pub fn extract(
        &self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<ExtractionResult, ExtractError> {
        if sequence.len() < self.pattern_length {
            return Err(ExtractError::ReadTooShort {
                read_len: sequence.len(),
                pattern_len: self.pattern_length,
            });
        }

        let (barcode_region, remaining_seq, barcode_qual, remaining_qual) = match self.prime_end {
            PrimeEnd::Five => (
                &sequence[..self.pattern_length],
                &sequence[self.pattern_length..],
                &quality[..self.pattern_length],
                &quality[self.pattern_length..],
            ),
            PrimeEnd::Three => (
                &sequence[sequence.len() - self.pattern_length..],
                &sequence[..sequence.len() - self.pattern_length],
                &quality[quality.len() - self.pattern_length..],
                &quality[..quality.len() - self.pattern_length],
            ),
        };

        let umi = extract_slice(barcode_region, self.umi_range.as_ref(), &self.umi_positions);
        let umi_quality = extract_slice(barcode_qual, self.umi_range.as_ref(), &self.umi_positions);
        let cell_barcode = extract_slice(
            barcode_region,
            self.cell_range.as_ref(),
            &self.cell_positions,
        );

        let (trimmed_sequence, trimmed_quality) = if self.sample_positions.is_empty() {
            (remaining_seq.to_vec(), remaining_qual.to_vec())
        } else {
            let sample_seq = extract_slice(
                barcode_region,
                self.sample_range.as_ref(),
                &self.sample_positions,
            );
            let sample_qual = extract_slice(
                barcode_qual,
                self.sample_range.as_ref(),
                &self.sample_positions,
            );
            match self.prime_end {
                PrimeEnd::Five => (
                    join_slices(&sample_seq, remaining_seq),
                    join_slices(&sample_qual, remaining_qual),
                ),
                PrimeEnd::Three => (
                    join_slices(remaining_seq, &sample_seq),
                    join_slices(remaining_qual, &sample_qual),
                ),
            }
        };

        Ok(ExtractionResult {
            umi,
            umi_quality,
            cell_barcode,
            trimmed_sequence,
            trimmed_quality,
        })
    }
}

/// Regex-method pattern using named capture groups.
///
/// Groups starting with `umi_` are extracted as UMI, `cell_` as cell barcode,
/// `discard_` as bases to remove. Everything else is kept in the output sequence.
#[derive(Debug, Clone)]
pub struct RegexPattern {
    pattern: regex::Regex,
}

impl RegexPattern {
    /// Parse a regex pattern string.
    ///
    /// # Errors
    /// Returns error if the regex is invalid or has no `umi_` or `cell_` groups.
    pub fn parse(pattern_str: &str) -> Result<Self, ExtractError> {
        let processed = preprocess_fuzzy(pattern_str)?;

        let pattern = regex::Regex::new(&processed)
            .map_err(|e| ExtractError::InvalidPattern(format!("invalid regex: {e}")))?;

        let has_barcode_group = pattern
            .capture_names()
            .flatten()
            .any(|name| name.starts_with("umi_") || name.starts_with("cell_"));

        if !has_barcode_group {
            return Err(ExtractError::InvalidPattern(
                "regex must contain at least one named group starting with 'umi_' or 'cell_'"
                    .into(),
            ));
        }

        Ok(Self { pattern })
    }

    /// Extract barcodes from a read's sequence and quality strings.
    ///
    /// # Errors
    /// Returns `RegexNoMatch` if the regex doesn't match the sequence.
    pub fn extract(
        &self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<ExtractionResult, ExtractError> {
        let seq_str = std::str::from_utf8(sequence)
            .map_err(|e| ExtractError::FastqParse(format!("non-UTF8 sequence: {e}")))?;

        let caps = self
            .pattern
            .captures(seq_str)
            .ok_or(ExtractError::RegexNoMatch)?;

        // Collect named group spans into (name, start, end) sorted by name
        let mut umi_spans: Vec<(&str, usize, usize)> = Vec::new();
        let mut cell_spans: Vec<(&str, usize, usize)> = Vec::new();
        let mut discard_spans: Vec<(usize, usize)> = Vec::new();

        for name in self.pattern.capture_names().flatten() {
            if let Some(m) = caps.name(name) {
                let span = (m.start(), m.end());
                if name.starts_with("umi_") {
                    umi_spans.push((name, span.0, span.1));
                } else if name.starts_with("cell_") {
                    cell_spans.push((name, span.0, span.1));
                } else if name.starts_with("discard_") {
                    discard_spans.push(span);
                }
            }
        }

        // Sort by group name for deterministic concatenation
        umi_spans.sort_by_key(|&(name, _, _)| name);
        cell_spans.sort_by_key(|&(name, _, _)| name);

        // Build extracted-position bitmask (O(n) lookup instead of O(n*m) Vec::contains)
        let mut extracted = vec![false; sequence.len()];
        for &(_, start, end) in &umi_spans {
            extracted[start..end].fill(true);
        }
        for &(_, start, end) in &cell_spans {
            extracted[start..end].fill(true);
        }
        for &(start, end) in &discard_spans {
            extracted[start..end].fill(true);
        }

        // Build UMI and cell by concatenating group values in sorted name order
        let mut umi = Vec::new();
        let mut umi_quality = Vec::new();
        for &(_, start, end) in &umi_spans {
            umi.extend_from_slice(&sequence[start..end]);
            umi_quality.extend_from_slice(&quality[start..end]);
        }

        let mut cell_barcode = Vec::new();
        for &(_, start, end) in &cell_spans {
            cell_barcode.extend_from_slice(&sequence[start..end]);
        }

        // Build trimmed sequence/quality: keep positions not in any extraction set
        let mut trimmed_sequence = Vec::new();
        let mut trimmed_quality = Vec::new();

        for (i, &is_extracted) in extracted.iter().enumerate() {
            if !is_extracted {
                trimmed_sequence.push(sequence[i]);
                trimmed_quality.push(quality[i]);
            }
        }

        Ok(ExtractionResult {
            umi,
            umi_quality,
            cell_barcode,
            trimmed_sequence,
            trimmed_quality,
        })
    }
}

/// Pre-process a regex string, replacing `CHAR{s<=N}` fuzzy quantifiers.
///
/// In Python's `regex` module, `{s<=N}` applies to the single preceding character
/// (not to an entire literal sequence). For N >= 1, `CHAR{s<=N}` matches any
/// single character, equivalent to `.`. For N == 0, it's an exact match (no-op).
fn preprocess_fuzzy(pattern_str: &str) -> Result<String, ExtractError> {
    let mut result = String::with_capacity(pattern_str.len());
    let bytes = pattern_str.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    while i < len {
        if bytes[i] == b'{'
            && i + 4 < len
            && bytes[i + 1] == b's'
            && bytes[i + 2] == b'<'
            && bytes[i + 3] == b'='
        {
            let num_start = i + 4;
            let mut num_end = num_start;
            while num_end < len && bytes[num_end].is_ascii_digit() {
                num_end += 1;
            }
            if num_end == num_start || num_end >= len || bytes[num_end] != b'}' {
                return Err(ExtractError::InvalidPattern(format!(
                    "malformed fuzzy quantifier at position {i}"
                )));
            }
            let max_subs: usize = std::str::from_utf8(&bytes[num_start..num_end])
                .expect("ASCII digits validated above")
                .parse()
                .expect("ASCII digits validated above");

            if result.is_empty() {
                return Err(ExtractError::InvalidPattern(format!(
                    "fuzzy quantifier at position {i} has no preceding character"
                )));
            }

            if max_subs >= 1 {
                // Replace the preceding character with '.' (any character)
                result.pop();
                result.push('.');
            }
            // For max_subs == 0, keep the character as-is (exact match)

            i = num_end + 1;
        } else {
            result.push(bytes[i] as char);
            i += 1;
        }
    }

    Ok(result)
}

/// If `positions` is a contiguous ascending sequence [a, a+1, ..., b-1], return Some(a..b).
fn as_contiguous_range(positions: &[usize]) -> Option<Range<usize>> {
    let start = *positions.first()?;
    let is_contiguous = positions
        .iter()
        .enumerate()
        .skip(1)
        .all(|(i, &pos)| pos == start + i);
    is_contiguous.then(|| start..start + positions.len())
}

fn extract_slice(source: &[u8], range: Option<&Range<usize>>, positions: &[usize]) -> Vec<u8> {
    range.map_or_else(
        || positions.iter().map(|&i| source[i]).collect(),
        |r| source[r.clone()].to_vec(),
    )
}

fn join_slices(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    result.extend_from_slice(a);
    result.extend_from_slice(b);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- StringPattern tests ---

    #[test]
    fn parse_valid_pattern() {
        let pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        assert_eq!(pat.umi_positions, vec![0, 1, 2, 7, 8]);
        assert_eq!(pat.sample_positions, vec![3, 4, 5, 6]);
        assert!(pat.cell_positions.is_empty());
        assert_eq!(pat.pattern_length, 9);
    }

    #[test]
    fn parse_pattern_with_cell() {
        let pat = StringPattern::parse("CCCNNNNXXXX", PrimeEnd::Five).unwrap();
        assert_eq!(pat.cell_positions, vec![0, 1, 2]);
        assert_eq!(pat.umi_positions, vec![3, 4, 5, 6]);
        assert_eq!(pat.sample_positions, vec![7, 8, 9, 10]);
    }

    #[test]
    fn parse_invalid_pattern() {
        assert!(StringPattern::parse("NNNZXXNN", PrimeEnd::Five).is_err());
        assert!(StringPattern::parse("", PrimeEnd::Five).is_err());
    }

    #[test]
    fn extract_5prime_nnnxxxxnn() {
        let pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        let seq = b"CAGGTTCAATCTCGGTGGGACCTC";
        let qual = b"1=DFFFFHHHHHJJJFGIJIJJIJ";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.umi, b"CAGAA");
        assert_eq!(result.umi_quality, b"1=DHH");
        assert!(result.cell_barcode.is_empty());
        assert_eq!(result.trimmed_sequence, b"GTTCTCTCGGTGGGACCTC");
        assert_eq!(result.trimmed_quality, b"FFFFHHHJJJFGIJIJJIJ");
    }

    #[test]
    fn extract_read_too_short() {
        let pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        assert!(pat.extract(b"ACGT", b"IIII").is_err());
    }

    #[test]
    fn extract_3prime() {
        let pat = StringPattern::parse("NNXX", PrimeEnd::Three).unwrap();
        let seq = b"ACGTAATTGG";
        let qual = b"IIIIIIIIII";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.umi, b"TT");
        assert_eq!(result.trimmed_sequence, b"ACGTAAGG");
    }

    // --- RegexPattern tests ---

    #[test]
    fn regex_parse_valid() {
        let pat = RegexPattern::parse(r"^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})").unwrap();
        assert!(pat.pattern.is_match("CAGGTTCAATCTCGGTGGGACCTC"));
    }

    #[test]
    fn regex_parse_no_barcode_groups() {
        assert!(RegexPattern::parse(r"^(.{3}).{4}(.{2})").is_err());
    }

    #[test]
    fn regex_parse_invalid_regex() {
        assert!(RegexPattern::parse(r"^(?P<umi_1>.{3").is_err());
    }

    #[test]
    fn regex_extract_equivalent_to_string() {
        // Regex ^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2}) should produce same result as NNNXXXXNN
        let string_pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        let regex_pat = RegexPattern::parse(r"^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})").unwrap();

        let seq = b"CAGGTTCAATCTCGGTGGGACCTC";
        let qual = b"1=DFFFFHHHHHJJJFGIJIJJIJ";

        let string_result = string_pat.extract(seq, qual).unwrap();
        let regex_result = regex_pat.extract(seq, qual).unwrap();

        assert_eq!(string_result.umi, regex_result.umi);
        assert_eq!(string_result.cell_barcode, regex_result.cell_barcode);
        assert_eq!(
            string_result.trimmed_sequence,
            regex_result.trimmed_sequence
        );
        assert_eq!(string_result.trimmed_quality, regex_result.trimmed_quality);
    }

    #[test]
    fn regex_extract_with_cell() {
        let pat =
            RegexPattern::parse(r"^(?P<cell_1>.{3})(?P<umi_1>.{4})(?P<discard_1>.{2})").unwrap();

        let seq = b"ABCDEFGHIJKLM";
        let qual = b"1234567890ABC";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.cell_barcode, b"ABC");
        assert_eq!(result.umi, b"DEFG");
        // Positions 0-8 extracted/discarded, remaining: JKLM (positions 9-12)
        assert_eq!(result.trimmed_sequence, b"JKLM");
        assert_eq!(result.trimmed_quality, b"0ABC");
    }

    #[test]
    fn regex_no_match() {
        let pat = RegexPattern::parse(r"^(?P<umi_1>ZZZZZ)").unwrap();
        let result = pat.extract(b"ACGTACGT", b"IIIIIIII");
        assert!(matches!(result, Err(ExtractError::RegexNoMatch)));
    }
}
