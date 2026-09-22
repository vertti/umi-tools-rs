use std::borrow::Cow;
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

/// Borrow unchanged slices while the input record remains in the parser buffer.
pub(crate) struct ExtractionView<'a> {
    pub umi: Cow<'a, [u8]>,
    pub umi_quality: Cow<'a, [u8]>,
    pub cell_barcode: Cow<'a, [u8]>,
    pub trimmed_sequence: Cow<'a, [u8]>,
    pub trimmed_quality: Cow<'a, [u8]>,
}

impl From<ExtractionResult> for ExtractionView<'_> {
    fn from(result: ExtractionResult) -> Self {
        Self {
            umi: Cow::Owned(result.umi),
            umi_quality: Cow::Owned(result.umi_quality),
            cell_barcode: Cow::Owned(result.cell_barcode),
            trimmed_sequence: Cow::Owned(result.trimmed_sequence),
            trimmed_quality: Cow::Owned(result.trimmed_quality),
        }
    }
}

impl ExtractionView<'_> {
    fn into_owned(self) -> ExtractionResult {
        ExtractionResult {
            umi: self.umi.into_owned(),
            umi_quality: self.umi_quality.into_owned(),
            cell_barcode: self.cell_barcode.into_owned(),
            trimmed_sequence: self.trimmed_sequence.into_owned(),
            trimmed_quality: self.trimmed_quality.into_owned(),
        }
    }
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
    pub(crate) fn extract_view<'a>(
        &self,
        sequence: &'a [u8],
        quality: &'a [u8],
    ) -> Result<ExtractionView<'a>, ExtractError> {
        match self {
            Self::String(pattern) => pattern.extract_view(sequence, quality),
            Self::Regex(pattern) => pattern.extract(sequence, quality).map(ExtractionView::from),
        }
    }

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
        self.extract_view(sequence, quality)
            .map(ExtractionView::into_owned)
    }

    fn extract_view<'a>(
        &self,
        sequence: &'a [u8],
        quality: &'a [u8],
    ) -> Result<ExtractionView<'a>, ExtractError> {
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
            (Cow::Borrowed(remaining_seq), Cow::Borrowed(remaining_qual))
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
                    Cow::Owned(join_slices(&sample_seq, remaining_seq)),
                    Cow::Owned(join_slices(&sample_qual, remaining_qual)),
                ),
                PrimeEnd::Three => (
                    Cow::Owned(join_slices(remaining_seq, &sample_seq)),
                    Cow::Owned(join_slices(remaining_qual, &sample_qual)),
                ),
            }
        };

        Ok(ExtractionView {
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
            .filter(|captures| captures.get(0).is_some_and(|matched| matched.start() == 0))
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
                .map_err(|_| {
                    ExtractError::InvalidPattern(format!(
                        "fuzzy quantifier at position {i} exceeds the supported integer range"
                    ))
                })?;

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

fn extract_slice<'a>(
    source: &'a [u8],
    range: Option<&Range<usize>>,
    positions: &[usize],
) -> Cow<'a, [u8]> {
    match range {
        Some(range) => Cow::Borrowed(&source[range.clone()]),
        None if positions.is_empty() => Cow::Borrowed(&[]),
        None => Cow::Owned(positions.iter().map(|&i| source[i]).collect()),
    }
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
    fn contiguous_extraction_borrows_all_five_fields() {
        let sequence = b"ACGTACGT";
        let quality = b"12345678";
        for (end, cell, umi, remaining) in [
            (PrimeEnd::Five, 0..2, 2..4, 4..8),
            (PrimeEnd::Three, 4..6, 6..8, 0..4),
        ] {
            let pattern = StringPattern::parse("CCNN", end).unwrap();
            let view = pattern.extract_view(sequence, quality).unwrap();
            for (actual, expected) in [
                (&view.cell_barcode, &sequence[cell]),
                (&view.umi, &sequence[umi.clone()]),
                (&view.umi_quality, &quality[umi]),
                (&view.trimmed_sequence, &sequence[remaining.clone()]),
                (&view.trimmed_quality, &quality[remaining]),
            ] {
                assert_eq!(actual.as_ref(), expected);
                assert!(
                    matches!(actual, Cow::Borrowed(_)),
                    "copied an unchanged slice"
                );
                assert_eq!(actual.as_ptr(), expected.as_ptr());
            }
        }
    }

    #[test]
    fn extraction_matches_position_reference_for_all_short_patterns() {
        for length in 1..=4u32 {
            for mut code in 0..3usize.pow(length) {
                let pattern: String = (0..length)
                    .map(|_| {
                        let ch = char::from(b"NCX"[code % 3]);
                        code /= 3;
                        ch
                    })
                    .collect();
                for end in [PrimeEnd::Five, PrimeEnd::Three] {
                    for extra in [0, 3] {
                        let size = pattern.len() + extra;
                        let sequence = &b"ACGTNAC"[..size];
                        let quality = &b"1234567"[..size];
                        let start = if end == PrimeEnd::Five { 0 } else { extra };
                        let marker = |index: usize| {
                            index
                                .checked_sub(start)
                                .and_then(|offset| pattern.as_bytes().get(offset))
                                .copied()
                                .unwrap_or(b'X')
                        };
                        let select = |source: &[u8], kind| -> Vec<u8> {
                            source
                                .iter()
                                .enumerate()
                                .filter(|&(index, _)| marker(index) == kind)
                                .map(|(_, &base)| base)
                                .collect()
                        };
                        let result = StringPattern::parse(&pattern, end)
                            .unwrap()
                            .extract_view(sequence, quality)
                            .unwrap();
                        assert_eq!(result.umi.as_ref(), select(sequence, b'N'));
                        assert_eq!(result.umi_quality.as_ref(), select(quality, b'N'));
                        assert_eq!(result.cell_barcode.as_ref(), select(sequence, b'C'));
                        assert_eq!(result.trimmed_sequence.as_ref(), select(sequence, b'X'));
                        assert_eq!(result.trimmed_quality.as_ref(), select(quality, b'X'));
                    }
                }
            }
        }
    }

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
    fn regex_matching_starts_at_the_beginning_of_the_read() {
        let pattern = RegexPattern::parse("(?P<umi_1>AA)").unwrap();
        assert!(matches!(
            pattern.extract(b"TTAA", b"IIII"),
            Err(ExtractError::RegexNoMatch)
        ));
        assert_eq!(pattern.extract(b"AATT", b"IIII").unwrap().umi, b"AA");

        let pattern = RegexPattern::parse(".*(?P<umi_1>AA)").unwrap();
        assert_eq!(pattern.extract(b"TTAA", b"IIII").unwrap().umi, b"AA");
    }

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
    fn oversized_fuzzy_quantifier_returns_an_error() {
        assert!(matches!(
            RegexPattern::parse("(?P<umi_1>A{s<=9999999999999999999999999999999})"),
            Err(ExtractError::InvalidPattern(_))
        ));
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
