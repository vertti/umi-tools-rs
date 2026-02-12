use crate::error::ExtractError;

/// Result of extracting barcodes from a single read's sequence and quality.
#[derive(Debug, Clone)]
pub struct ExtractionResult {
    pub umi: Vec<u8>,
    pub cell_barcode: Vec<u8>,
    pub trimmed_sequence: Vec<u8>,
    pub trimmed_quality: Vec<u8>,
}

/// Which end of the read to extract the barcode from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimeEnd {
    Five,
    Three,
}

/// A parsed barcode pattern that knows how to extract UMI/cell/sample bases from a read.
///
/// Pattern characters:
/// - `N` — UMI base (extracted to read name)
/// - `C` — Cell barcode base (extracted to read name)
/// - `X` — Sample/discard base (stays in output sequence, removed from barcode region)
#[derive(Debug, Clone)]
pub struct BarcodePattern {
    umi_positions: Vec<usize>,
    cell_positions: Vec<usize>,
    sample_positions: Vec<usize>,
    pattern_length: usize,
    prime_end: PrimeEnd,
}

impl BarcodePattern {
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

        Ok(Self {
            umi_positions,
            cell_positions,
            sample_positions,
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

        let umi = gather_positions(barcode_region, &self.umi_positions);
        let cell_barcode = gather_positions(barcode_region, &self.cell_positions);
        let sample_seq = gather_positions(barcode_region, &self.sample_positions);
        let sample_qual = gather_positions(barcode_qual, &self.sample_positions);

        let (trimmed_sequence, trimmed_quality) = match self.prime_end {
            PrimeEnd::Five => (
                join_slices(&sample_seq, remaining_seq),
                join_slices(&sample_qual, remaining_qual),
            ),
            PrimeEnd::Three => (
                join_slices(remaining_seq, &sample_seq),
                join_slices(remaining_qual, &sample_qual),
            ),
        };

        Ok(ExtractionResult {
            umi,
            cell_barcode,
            trimmed_sequence,
            trimmed_quality,
        })
    }
}

fn gather_positions(source: &[u8], positions: &[usize]) -> Vec<u8> {
    positions.iter().map(|&i| source[i]).collect()
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

    #[test]
    fn parse_valid_pattern() {
        let pat = BarcodePattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        assert_eq!(pat.umi_positions, vec![0, 1, 2, 7, 8]);
        assert_eq!(pat.sample_positions, vec![3, 4, 5, 6]);
        assert!(pat.cell_positions.is_empty());
        assert_eq!(pat.pattern_length, 9);
    }

    #[test]
    fn parse_pattern_with_cell() {
        let pat = BarcodePattern::parse("CCCNNNNXXXX", PrimeEnd::Five).unwrap();
        assert_eq!(pat.cell_positions, vec![0, 1, 2]);
        assert_eq!(pat.umi_positions, vec![3, 4, 5, 6]);
        assert_eq!(pat.sample_positions, vec![7, 8, 9, 10]);
    }

    #[test]
    fn parse_invalid_pattern() {
        assert!(BarcodePattern::parse("NNNZXXNN", PrimeEnd::Five).is_err());
        assert!(BarcodePattern::parse("", PrimeEnd::Five).is_err());
    }

    #[test]
    fn extract_5prime_nnnxxxxnn() {
        let pat = BarcodePattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        let seq = b"CAGGTTCAATCTCGGTGGGACCTC";
        let qual = b"1=DFFFFHHHHHJJJFGIJIJJIJ";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.umi, b"CAGAA");
        assert!(result.cell_barcode.is_empty());
        assert_eq!(result.trimmed_sequence, b"GTTCTCTCGGTGGGACCTC");
        assert_eq!(result.trimmed_quality, b"FFFFHHHJJJFGIJIJJIJ");
    }

    #[test]
    fn extract_read_too_short() {
        let pat = BarcodePattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        assert!(pat.extract(b"ACGT", b"IIII").is_err());
    }

    #[test]
    fn extract_3prime() {
        let pat = BarcodePattern::parse("NNXX", PrimeEnd::Three).unwrap();
        let seq = b"ACGTAATTGG";
        let qual = b"IIIIIIIIII";

        let result = pat.extract(seq, qual).unwrap();

        // Last 4 bases: TTGG
        // UMI positions 0,1 → TT
        // Sample positions 2,3 → GG
        // Remaining: ACGTAA
        // Trimmed = remaining + sample = ACGTAAGG
        assert_eq!(result.umi, b"TT");
        assert_eq!(result.trimmed_sequence, b"ACGTAAGG");
    }
}
