use std::collections::HashSet;
use std::io::{BufWriter, Write};

use needletail::parser::{FastqReader, FastxReader, SequenceRecord};

use crate::error::ExtractError;
use crate::pattern::BarcodePattern;

/// Quality score encoding scheme.
#[derive(Debug, Clone, Copy, Default)]
pub enum QualityEncoding {
    #[default]
    Phred33,
    Phred64,
    Solexa,
}

impl QualityEncoding {
    #[must_use]
    pub const fn offset(self) -> u8 {
        match self {
            Self::Phred33 => 33,
            Self::Phred64 => 64,
            Self::Solexa => 59,
        }
    }
}

/// Configuration for the extract command.
#[derive(Debug, Clone)]
pub struct ExtractConfig {
    pub pattern: Option<BarcodePattern>,
    pub pattern2: Option<BarcodePattern>,
    pub umi_separator: u8,
    pub quality_filter_threshold: Option<u8>,
    pub quality_encoding: QualityEncoding,
    pub whitelist: Option<HashSet<Vec<u8>>>,
}

/// Statistics from an extraction run.
#[derive(Debug, Default)]
pub struct ExtractStats {
    pub input_reads: u64,
    pub output_reads: u64,
    pub too_short: u64,
    pub no_match: u64,
    pub quality_filtered: u64,
    pub whitelist_filtered: u64,
}

/// Extract UMIs from FASTQ reads, writing modified reads to `output`.
///
/// # Errors
/// Returns error on I/O or parse failures.
pub fn extract_reads<R: std::io::Read + Send, W: Write>(
    config: &ExtractConfig,
    input: R,
    output: W,
) -> Result<ExtractStats, ExtractError> {
    let pattern = config.pattern.as_ref().ok_or_else(|| {
        ExtractError::InvalidPattern("no pattern provided for single-end extraction".into())
    })?;

    let mut stats = ExtractStats::default();
    let mut writer = BufWriter::new(output);
    let mut reader = FastqReader::new(input);

    while let Some(result) = reader.next() {
        let record = result.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
        stats.input_reads += 1;

        match process_record(&record, pattern, config.umi_separator) {
            Ok(processed) => {
                if let Some(threshold) = config.quality_filter_threshold {
                    let offset = config.quality_encoding.offset();
                    if processed
                        .umi_quality
                        .iter()
                        .any(|&q| q.saturating_sub(offset) < threshold)
                    {
                        stats.quality_filtered += 1;
                        continue;
                    }
                }
                write_fastq_record(&mut writer, &processed.id, &processed.seq, &processed.qual)?;
                stats.output_reads += 1;
            }
            Err(ExtractError::ReadTooShort { .. }) => {
                stats.too_short += 1;
            }
            Err(ExtractError::RegexNoMatch) => {
                stats.no_match += 1;
            }
            Err(e) => return Err(e),
        }
    }

    writer.flush()?;
    Ok(stats)
}

struct ProcessedRecord {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    umi_quality: Vec<u8>,
}

fn process_record(
    record: &SequenceRecord,
    pattern: &BarcodePattern,
    umi_separator: u8,
) -> Result<ProcessedRecord, ExtractError> {
    let seq = record.seq();
    let qual = record
        .qual()
        .ok_or_else(|| ExtractError::FastqParse("missing quality scores in FASTQ record".into()))?;

    let result = pattern.extract(&seq, qual)?;

    let id = build_read_name(
        record.id(),
        &result.cell_barcode,
        &result.umi,
        umi_separator,
    );

    Ok(ProcessedRecord {
        id,
        seq: result.trimmed_sequence,
        qual: result.trimmed_quality,
        umi_quality: result.umi_quality,
    })
}

/// Build a new read identifier with barcode(s) inserted after the read name.
///
/// Splits header at first space: `NAME COMMENT` → `NAME{sep}UMI COMMENT`
fn build_read_name(header: &[u8], cell: &[u8], umi: &[u8], separator: u8) -> Vec<u8> {
    let (name, comment) = header
        .iter()
        .position(|&b| b == b' ')
        .map_or((header, None), |pos| (&header[..pos], Some(&header[pos..])));

    let mut out = Vec::with_capacity(header.len() + 1 + cell.len() + 1 + umi.len());
    out.extend_from_slice(name);
    if !cell.is_empty() {
        out.push(separator);
        out.extend_from_slice(cell);
    }
    out.push(separator);
    out.extend_from_slice(umi);
    if let Some(c) = comment {
        out.extend_from_slice(c);
    }
    out
}

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

/// Extract UMIs from paired-end FASTQ reads (read2-only pattern mode).
///
/// Pattern is applied to read2 only. UMI is appended to both read names.
/// Read1 is written untrimmed to `output1`, read2 is written trimmed to `output2`.
///
/// # Errors
/// Returns error on I/O failures, parse errors, or mismatched read counts.
pub fn extract_reads_paired<R1, R2, W1, W2>(
    config: &ExtractConfig,
    input1: R1,
    input2: R2,
    output1: W1,
    output2: W2,
) -> Result<ExtractStats, ExtractError>
where
    R1: std::io::Read + Send,
    R2: std::io::Read + Send,
    W1: Write,
    W2: Write,
{
    let pattern2 = config.pattern2.as_ref().ok_or_else(|| {
        ExtractError::InvalidPattern("no pattern2 provided for paired-end extraction".into())
    })?;

    let mut stats = ExtractStats::default();
    let mut writer1 = BufWriter::new(output1);
    let mut writer2 = BufWriter::new(output2);
    let mut reader1 = FastqReader::new(input1);
    let mut reader2 = FastqReader::new(input2);

    loop {
        let rec1 = reader1.next();
        let rec2 = reader2.next();

        match (rec1, rec2) {
            (Some(r1), Some(r2)) => {
                let r1 = r1.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                let r2 = r2.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                stats.input_reads += 1;

                let r2_seq = r2.seq();
                let r2_qual = r2.qual().ok_or_else(|| {
                    ExtractError::FastqParse("missing quality scores in read2".into())
                })?;

                let extraction = match pattern2.extract(&r2_seq, r2_qual) {
                    Ok(result) => result,
                    Err(ExtractError::ReadTooShort { .. }) => {
                        stats.too_short += 1;
                        continue;
                    }
                    Err(ExtractError::RegexNoMatch) => {
                        stats.no_match += 1;
                        continue;
                    }
                    Err(e) => return Err(e),
                };

                if let Some(threshold) = config.quality_filter_threshold {
                    let offset = config.quality_encoding.offset();
                    if extraction
                        .umi_quality
                        .iter()
                        .any(|&q| q.saturating_sub(offset) < threshold)
                    {
                        stats.quality_filtered += 1;
                        continue;
                    }
                }

                let r1_id = build_read_name(
                    r1.id(),
                    &extraction.cell_barcode,
                    &extraction.umi,
                    config.umi_separator,
                );
                let r2_id = build_read_name(
                    r2.id(),
                    &extraction.cell_barcode,
                    &extraction.umi,
                    config.umi_separator,
                );

                // Read1: untrimmed, with new read name
                let r1_seq = r1.seq();
                let r1_qual = r1.qual().ok_or_else(|| {
                    ExtractError::FastqParse("missing quality scores in read1".into())
                })?;
                write_fastq_record(&mut writer1, &r1_id, &r1_seq, r1_qual)?;

                // Read2: trimmed, with new read name
                write_fastq_record(
                    &mut writer2,
                    &r2_id,
                    &extraction.trimmed_sequence,
                    &extraction.trimmed_quality,
                )?;

                stats.output_reads += 1;
            }
            (None, None) => break,
            _ => {
                return Err(ExtractError::FastqParse(
                    "read1 and read2 files have different numbers of records".into(),
                ));
            }
        }
    }

    writer1.flush()?;
    writer2.flush()?;
    Ok(stats)
}

/// Extract UMIs from paired-end FASTQ reads (read1-pattern mode with read2 output).
///
/// Pattern is applied to read1 to extract cell barcode + UMI. Read2 is written
/// untrimmed to `output` with the cell+UMI appended to read2's header.
/// Reads whose cell barcode is not in the whitelist (if provided) are discarded.
///
/// # Errors
/// Returns error on I/O failures, parse errors, or mismatched read counts.
pub fn extract_reads_paired_r1_pattern<R1, R2, W>(
    config: &ExtractConfig,
    input1: R1,
    input2: R2,
    output: W,
) -> Result<ExtractStats, ExtractError>
where
    R1: std::io::Read + Send,
    R2: std::io::Read + Send,
    W: Write,
{
    let pattern = config.pattern.as_ref().ok_or_else(|| {
        ExtractError::InvalidPattern(
            "no pattern provided for paired-end read1-pattern extraction".into(),
        )
    })?;

    let mut stats = ExtractStats::default();
    let mut writer = BufWriter::new(output);
    let mut reader1 = FastqReader::new(input1);
    let mut reader2 = FastqReader::new(input2);

    loop {
        let rec1 = reader1.next();
        let rec2 = reader2.next();

        match (rec1, rec2) {
            (Some(r1), Some(r2)) => {
                let r1 = r1.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                let r2 = r2.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                stats.input_reads += 1;

                let r1_seq = r1.seq();
                let r1_qual = r1.qual().ok_or_else(|| {
                    ExtractError::FastqParse("missing quality scores in read1".into())
                })?;

                let extraction = match pattern.extract(&r1_seq, r1_qual) {
                    Ok(result) => result,
                    Err(ExtractError::ReadTooShort { .. }) => {
                        stats.too_short += 1;
                        continue;
                    }
                    Err(ExtractError::RegexNoMatch) => {
                        stats.no_match += 1;
                        continue;
                    }
                    Err(e) => return Err(e),
                };

                if let Some(threshold) = config.quality_filter_threshold {
                    let offset = config.quality_encoding.offset();
                    if extraction
                        .umi_quality
                        .iter()
                        .any(|&q| q.saturating_sub(offset) < threshold)
                    {
                        stats.quality_filtered += 1;
                        continue;
                    }
                }

                if let Some(ref whitelist) = config.whitelist
                    && !whitelist.contains(&extraction.cell_barcode)
                {
                    stats.whitelist_filtered += 1;
                    continue;
                }

                let r2_id = build_read_name(
                    r2.id(),
                    &extraction.cell_barcode,
                    &extraction.umi,
                    config.umi_separator,
                );

                let r2_seq = r2.seq();
                let r2_qual = r2.qual().ok_or_else(|| {
                    ExtractError::FastqParse("missing quality scores in read2".into())
                })?;
                write_fastq_record(&mut writer, &r2_id, &r2_seq, r2_qual)?;

                stats.output_reads += 1;
            }
            (None, None) => break,
            _ => {
                return Err(ExtractError::FastqParse(
                    "read1 and read2 files have different numbers of records".into(),
                ));
            }
        }
    }

    writer.flush()?;
    Ok(stats)
}
