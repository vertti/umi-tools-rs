use std::collections::{HashMap, HashSet};
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
    pub correction_map: Option<HashMap<Vec<u8>, Vec<u8>>>,
    pub blacklist: Option<HashSet<Vec<u8>>>,
    pub ignore_read_pair_suffixes: bool,
    pub reconcile_pairs: bool,
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
        false,
    );

    Ok(ProcessedRecord {
        id,
        seq: result.trimmed_sequence,
        qual: result.trimmed_quality,
        umi_quality: result.umi_quality,
    })
}

/// Return the name portion of a FASTQ header (before first space).
fn read_name(header: &[u8]) -> &[u8] {
    header
        .iter()
        .position(|&b| b == b' ')
        .map_or(header, |pos| &header[..pos])
}

/// Strip trailing `/1` or `/2` from a read name.
fn strip_pair_suffix(name: &[u8]) -> &[u8] {
    if name.len() >= 2 && name[name.len() - 2] == b'/' {
        let last = name[name.len() - 1];
        if last == b'1' || last == b'2' {
            return &name[..name.len() - 2];
        }
    }
    name
}

/// Build a new read identifier with barcode(s) inserted after the read name.
///
/// Splits header at first space: `NAME COMMENT` → `NAME{sep}UMI COMMENT`.
/// If `strip_suffixes` is true, strips trailing `/1` or `/2` from the name.
fn build_read_name(
    header: &[u8],
    cell: &[u8],
    umi: &[u8],
    separator: u8,
    strip_suffixes: bool,
) -> Vec<u8> {
    let (name, comment) = header
        .iter()
        .position(|&b| b == b' ')
        .map_or((header, None), |pos| (&header[..pos], Some(&header[pos..])));

    let name = if strip_suffixes {
        strip_pair_suffix(name)
    } else {
        name
    };

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
                    false,
                );
                let r2_id = build_read_name(
                    r2.id(),
                    &extraction.cell_barcode,
                    &extraction.umi,
                    config.umi_separator,
                    false,
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

/// Process a single read pair in the r1-pattern extraction path.
///
/// Returns `true` if the pair produced output, `false` if filtered/skipped.
fn process_r1_pattern_pair<W: Write>(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    pattern: &BarcodePattern,
    config: &ExtractConfig,
    stats: &mut ExtractStats,
    writer: &mut W,
) -> Result<bool, ExtractError> {
    let r1_seq = r1.seq();
    let r1_qual = r1
        .qual()
        .ok_or_else(|| ExtractError::FastqParse("missing quality scores in read1".into()))?;

    let extraction = match pattern.extract(&r1_seq, r1_qual) {
        Ok(result) => result,
        Err(ExtractError::ReadTooShort { .. }) => {
            stats.too_short += 1;
            return Ok(false);
        }
        Err(ExtractError::RegexNoMatch) => {
            stats.no_match += 1;
            return Ok(false);
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
            return Ok(false);
        }
    }

    if let Some(ref blacklist) = config.blacklist
        && blacklist.contains(&extraction.cell_barcode)
    {
        stats.whitelist_filtered += 1;
        return Ok(false);
    }

    let cell_barcode = if let Some(ref whitelist) = config.whitelist {
        if whitelist.contains(&extraction.cell_barcode) {
            extraction.cell_barcode.clone()
        } else if let Some(ref correction_map) = config.correction_map
            && let Some(corrected) = correction_map.get(&extraction.cell_barcode)
        {
            corrected.clone()
        } else {
            stats.whitelist_filtered += 1;
            return Ok(false);
        }
    } else {
        extraction.cell_barcode.clone()
    };

    let r2_id = build_read_name(
        r2.id(),
        &cell_barcode,
        &extraction.umi,
        config.umi_separator,
        config.ignore_read_pair_suffixes,
    );

    let r2_seq = r2.seq();
    let r2_qual = r2
        .qual()
        .ok_or_else(|| ExtractError::FastqParse("missing quality scores in read2".into()))?;
    write_fastq_record(writer, &r2_id, &r2_seq, r2_qual)?;

    stats.output_reads += 1;
    Ok(true)
}

/// Write original untrimmed reads to filtered output files (headers NOT modified).
fn write_filtered_pair<W: Write>(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    filt_writer1: &mut Option<BufWriter<W>>,
    filt_writer2: &mut Option<BufWriter<W>>,
) -> Result<(), ExtractError> {
    if let Some(fw) = filt_writer1.as_mut() {
        let r1_qual = r1
            .qual()
            .ok_or_else(|| ExtractError::FastqParse("missing quality scores in read1".into()))?;
        write_fastq_record(fw, r1.id(), &r1.seq(), r1_qual)?;
    }
    if let Some(fw) = filt_writer2.as_mut() {
        let r2_qual = r2
            .qual()
            .ok_or_else(|| ExtractError::FastqParse("missing quality scores in read2".into()))?;
        write_fastq_record(fw, r2.id(), &r2.seq(), r2_qual)?;
    }
    Ok(())
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
    filtered_out1: Option<Box<dyn Write>>,
    filtered_out2: Option<Box<dyn Write>>,
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
    let mut filt_writer1 = filtered_out1.map(BufWriter::new);
    let mut filt_writer2 = filtered_out2.map(BufWriter::new);
    let mut reader1 = FastqReader::new(input1);
    let mut reader2 = FastqReader::new(input2);

    if config.reconcile_pairs {
        // Reconcile mode: read1 is a pre-filtered subset, read2 is the full set.
        // Both files maintain original sequencing order. For each read1 record,
        // advance read2 until a matching read name is found; skip unmatched read2s.
        while let Some(r1_result) = reader1.next() {
            let r1 = r1_result.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
            stats.input_reads += 1;
            let r1_name = read_name(r1.id());

            loop {
                match reader2.next() {
                    Some(r2_result) => {
                        let r2 = r2_result.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                        if read_name(r2.id()) == r1_name {
                            let kept = process_r1_pattern_pair(
                                &r1,
                                &r2,
                                pattern,
                                config,
                                &mut stats,
                                &mut writer,
                            )?;
                            if !kept {
                                write_filtered_pair(
                                    &r1,
                                    &r2,
                                    &mut filt_writer1,
                                    &mut filt_writer2,
                                )?;
                            }
                            break;
                        }
                    }
                    None => {
                        return Err(ExtractError::FastqParse(format!(
                            "read2 exhausted before finding match for read1: {}",
                            String::from_utf8_lossy(r1_name)
                        )));
                    }
                }
            }
        }
    } else {
        // Lockstep mode: read1 and read2 must have matching records in order.
        loop {
            let rec1 = reader1.next();
            let rec2 = reader2.next();

            match (rec1, rec2) {
                (Some(r1), Some(r2)) => {
                    let r1 = r1.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                    let r2 = r2.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
                    stats.input_reads += 1;
                    let kept = process_r1_pattern_pair(
                        &r1,
                        &r2,
                        pattern,
                        config,
                        &mut stats,
                        &mut writer,
                    )?;
                    if !kept {
                        write_filtered_pair(&r1, &r2, &mut filt_writer1, &mut filt_writer2)?;
                    }
                }
                (None, None) => break,
                _ => {
                    return Err(ExtractError::FastqParse(
                        "read1 and read2 files have different numbers of records".into(),
                    ));
                }
            }
        }
    }

    writer.flush()?;
    if let Some(fw) = filt_writer1.as_mut() {
        fw.flush()?;
    }
    if let Some(fw) = filt_writer2.as_mut() {
        fw.flush()?;
    }
    Ok(stats)
}
