use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};

use needletail::parser::{FastqReader, FastxReader, SequenceRecord};

use crate::error::ExtractError;
use crate::pattern::{BarcodePattern, ExtractionResult};

/// Returns `true` if any base in `umi_quality` falls below the threshold after
/// subtracting the encoding offset.
fn fails_quality_filter(umi_quality: &[u8], threshold: u8, offset: u8) -> bool {
    umi_quality
        .iter()
        .any(|&q| q.saturating_sub(offset) < threshold)
}

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
    pub quality_filter_mask: Option<u8>,
    pub either_read_resolve: EitherReadResolve,
    pub subset_reads: Option<u64>,
}

/// What to do when both reads of a pair match in `--either-read` mode.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum EitherReadResolve {
    #[default]
    Discard,
    /// Keep the read whose UMI has the higher minimum base quality; read1 wins ties.
    Quality,
}

impl ExtractConfig {
    /// UMI as written to the read name: with `--quality-filter-mask`, bases below the
    /// threshold become `N`.
    fn final_umi(&self, umi: &[u8], quality: &[u8]) -> Vec<u8> {
        let Some(mask) = self.quality_filter_mask else {
            return umi.to_vec();
        };
        let offset = self.quality_encoding.offset();
        umi.iter()
            .zip(quality)
            .map(|(&base, &q)| {
                if q.saturating_sub(offset) < mask {
                    b'N'
                } else {
                    base
                }
            })
            .collect()
    }

    /// `--subset-reads`: `umi_tools` stops once the input count exceeds the limit.
    fn past_subset(&self, input_reads: u64) -> bool {
        self.subset_reads.is_some_and(|limit| input_reads > limit)
    }
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
    pub both_matched: u64,
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
    let mut writer = BufWriter::with_capacity(64 * 1024, output);
    let mut reader = FastqReader::new(input);

    while let Some(result) = reader.next() {
        let record = result.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
        stats.input_reads += 1;
        if config.past_subset(stats.input_reads) {
            break;
        }

        match process_record(&record, pattern) {
            Ok(processed) => {
                if let Some(threshold) = config.quality_filter_threshold
                    && fails_quality_filter(
                        &processed.umi_quality,
                        threshold,
                        config.quality_encoding.offset(),
                    )
                {
                    stats.quality_filtered += 1;
                    continue;
                }
                let umi = config.final_umi(&processed.umi, &processed.umi_quality);
                let id = build_read_name(
                    record.id(),
                    &processed.cell,
                    &umi,
                    config.umi_separator,
                    false,
                );
                write_fastq_record(&mut writer, &id, &processed.seq, &processed.qual)?;
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
    cell: Vec<u8>,
    umi: Vec<u8>,
    umi_quality: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

fn process_record(
    record: &SequenceRecord,
    pattern: &BarcodePattern,
) -> Result<ProcessedRecord, ExtractError> {
    let seq = record.seq();
    let qual = record
        .qual()
        .ok_or_else(|| ExtractError::FastqParse("missing quality scores in FASTQ record".into()))?;

    let result = pattern.extract(&seq, qual)?;

    Ok(ProcessedRecord {
        cell: result.cell_barcode,
        umi: result.umi,
        umi_quality: result.umi_quality,
        seq: result.trimmed_sequence,
        qual: result.trimmed_quality,
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

/// Extract UMIs from paired-end FASTQ reads with a read2 pattern.
///
/// When a read1 pattern is also supplied, concatenate barcodes in read1/read2
/// order and trim both reads. Otherwise, leave read1's sequence unchanged.
/// Append the combined barcodes to both read names.
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
    let mut writer1 = BufWriter::with_capacity(64 * 1024, output1);
    let mut writer2 = BufWriter::with_capacity(64 * 1024, output2);
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
                if config.past_subset(stats.input_reads) {
                    break;
                }

                let (mut extraction1, extraction) =
                    match process_paired_records(&r1, &r2, config.pattern.as_ref(), pattern2) {
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

                let (umi, umi_quality, cell_barcode) = if let Some(first) = extraction1.as_mut() {
                    first.umi.extend_from_slice(&extraction.umi);
                    first.umi_quality.extend_from_slice(&extraction.umi_quality);
                    first.cell.extend_from_slice(&extraction.cell);
                    (&first.umi, &first.umi_quality, &first.cell)
                } else {
                    (&extraction.umi, &extraction.umi_quality, &extraction.cell)
                };

                if let Some(threshold) = config.quality_filter_threshold
                    && fails_quality_filter(
                        umi_quality,
                        threshold,
                        config.quality_encoding.offset(),
                    )
                {
                    stats.quality_filtered += 1;
                    continue;
                }

                let umi = config.final_umi(umi, umi_quality);
                let r1_id = build_read_name(
                    r1.id(),
                    cell_barcode,
                    &umi,
                    config.umi_separator,
                    config.ignore_read_pair_suffixes,
                );
                let r2_id = build_read_name(
                    r2.id(),
                    cell_barcode,
                    &umi,
                    config.umi_separator,
                    config.ignore_read_pair_suffixes,
                );

                if let Some(first) = extraction1 {
                    write_fastq_record(&mut writer1, &r1_id, &first.seq, &first.qual)?;
                } else {
                    let r1_seq = r1.seq();
                    let r1_qual = r1.qual().ok_or_else(|| {
                        ExtractError::FastqParse("missing quality scores in read1".into())
                    })?;
                    write_fastq_record(&mut writer1, &r1_id, &r1_seq, r1_qual)?;
                }

                // Read2: trimmed, with new read name
                write_fastq_record(&mut writer2, &r2_id, &extraction.seq, &extraction.qual)?;

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

fn process_paired_records(
    r1: &SequenceRecord,
    r2: &SequenceRecord,
    pattern1: Option<&BarcodePattern>,
    pattern2: &BarcodePattern,
) -> Result<(Option<ProcessedRecord>, ProcessedRecord), ExtractError> {
    let first = pattern1
        .map(|pattern| process_record(r1, pattern))
        .transpose()?;
    let second = process_record(r2, pattern2)?;
    Ok((first, second))
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

    if let Some(threshold) = config.quality_filter_threshold
        && fails_quality_filter(
            &extraction.umi_quality,
            threshold,
            config.quality_encoding.offset(),
        )
    {
        stats.quality_filtered += 1;
        return Ok(false);
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

    let umi = config.final_umi(&extraction.umi, &extraction.umi_quality);
    let r2_id = build_read_name(
        r2.id(),
        &cell_barcode,
        &umi,
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
    let mut writer = BufWriter::with_capacity(64 * 1024, output);
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
            if config.past_subset(stats.input_reads) {
                break;
            }
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
                    if config.past_subset(stats.input_reads) {
                        break;
                    }
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

/// Try extracting with a pattern, returning `None` for recoverable failures (too short, no match).
fn try_extract(
    pattern: &BarcodePattern,
    seq: &[u8],
    qual: &[u8],
) -> Result<Option<ExtractionResult>, ExtractError> {
    match pattern.extract(seq, qual) {
        Ok(result) => Ok(Some(result)),
        Err(ExtractError::ReadTooShort { .. } | ExtractError::RegexNoMatch) => Ok(None),
        Err(e) => Err(e),
    }
}

/// Extract UMIs from paired-end FASTQ reads in either-read mode.
///
/// Both patterns are tried on their respective reads. If exactly one matches,
/// the UMI is taken from that read (and only that read is trimmed). If both
/// match, the pair is discarded (default `--either-read-resolve=discard`).
/// If neither matches, the pair is discarded as `no_match`.
///
/// # Errors
/// Returns error on I/O failures, parse errors, or mismatched read counts.
#[allow(clippy::too_many_lines)]
pub fn extract_reads_either_read<R1, R2, W1, W2>(
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
    let pattern1 = config.pattern.as_ref().ok_or_else(|| {
        ExtractError::InvalidPattern("no pattern provided for either-read extraction".into())
    })?;
    let pattern2 = config.pattern2.as_ref().ok_or_else(|| {
        ExtractError::InvalidPattern("no pattern2 provided for either-read extraction".into())
    })?;

    let mut stats = ExtractStats::default();
    let mut writer1 = BufWriter::with_capacity(64 * 1024, output1);
    let mut writer2 = BufWriter::with_capacity(64 * 1024, output2);
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
                if config.past_subset(stats.input_reads) {
                    break;
                }

                let r1_seq = r1.seq();
                let r1_qual = r1.qual().ok_or_else(|| {
                    ExtractError::FastqParse("missing quality scores in read1".into())
                })?;
                let r2_seq = r2.seq();
                let r2_qual = r2.qual().ok_or_else(|| {
                    ExtractError::FastqParse("missing quality scores in read2".into())
                })?;

                let r1_result = try_extract(pattern1, &r1_seq, r1_qual)?;
                let r2_result = try_extract(pattern2, &r2_seq, r2_qual)?;

                match (r1_result, r2_result) {
                    (Some(extraction1), Some(extraction2)) => {
                        if config.either_read_resolve == EitherReadResolve::Discard {
                            stats.both_matched += 1;
                            continue;
                        }
                        let offset = config.quality_encoding.offset();
                        let min_quality = |quality: &[u8]| {
                            quality
                                .iter()
                                .map(|&q| q.saturating_sub(offset))
                                .min()
                                .unwrap_or(0)
                        };
                        let chosen = if min_quality(&extraction1.umi_quality)
                            >= min_quality(&extraction2.umi_quality)
                        {
                            extraction1
                        } else {
                            extraction2
                        };
                        if let Some(threshold) = config.quality_filter_threshold
                            && fails_quality_filter(&chosen.umi_quality, threshold, offset)
                        {
                            stats.quality_filtered += 1;
                            continue;
                        }
                        let umi = config.final_umi(&chosen.umi, &chosen.umi_quality);
                        let new_id = build_read_name(
                            r1.id(),
                            &chosen.cell_barcode,
                            &umi,
                            config.umi_separator,
                            false,
                        );
                        // umi_tools leaves both sequences untrimmed when both reads match.
                        write_fastq_record(&mut writer1, &new_id, &r1_seq, r1_qual)?;
                        write_fastq_record(&mut writer2, &new_id, &r2_seq, r2_qual)?;
                        stats.output_reads += 1;
                    }
                    (Some(extraction), None) => {
                        if let Some(threshold) = config.quality_filter_threshold
                            && fails_quality_filter(
                                &extraction.umi_quality,
                                threshold,
                                config.quality_encoding.offset(),
                            )
                        {
                            stats.quality_filtered += 1;
                            continue;
                        }

                        let umi = config.final_umi(&extraction.umi, &extraction.umi_quality);
                        // Both headers built from read1 (matches Python umi-tools behavior)
                        let new_id = build_read_name(
                            r1.id(),
                            &extraction.cell_barcode,
                            &umi,
                            config.umi_separator,
                            false,
                        );

                        // Read1: trimmed
                        write_fastq_record(
                            &mut writer1,
                            &new_id,
                            &extraction.trimmed_sequence,
                            &extraction.trimmed_quality,
                        )?;
                        // Read2: untrimmed
                        write_fastq_record(&mut writer2, &new_id, &r2_seq, r2_qual)?;

                        stats.output_reads += 1;
                    }
                    (None, Some(extraction)) => {
                        if let Some(threshold) = config.quality_filter_threshold
                            && fails_quality_filter(
                                &extraction.umi_quality,
                                threshold,
                                config.quality_encoding.offset(),
                            )
                        {
                            stats.quality_filtered += 1;
                            continue;
                        }

                        let umi = config.final_umi(&extraction.umi, &extraction.umi_quality);
                        // Both headers built from read1 (matches Python umi-tools behavior)
                        let new_id = build_read_name(
                            r1.id(),
                            &extraction.cell_barcode,
                            &umi,
                            config.umi_separator,
                            false,
                        );

                        // Read1: untrimmed
                        write_fastq_record(&mut writer1, &new_id, &r1_seq, r1_qual)?;
                        // Read2: trimmed
                        write_fastq_record(
                            &mut writer2,
                            &new_id,
                            &extraction.trimmed_sequence,
                            &extraction.trimmed_quality,
                        )?;

                        stats.output_reads += 1;
                    }
                    (None, None) => {
                        stats.no_match += 1;
                    }
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

    writer1.flush()?;
    writer2.flush()?;
    Ok(stats)
}
