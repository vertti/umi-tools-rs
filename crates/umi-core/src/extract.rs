use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};

use needletail::parser::{FastqReader, FastxReader, SequenceRecord};

use crate::error::ExtractError;
use crate::fastq::write_fastq_record;
use crate::pattern::{BarcodePattern, ExtractionView};

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
    pub umi_separator: Vec<u8>,
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
    fn mask_umi(&self, umi: &mut Cow<'_, [u8]>, quality: &[u8]) {
        let Some(mask) = self.quality_filter_mask else {
            return;
        };
        let offset = self.quality_encoding.offset();
        for (index, &q) in quality.iter().enumerate() {
            if q.saturating_sub(offset) < mask && umi[index] != b'N' {
                umi.to_mut()[index] = b'N';
            }
        }
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

/// Return the name portion of a FASTQ header (before first space).
fn read_name(header: &[u8]) -> &[u8] {
    header
        .iter()
        .position(|&b| b == b' ')
        .map_or(header, |pos| &header[..pos])
}

fn pair_name(header: &[u8], strip_suffixes: bool) -> &[u8] {
    let name = read_name(header);
    if strip_suffixes {
        strip_pair_suffix(name)
    } else {
        name
    }
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
    separator: &[u8],
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

    let mut out = Vec::with_capacity(header.len() + 2 * separator.len() + cell.len() + umi.len());
    out.extend_from_slice(name);
    if !cell.is_empty() {
        out.extend_from_slice(separator);
        out.extend_from_slice(cell);
    }
    out.extend_from_slice(separator);
    out.extend_from_slice(umi);
    if let Some(c) = comment {
        out.extend_from_slice(c);
    }
    out
}

/// Whether to combine barcodes from the supplied patterns or choose either read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtractMode {
    Combine,
    EitherRead,
}

/// Output destinations are independent of how barcodes are extracted.
#[derive(Default)]
pub struct ExtractOutputs<'a> {
    pub read1: Option<Box<dyn Write + 'a>>,
    pub read2: Option<Box<dyn Write + 'a>>,
    pub filtered1: Option<Box<dyn Write + 'a>>,
    pub filtered2: Option<Box<dyn Write + 'a>>,
}

impl<'a> ExtractOutputs<'a> {
    fn paired<W1: Write + 'a, W2: Write + 'a>(read1: W1, read2: W2) -> Self {
        Self {
            read1: Some(Box::new(read1)),
            read2: Some(Box::new(read2)),
            ..Self::default()
        }
    }
}

type Output<'a> = Option<BufWriter<Box<dyn Write + 'a>>>;

struct Writers<'a> {
    read1: Output<'a>,
    read2: Output<'a>,
    filtered1: Output<'a>,
    filtered2: Output<'a>,
}

impl<'a> From<ExtractOutputs<'a>> for Writers<'a> {
    fn from(outputs: ExtractOutputs<'a>) -> Self {
        let buffer = |writer| BufWriter::with_capacity(64 * 1024, writer);
        Self {
            read1: outputs.read1.map(buffer),
            read2: outputs.read2.map(buffer),
            filtered1: outputs.filtered1.map(buffer),
            filtered2: outputs.filtered2.map(buffer),
        }
    }
}

impl Writers<'_> {
    fn flush(&mut self) -> Result<(), ExtractError> {
        for writer in [
            &mut self.read1,
            &mut self.read2,
            &mut self.filtered1,
            &mut self.filtered2,
        ]
        .into_iter()
        .flatten()
        {
            writer.flush()?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct Barcodes<'a> {
    cell: Cow<'a, [u8]>,
    umi: Cow<'a, [u8]>,
    quality: Cow<'a, [u8]>,
}

struct TrimmedRead<'a> {
    sequence: Cow<'a, [u8]>,
    quality: Cow<'a, [u8]>,
}

struct Extraction<'a> {
    barcodes: Barcodes<'a>,
    read1: Option<TrimmedRead<'a>>,
    read2: Option<TrimmedRead<'a>>,
}

fn split_extraction(result: ExtractionView<'_>) -> (Barcodes<'_>, Option<TrimmedRead<'_>>) {
    (
        Barcodes {
            cell: result.cell_barcode,
            umi: result.umi,
            quality: result.umi_quality,
        },
        Some(TrimmedRead {
            sequence: result.trimmed_sequence,
            quality: result.trimmed_quality,
        }),
    )
}

fn extract_record<'a>(
    record: &'a SequenceRecord,
    pattern: &BarcodePattern,
) -> Result<ExtractionView<'a>, ExtractError> {
    let quality = record
        .qual()
        .ok_or_else(|| ExtractError::FastqParse("missing quality scores in FASTQ record".into()))?;
    pattern.extract_view(record.raw_seq(), quality)
}

fn append_barcode<'a>(first: &mut Cow<'a, [u8]>, second: Cow<'a, [u8]>) {
    if first.is_empty() {
        *first = second;
    } else if !second.is_empty() {
        first.to_mut().extend_from_slice(&second);
    }
}

fn extract_combined<'a>(
    config: &ExtractConfig,
    r1: &'a SequenceRecord,
    r2: Option<&'a SequenceRecord>,
) -> Result<Extraction<'a>, ExtractError> {
    let (mut barcodes, read1) = config
        .pattern
        .as_ref()
        .map(|pattern| extract_record(r1, pattern).map(split_extraction))
        .transpose()?
        .unwrap_or_default();
    let (second, read2) = config
        .pattern2
        .as_ref()
        .zip(r2)
        .map(|(pattern, record)| extract_record(record, pattern).map(split_extraction))
        .transpose()?
        .unwrap_or_default();
    append_barcode(&mut barcodes.cell, second.cell);
    append_barcode(&mut barcodes.umi, second.umi);
    append_barcode(&mut barcodes.quality, second.quality);
    Ok(Extraction {
        barcodes,
        read1,
        read2,
    })
}

fn try_extract<'a>(
    record: &'a SequenceRecord,
    pattern: &BarcodePattern,
) -> Result<Option<ExtractionView<'a>>, ExtractError> {
    match extract_record(record, pattern) {
        Ok(result) => Ok(Some(result)),
        Err(ExtractError::ReadTooShort { .. } | ExtractError::RegexNoMatch) => Ok(None),
        Err(e) => Err(e),
    }
}

fn extract_either<'a>(
    config: &ExtractConfig,
    r1: &'a SequenceRecord,
    r2: &'a SequenceRecord,
    stats: &mut ExtractStats,
) -> Result<Option<Extraction<'a>>, ExtractError> {
    let first = try_extract(
        r1,
        config.pattern.as_ref().expect("validated read1 pattern"),
    )?;
    let second = try_extract(
        r2,
        config.pattern2.as_ref().expect("validated read2 pattern"),
    )?;
    let extraction = match (first, second) {
        (Some(first), Some(second)) => {
            if config.either_read_resolve == EitherReadResolve::Discard {
                stats.both_matched += 1;
                return Ok(None);
            }
            let min_quality = |result: &ExtractionView| {
                result
                    .umi_quality
                    .iter()
                    .map(|q| q.saturating_sub(config.quality_encoding.offset()))
                    .min()
                    .unwrap_or(0)
            };
            let chosen = if min_quality(&first) >= min_quality(&second) {
                first
            } else {
                second
            };
            let (barcodes, _) = split_extraction(chosen);
            // Python leaves both sequences intact when both patterns match.
            Extraction {
                barcodes,
                read1: None,
                read2: None,
            }
        }
        (Some(first), None) => {
            let (barcodes, read1) = split_extraction(first);
            Extraction {
                barcodes,
                read1,
                read2: None,
            }
        }
        (None, Some(second)) => {
            let (barcodes, read2) = split_extraction(second);
            Extraction {
                barcodes,
                read1: None,
                read2,
            }
        }
        (None, None) => {
            stats.no_match += 1;
            return Ok(None);
        }
    };
    Ok(Some(extraction))
}

/// All extraction modes apply filters in the same order to the selected barcodes.
fn filter_barcodes(
    config: &ExtractConfig,
    barcodes: &mut Barcodes,
    stats: &mut ExtractStats,
) -> bool {
    if let Some(threshold) = config.quality_filter_threshold
        && fails_quality_filter(
            &barcodes.quality,
            threshold,
            config.quality_encoding.offset(),
        )
    {
        stats.quality_filtered += 1;
        return false;
    }
    config.mask_umi(&mut barcodes.umi, &barcodes.quality);
    let blacklisted = |cell: &[u8]| {
        config
            .blacklist
            .as_ref()
            .is_some_and(|list| list.contains(cell))
    };
    if blacklisted(&barcodes.cell) {
        stats.whitelist_filtered += 1;
        return false;
    }
    if let Some(whitelist) = &config.whitelist
        && !whitelist.contains(barcodes.cell.as_ref())
    {
        if let Some(corrected) = config
            .correction_map
            .as_ref()
            .and_then(|map| map.get(barcodes.cell.as_ref()))
        {
            barcodes.cell.to_mut().clone_from(corrected);
        } else {
            stats.whitelist_filtered += 1;
            return false;
        }
    }
    if blacklisted(&barcodes.cell) {
        stats.whitelist_filtered += 1;
        return false;
    }
    true
}

fn extract_and_filter<'a>(
    config: &ExtractConfig,
    mode: ExtractMode,
    r1: &'a SequenceRecord,
    r2: Option<&'a SequenceRecord>,
    stats: &mut ExtractStats,
) -> Result<Option<Extraction<'a>>, ExtractError> {
    let result = match mode {
        ExtractMode::Combine => extract_combined(config, r1, r2).map(Some),
        ExtractMode::EitherRead => {
            extract_either(config, r1, r2.expect("validated paired input"), stats)
        }
    };
    let mut extraction = match result {
        Ok(Some(extraction)) => extraction,
        Ok(None) => return Ok(None),
        Err(ExtractError::ReadTooShort { .. }) => {
            stats.too_short += 1;
            return Ok(None);
        }
        Err(ExtractError::RegexNoMatch) => {
            stats.no_match += 1;
            return Ok(None);
        }
        Err(error) => return Err(error),
    };
    Ok(filter_barcodes(config, &mut extraction.barcodes, stats).then_some(extraction))
}

fn write_output(
    writer: &mut Output<'_>,
    record: &SequenceRecord,
    header: &[u8],
    trimmed: Option<&TrimmedRead>,
) -> Result<(), ExtractError> {
    if let Some(writer) = writer {
        if let Some(trimmed) = trimmed {
            write_fastq_record(writer, header, &trimmed.sequence, &trimmed.quality)?;
        } else {
            let quality = record.qual().ok_or_else(|| {
                ExtractError::FastqParse("missing quality scores in FASTQ record".into())
            })?;
            write_fastq_record(writer, header, &record.seq(), quality)?;
        }
    }
    Ok(())
}

fn filtered_header(header: &[u8], strip_suffixes: bool) -> Vec<u8> {
    let name = read_name(header);
    let mut result = if strip_suffixes {
        strip_pair_suffix(name)
    } else {
        name
    }
    .to_vec();
    result.extend_from_slice(&header[name.len()..]);
    result
}

fn process_records(
    config: &ExtractConfig,
    mode: ExtractMode,
    r1: &SequenceRecord,
    r2: Option<&SequenceRecord>,
    stats: &mut ExtractStats,
    writers: &mut Writers<'_>,
) -> Result<(), ExtractError> {
    let Some(extraction) = extract_and_filter(config, mode, r1, r2, stats)? else {
        write_output(
            &mut writers.filtered1,
            r1,
            &filtered_header(r1.id(), config.ignore_read_pair_suffixes),
            None,
        )?;
        if let Some(r2) = r2 {
            write_output(
                &mut writers.filtered2,
                r2,
                &filtered_header(r2.id(), config.ignore_read_pair_suffixes),
                None,
            )?;
        }
        return Ok(());
    };
    let header = |record: &SequenceRecord| {
        build_read_name(
            record.id(),
            &extraction.barcodes.cell,
            &extraction.barcodes.umi,
            &config.umi_separator,
            config.ignore_read_pair_suffixes,
        )
    };
    let id1 = header(r1);
    write_output(&mut writers.read1, r1, &id1, extraction.read1.as_ref())?;
    if let Some(r2) = r2 {
        // Python's either-read mode copies the entire read1 header to read2.
        let id2 = if mode == ExtractMode::EitherRead {
            id1
        } else {
            header(r2)
        };
        write_output(&mut writers.read2, r2, &id2, extraction.read2.as_ref())?;
    }
    stats.output_reads += 1;
    Ok(())
}

fn validate(config: &ExtractConfig, mode: ExtractMode, paired: bool) -> Result<(), ExtractError> {
    let error = if config.pattern.is_none() && config.pattern2.is_none() {
        Some("at least one barcode pattern is required")
    } else if config.pattern2.is_some() && !paired {
        Some("a read2 pattern requires paired input")
    } else if mode == ExtractMode::EitherRead
        && (!paired || config.pattern.is_none() || config.pattern2.is_none())
    {
        Some("either-read extraction requires paired input and both patterns")
    } else {
        None
    };
    if let Some(error) = error {
        return Err(ExtractError::InvalidPattern(error.into()));
    }
    Ok(())
}

/// Extract, filter and write single or paired reads through one pipeline.
///
/// Pair matching and reconciliation happen before extraction. Output destinations
/// do not change how barcodes are selected or how filters are applied.
///
/// # Errors
/// Returns errors for invalid options, malformed FASTQ, unmatched pairs or I/O failures.
pub fn extract_with_outputs<R1: std::io::Read + Send, R2: std::io::Read + Send>(
    config: &ExtractConfig,
    mode: ExtractMode,
    input1: R1,
    input2: Option<R2>,
    outputs: ExtractOutputs<'_>,
) -> Result<ExtractStats, ExtractError> {
    validate(config, mode, input2.is_some())?;
    let mut reader1 = FastqReader::new(input1);
    let mut reader2 = input2.map(FastqReader::new);
    let mut writers = Writers::from(outputs);
    let mut stats = ExtractStats::default();
    while let Some(record) = reader1.next() {
        let r1 = record.map_err(|error| ExtractError::FastqParse(error.to_string()))?;
        stats.input_reads += 1;
        if config.past_subset(stats.input_reads) {
            break;
        }
        if let Some(reader2) = reader2.as_mut() {
            let name1 = pair_name(r1.id(), config.ignore_read_pair_suffixes);
            loop {
                let r2 = reader2
                    .next()
                    .ok_or_else(|| {
                        ExtractError::FastqParse(
                            "read2 exhausted before finding match for read1".into(),
                        )
                    })?
                    .map_err(|error| ExtractError::FastqParse(error.to_string()))?;
                let name2 = pair_name(r2.id(), config.ignore_read_pair_suffixes);
                if name2 != name1 {
                    if config.reconcile_pairs {
                        continue;
                    }
                    return Err(ExtractError::FastqParse(format!(
                        "read pairs do not match: {} != {}",
                        String::from_utf8_lossy(name1),
                        String::from_utf8_lossy(name2),
                    )));
                }
                process_records(config, mode, &r1, Some(&r2), &mut stats, &mut writers)?;
                break;
            }
        } else {
            process_records(config, mode, &r1, None, &mut stats, &mut writers)?;
        }
    }
    if !config.reconcile_pairs
        && !config.past_subset(stats.input_reads)
        && reader2
            .as_mut()
            .is_some_and(|reader| reader.next().is_some())
    {
        return Err(ExtractError::FastqParse(
            "read1 and read2 files have different numbers of records".into(),
        ));
    }
    writers.flush()?;
    Ok(stats)
}

// Compatibility entry points only adapt arguments; keep processing in extract_with_outputs.

/// Extract single-end reads using the shared pipeline.
///
/// # Errors
/// Returns errors for invalid options, malformed FASTQ or I/O failures.
pub fn extract_reads<R: std::io::Read + Send, W: Write>(
    config: &ExtractConfig,
    input: R,
    output: W,
) -> Result<ExtractStats, ExtractError> {
    extract_with_outputs(
        config,
        ExtractMode::Combine,
        input,
        None::<std::io::Empty>,
        ExtractOutputs {
            read1: Some(Box::new(output)),
            ..ExtractOutputs::default()
        },
    )
}

/// Extract paired reads, combining barcodes from the configured patterns.
///
/// # Errors
/// Returns errors for invalid options, malformed FASTQ, unmatched pairs or I/O failures.
pub fn extract_reads_paired<
    R1: std::io::Read + Send,
    R2: std::io::Read + Send,
    W1: Write,
    W2: Write,
>(
    config: &ExtractConfig,
    input1: R1,
    input2: R2,
    output1: W1,
    output2: W2,
) -> Result<ExtractStats, ExtractError> {
    extract_with_outputs(
        config,
        ExtractMode::Combine,
        input1,
        Some(input2),
        ExtractOutputs::paired(output1, output2),
    )
}

/// Extract paired reads and write read2 to the primary output.
///
/// # Errors
/// Returns errors for invalid options, malformed FASTQ, unmatched pairs or I/O failures.
pub fn extract_reads_paired_r1_pattern<
    R1: std::io::Read + Send,
    R2: std::io::Read + Send,
    W: Write,
>(
    config: &ExtractConfig,
    input1: R1,
    input2: R2,
    output: W,
    filtered_out1: Option<Box<dyn Write>>,
    filtered_out2: Option<Box<dyn Write>>,
) -> Result<ExtractStats, ExtractError> {
    extract_with_outputs(
        config,
        ExtractMode::Combine,
        input1,
        Some(input2),
        ExtractOutputs {
            read1: None,
            read2: Some(Box::new(output)),
            filtered1: filtered_out1,
            filtered2: filtered_out2,
        },
    )
}

/// Extract paired reads by choosing either matching pattern.
///
/// # Errors
/// Returns errors for invalid options, malformed FASTQ, unmatched pairs or I/O failures.
pub fn extract_reads_either_read<
    R1: std::io::Read + Send,
    R2: std::io::Read + Send,
    W1: Write,
    W2: Write,
>(
    config: &ExtractConfig,
    input1: R1,
    input2: R2,
    output1: W1,
    output2: W2,
) -> Result<ExtractStats, ExtractError> {
    extract_with_outputs(
        config,
        ExtractMode::EitherRead,
        input1,
        Some(input2),
        ExtractOutputs::paired(output1, output2),
    )
}
