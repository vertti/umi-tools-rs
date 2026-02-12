use std::io::{BufWriter, Write};

use needletail::parser::{FastqReader, FastxReader, SequenceRecord};

use crate::error::ExtractError;
use crate::pattern::BarcodePattern;

/// Configuration for the extract command.
#[derive(Debug, Clone)]
pub struct ExtractConfig {
    pub pattern: BarcodePattern,
    pub umi_separator: u8,
}

/// Statistics from an extraction run.
#[derive(Debug, Default)]
pub struct ExtractStats {
    pub input_reads: u64,
    pub output_reads: u64,
    pub too_short: u64,
    pub no_match: u64,
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
    let mut stats = ExtractStats::default();
    let mut writer = BufWriter::new(output);
    let mut reader = FastqReader::new(input);

    while let Some(result) = reader.next() {
        let record = result.map_err(|e| ExtractError::FastqParse(e.to_string()))?;
        stats.input_reads += 1;

        match process_record(&record, config) {
            Ok(processed) => {
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
}

fn process_record(
    record: &SequenceRecord,
    config: &ExtractConfig,
) -> Result<ProcessedRecord, ExtractError> {
    let seq = record.seq();
    let qual = record
        .qual()
        .ok_or_else(|| ExtractError::FastqParse("missing quality scores in FASTQ record".into()))?;

    let result = config.pattern.extract(&seq, qual)?;

    let id = build_read_name(
        record.id(),
        &result.cell_barcode,
        &result.umi,
        config.umi_separator,
    );

    Ok(ProcessedRecord {
        id,
        seq: result.trimmed_sequence,
        qual: result.trimmed_quality,
    })
}

/// Build a new read identifier with barcode(s) appended.
///
/// Format: `READ_ID{sep}UMI` or `READ_ID{sep}CELL{sep}UMI`
fn build_read_name(id: &[u8], cell: &[u8], umi: &[u8], separator: u8) -> Vec<u8> {
    let mut name = Vec::with_capacity(id.len() + 1 + cell.len() + 1 + umi.len());
    name.extend_from_slice(id);
    if !cell.is_empty() {
        name.push(separator);
        name.extend_from_slice(cell);
    }
    name.push(separator);
    name.extend_from_slice(umi);
    name
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
