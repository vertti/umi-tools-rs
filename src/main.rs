use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use umi_core::extract::{
    ExtractConfig, QualityEncoding, extract_reads, extract_reads_paired,
    extract_reads_paired_r1_pattern,
};
use umi_core::pattern::{BarcodePattern, PrimeEnd, RegexPattern, StringPattern};

#[derive(Parser)]
#[command(name = "umi-tools-rs", version, about = "Fast UMI tools in Rust")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract UMI from FASTQ reads
    Extract {
        /// Barcode pattern for read1 (e.g. NNNXXXXNN). N=UMI, C=cell, X=discard.
        #[arg(long = "bc-pattern")]
        bc_pattern: Option<String>,

        /// Barcode pattern for read2 (paired-end mode)
        #[arg(long = "bc-pattern2")]
        bc_pattern2: Option<String>,

        /// Extraction method: "string" for fixed-position, "regex" for named capture groups
        #[arg(long = "extract-method", default_value = "string")]
        extract_method: String,

        /// Input FASTQ file (default: stdin). Gzip detected by .gz extension.
        #[arg(short = 'I', long = "stdin")]
        input: Option<String>,

        /// Output FASTQ file (default: stdout). Gzip if .gz extension.
        #[arg(short = 'S', long = "stdout")]
        output: Option<String>,

        /// Read2 input FASTQ file (paired-end mode)
        #[arg(long = "read2-in")]
        read2_in: Option<String>,

        /// Read2 output FASTQ file (paired-end mode)
        #[arg(long = "read2-out")]
        read2_out: Option<String>,

        /// Write read2 to stdout (paired-end mode, pattern on read1)
        #[arg(long = "read2-stdout")]
        read2_stdout: bool,

        /// Whitelist file of accepted cell barcodes (one per line)
        #[arg(long = "whitelist")]
        whitelist: Option<String>,

        /// Extract from 3' end instead of 5'
        #[arg(long = "3prime")]
        prime3: bool,

        /// UMI separator character in read name
        #[arg(long = "umi-separator", default_value = "_")]
        umi_separator: String,

        /// Minimum per-base quality score for UMI bases (reads below are discarded)
        #[arg(long = "quality-filter-threshold")]
        quality_filter_threshold: Option<u8>,

        /// Quality encoding scheme: phred33, phred64, solexa
        #[arg(long = "quality-encoding", default_value = "phred33")]
        quality_encoding: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Extract {
            bc_pattern,
            bc_pattern2,
            extract_method,
            input,
            output,
            read2_in,
            read2_out,
            read2_stdout,
            whitelist,
            prime3,
            umi_separator,
            quality_filter_threshold,
            quality_encoding,
        } => {
            let is_paired = read2_in.is_some();
            if !is_paired && bc_pattern.is_none() {
                bail!("--bc-pattern is required for single-end extraction");
            }
            if is_paired && bc_pattern.is_none() && bc_pattern2.is_none() {
                bail!("at least one of --bc-pattern or --bc-pattern2 is required");
            }

            run_extract(
                bc_pattern.as_deref(),
                bc_pattern2.as_deref(),
                &extract_method,
                input.as_deref(),
                output.as_deref(),
                read2_in.as_deref(),
                read2_out.as_deref(),
                read2_stdout,
                whitelist.as_deref(),
                prime3,
                &umi_separator,
                quality_filter_threshold,
                &quality_encoding,
            )
        }
    }
}

fn parse_pattern(raw: &str, extract_method: &str, prime3: bool) -> Result<BarcodePattern> {
    match extract_method {
        "string" => {
            let prime_end = if prime3 {
                PrimeEnd::Three
            } else {
                PrimeEnd::Five
            };
            Ok(BarcodePattern::String(
                StringPattern::parse(raw, prime_end).context("failed to parse barcode pattern")?,
            ))
        }
        "regex" => Ok(BarcodePattern::Regex(
            RegexPattern::parse(raw).context("failed to parse regex pattern")?,
        )),
        other => bail!("unknown extract method '{other}'; expected 'string' or 'regex'"),
    }
}

fn open_input(path: Option<&str>) -> Result<Box<dyn Read + Send>> {
    match path {
        Some(p) => {
            let file = File::open(p).with_context(|| format!("failed to open input file: {p}"))?;
            if is_gzipped(p) {
                Ok(Box::new(MultiGzDecoder::new(file)))
            } else {
                Ok(Box::new(file))
            }
        }
        None => Ok(Box::new(io::stdin())),
    }
}

fn open_output(path: Option<&str>) -> Result<Box<dyn Write>> {
    match path {
        Some(p) => {
            let file =
                File::create(p).with_context(|| format!("failed to create output file: {p}"))?;
            if is_gzipped(p) {
                Ok(Box::new(GzEncoder::new(file, Compression::default())))
            } else {
                Ok(Box::new(file))
            }
        }
        None => Ok(Box::new(io::stdout().lock())),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_extract(
    bc_pattern: Option<&str>,
    bc_pattern2: Option<&str>,
    extract_method: &str,
    input_path: Option<&str>,
    output_path: Option<&str>,
    read2_in_path: Option<&str>,
    read2_out_path: Option<&str>,
    read2_stdout: bool,
    whitelist_path: Option<&str>,
    prime3: bool,
    umi_separator: &str,
    quality_filter_threshold: Option<u8>,
    quality_encoding: &str,
) -> Result<()> {
    let pattern = bc_pattern
        .map(|p| parse_pattern(p, extract_method, prime3))
        .transpose()?;
    let pattern2 = bc_pattern2
        .map(|p| parse_pattern(p, extract_method, prime3))
        .transpose()?;

    let sep_byte = umi_separator.as_bytes().first().copied().unwrap_or(b'_');

    let qe = match quality_encoding {
        "phred33" => QualityEncoding::Phred33,
        "phred64" => QualityEncoding::Phred64,
        "solexa" => QualityEncoding::Solexa,
        other => {
            bail!("unknown quality encoding '{other}'; expected 'phred33', 'phred64', or 'solexa'")
        }
    };

    let whitelist = whitelist_path.map(load_whitelist).transpose()?;

    let config = ExtractConfig {
        pattern,
        pattern2,
        umi_separator: sep_byte,
        quality_filter_threshold,
        quality_encoding: qe,
        whitelist,
    };

    let reader1 = open_input(input_path)?;

    let stats = if let Some(r2_path) = read2_in_path {
        let reader2 = open_input(Some(r2_path))?;
        if read2_stdout {
            let writer = open_output(output_path)?;
            extract_reads_paired_r1_pattern(&config, reader1, reader2, writer)
                .context("paired-end extraction failed")?
        } else {
            let writer1 = open_output(output_path)?;
            let writer2 = open_output(read2_out_path)
                .context("--read2-out is required when --read2-in is specified")?;
            extract_reads_paired(&config, reader1, reader2, writer1, writer2)
                .context("paired-end extraction failed")?
        }
    } else {
        let writer1 = open_output(output_path)?;
        extract_reads(&config, reader1, writer1).context("extraction failed")?
    };

    eprintln!(
        "Reads input: {}, output: {}, too short: {}, no match: {}, quality filtered: {}, whitelist filtered: {}",
        stats.input_reads,
        stats.output_reads,
        stats.too_short,
        stats.no_match,
        stats.quality_filtered,
        stats.whitelist_filtered,
    );

    Ok(())
}

fn load_whitelist(path: &str) -> Result<HashSet<Vec<u8>>> {
    let file =
        File::open(path).with_context(|| format!("failed to open whitelist file: {path}"))?;
    let reader = io::BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("failed to read whitelist file: {path}"))?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            set.insert(trimmed.as_bytes().to_vec());
        }
    }
    Ok(set)
}

fn is_gzipped(path: &str) -> bool {
    Path::new(path)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}
