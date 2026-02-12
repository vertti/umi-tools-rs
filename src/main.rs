use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use umi_core::extract::{ExtractConfig, QualityEncoding, extract_reads};
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
        /// Barcode pattern (e.g. NNNXXXXNN). N=UMI, C=cell, X=discard.
        #[arg(long = "bc-pattern")]
        bc_pattern: String,

        /// Extraction method: "string" for fixed-position, "regex" for named capture groups
        #[arg(long = "extract-method", default_value = "string")]
        extract_method: String,

        /// Input FASTQ file (default: stdin). Gzip detected by .gz extension.
        #[arg(short = 'I', long = "stdin")]
        input: Option<String>,

        /// Output FASTQ file (default: stdout). Gzip if .gz extension.
        #[arg(short = 'S', long = "stdout")]
        output: Option<String>,

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
            extract_method,
            input,
            output,
            prime3,
            umi_separator,
            quality_filter_threshold,
            quality_encoding,
        } => run_extract(
            &bc_pattern,
            &extract_method,
            input.as_deref(),
            output.as_deref(),
            prime3,
            &umi_separator,
            quality_filter_threshold,
            &quality_encoding,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_extract(
    bc_pattern: &str,
    extract_method: &str,
    input_path: Option<&str>,
    output_path: Option<&str>,
    prime3: bool,
    umi_separator: &str,
    quality_filter_threshold: Option<u8>,
    quality_encoding: &str,
) -> Result<()> {
    let pattern = match extract_method {
        "string" => {
            let prime_end = if prime3 {
                PrimeEnd::Three
            } else {
                PrimeEnd::Five
            };
            BarcodePattern::String(
                StringPattern::parse(bc_pattern, prime_end)
                    .context("failed to parse barcode pattern")?,
            )
        }
        "regex" => BarcodePattern::Regex(
            RegexPattern::parse(bc_pattern).context("failed to parse regex pattern")?,
        ),
        other => bail!("unknown extract method '{other}'; expected 'string' or 'regex'"),
    };

    let sep_byte = umi_separator.as_bytes().first().copied().unwrap_or(b'_');

    let qe = match quality_encoding {
        "phred33" => QualityEncoding::Phred33,
        "phred64" => QualityEncoding::Phred64,
        "solexa" => QualityEncoding::Solexa,
        other => {
            bail!("unknown quality encoding '{other}'; expected 'phred33', 'phred64', or 'solexa'")
        }
    };

    let config = ExtractConfig {
        pattern,
        umi_separator: sep_byte,
        quality_filter_threshold,
        quality_encoding: qe,
    };

    let reader: Box<dyn Read + Send> = match input_path {
        Some(path) => {
            let file =
                File::open(path).with_context(|| format!("failed to open input file: {path}"))?;
            if is_gzipped(path) {
                Box::new(MultiGzDecoder::new(file))
            } else {
                Box::new(file)
            }
        }
        None => Box::new(io::stdin()),
    };

    let writer: Box<dyn Write> = match output_path {
        Some(path) => {
            let file = File::create(path)
                .with_context(|| format!("failed to create output file: {path}"))?;
            if is_gzipped(path) {
                Box::new(GzEncoder::new(file, Compression::default()))
            } else {
                Box::new(file)
            }
        }
        None => Box::new(io::stdout().lock()),
    };

    let stats = extract_reads(&config, reader, writer).context("extraction failed")?;

    eprintln!(
        "Reads input: {}, output: {}, too short: {}, no match: {}, quality filtered: {}",
        stats.input_reads,
        stats.output_reads,
        stats.too_short,
        stats.no_match,
        stats.quality_filtered
    );

    Ok(())
}

fn is_gzipped(path: &str) -> bool {
    Path::new(path)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}
