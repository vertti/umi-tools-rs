use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use umi_core::extract::{ExtractConfig, extract_reads};
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
        } => run_extract(
            &bc_pattern,
            &extract_method,
            input.as_deref(),
            output.as_deref(),
            prime3,
            &umi_separator,
        ),
    }
}

fn run_extract(
    bc_pattern: &str,
    extract_method: &str,
    input_path: Option<&str>,
    output_path: Option<&str>,
    prime3: bool,
    umi_separator: &str,
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

    let config = ExtractConfig {
        pattern,
        umi_separator: sep_byte,
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
        "Reads input: {}, output: {}, too short: {}, no match: {}",
        stats.input_reads, stats.output_reads, stats.too_short, stats.no_match
    );

    Ok(())
}

fn is_gzipped(path: &str) -> bool {
    Path::new(path)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}
