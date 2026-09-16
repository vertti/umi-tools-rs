use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Read, Write};
use std::path::Path;
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, PoisonError};

use anyhow::{Context, Result, bail};
use clap::{ArgAction, Parser, Subcommand};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use umi_core::alignment_io::{AlignmentFormat, determine_format};
use umi_core::count::{CountConfig, CountTabConfig, run_count, run_count_tab};
use umi_core::dedup::{
    DedupConfig, DedupMethod, MultimappingDetection, PositionOptions, run_dedup,
};
use umi_core::extract::{
    ExtractConfig, QualityEncoding, extract_reads, extract_reads_either_read, extract_reads_paired,
    extract_reads_paired_r1_pattern,
};
use umi_core::group::{ChimericPairs, GroupConfig, UnmappedHandling, run_group};
use umi_core::pattern::{BarcodePattern, PrimeEnd, RegexPattern, StringPattern};
use umi_core::whitelist::{EdAboveThreshold, KneeMethod, WhitelistConfig, run_whitelist};

#[derive(Parser)]
#[command(
    name = "umi-tools-rs",
    version,
    about = "Fast UMI tools in Rust",
    infer_long_args = true,
    propagate_version = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract UMI from FASTQ reads
    Extract(ExtractArgs),

    /// Build a whitelist of valid cell barcodes from FASTQ
    Whitelist(WhitelistArgs),

    /// Group PCR duplicates in BAM by UMI and mapping position
    Group(GroupArgs),

    /// Deduplicate BAM reads based on UMI and mapping position
    Dedup(DedupArgs),

    /// Count UMI-deduplicated reads per gene from BAM
    Count(CountArgs),

    /// Count UMI-deduplicated reads per gene from tab-delimited input
    #[command(name = "count_tab")]
    CountTab(CountTabArgs),
}

#[derive(clap::Args)]
#[allow(clippy::struct_excessive_bools)]
struct ExtractArgs {
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

    /// Strip /1 and /2 suffixes from read names before appending UMI
    #[arg(long = "ignore-read-pair-suffixes")]
    ignore_read_pair_suffixes: bool,

    /// Reconcile read pairs when read1 is a pre-filtered subset of read2
    #[arg(long = "reconcile-pairs")]
    reconcile_pairs: bool,

    /// Error-correct cell barcodes using whitelist correction map
    #[arg(long = "error-correct-cell")]
    error_correct_cell: bool,

    /// Blacklist file of rejected cell barcodes (one per line)
    #[arg(long = "blacklist")]
    blacklist: Option<String>,

    /// Output file for filtered read1 (reads that fail any filter)
    #[arg(long = "filtered-out")]
    filtered_out: Option<String>,

    /// Output file for filtered read2 (reads that fail any filter)
    #[arg(long = "filtered-out2")]
    filtered_out2: Option<String>,

    /// Either-read mode: try pattern on both reads, use whichever matches
    #[arg(long = "either-read")]
    either_read: bool,

    /// No randomness in this command; accepted for umi-tools compatibility
    #[arg(long = "random-seed")]
    _random_seed: Option<u64>,

    #[command(flatten)]
    common: CommonArgs,
}

#[derive(clap::Args)]
struct WhitelistArgs {
    /// Barcode pattern (e.g. CCCCCCNNNNNNNNNN). N=UMI, C=cell, X=discard.
    #[arg(long = "bc-pattern")]
    bc_pattern: String,

    /// Extraction method: "string" for fixed-position, "regex" for named capture groups
    #[arg(long = "extract-method", default_value = "string")]
    extract_method: String,

    /// Input FASTQ file (default: stdin). Gzip detected by .gz extension.
    #[arg(short = 'I', long = "stdin")]
    input: Option<String>,

    /// Output TSV file (default: stdout).
    #[arg(short = 'S', long = "stdout")]
    output: Option<String>,

    /// Extract from 3' end instead of 5'
    #[arg(long = "3prime")]
    prime3: bool,

    /// Knee detection method: "distance" or "density"
    #[arg(long = "knee-method", default_value = "distance")]
    knee_method: String,

    /// Force whitelist to include this many cells
    #[arg(long = "set-cell-number")]
    set_cell_number: Option<usize>,

    /// Expected number of cells (hint for density method)
    #[arg(long = "expect-cells")]
    expect_cells: Option<usize>,

    /// Maximum Hamming distance for error correction (default: 1)
    #[arg(long = "error-correct-threshold", default_value = "1")]
    error_correct_threshold: usize,

    /// Handle whitelist barcodes within edit distance of higher-count whitelist barcode
    #[arg(long = "ed-above-threshold")]
    ed_above_threshold: Option<String>,

    /// Prefix for knee plots. Accepted for umi-tools compatibility; plots are not generated
    #[arg(long = "plot-prefix")]
    plot_prefix: Option<String>,

    /// Output file for reads that failed barcode extraction
    #[arg(long = "filtered-out")]
    filtered_out: Option<String>,

    /// Max reads to process (default: `100_000_000`)
    #[arg(long = "subset-reads", default_value = "100000000")]
    subset_reads: usize,

    /// No randomness in this command; accepted for umi-tools compatibility
    #[arg(long = "random-seed")]
    _random_seed: Option<u64>,

    #[command(flatten)]
    common: CommonArgs,
}

#[derive(clap::Args)]
#[allow(clippy::struct_excessive_bools)]
struct GroupArgs {
    /// Input BAM file
    #[arg(short = 'I', long = "stdin")]
    input: Option<String>,

    /// Grouping method: unique, percentile, cluster, adjacency, directional
    #[arg(long = "method", default_value = "directional")]
    method: String,

    /// Ignore UMI — group by position only
    #[arg(long = "ignore-umi")]
    ignore_umi: bool,

    /// Output file for the tagged alignments (default: stdout; requires --output-bam)
    #[arg(short = 'S', long = "stdout")]
    output: Option<String>,

    #[command(flatten)]
    input_format: InputFormatArgs,

    #[command(flatten)]
    output_format: OutputFormatArgs,

    /// Random seed for reproducible tie-breaking
    #[arg(long = "random-seed", default_value = "0")]
    random_seed: u64,

    /// UMI separator in read name
    #[arg(long = "umi-separator", default_value = "_")]
    umi_separator: String,

    /// Only process reads on this chromosome
    #[arg(long = "chrom")]
    chrom: Option<String>,

    /// Output TSV file with group assignments
    #[arg(long = "group-out")]
    group_out: Option<String>,

    /// Write tagged BAM/SAM to stdout
    #[arg(long = "output-bam")]
    output_bam: bool,

    /// Skip coordinate sorting of output
    #[arg(long = "no-sort-output")]
    no_sort_output: bool,

    /// Random subset of reads to process (0.0-1.0)
    #[arg(long = "subset")]
    subset: Option<f32>,

    #[command(flatten)]
    position: PositionArgs,

    /// Minimum mapping quality for a read to be retained
    #[arg(long = "mapping-quality", default_value = "0")]
    mapping_quality: u8,

    /// Accepted for umi-tools compatibility; group keeps every read, so the tag is never consulted
    #[arg(long = "multimapping-detection-method", value_parser = ["NH", "X0", "XT"])]
    multimapping_detection_method: Option<String>,

    /// Buffer a whole contig before grouping instead of a 1 kb window; uses more memory
    #[arg(long = "buffer-whole-contig", alias = "whole-contig", action = ArgAction::SetTrue, overrides_with = "buffer_whole_contig")]
    buffer_whole_contig: bool,

    /// Include unmapped reads in output (alias for --unmapped-reads=output)
    #[arg(long = "output-unmapped")]
    output_unmapped: bool,

    /// Enable paired-end grouping
    #[arg(long = "paired")]
    paired: bool,

    /// How to handle chimeric read pairs: discard, output, use
    #[arg(long = "chimeric-pairs")]
    chimeric_pairs: Option<String>,

    /// How to handle unmapped reads: discard, output, use
    #[arg(long = "unmapped-reads", default_value = "discard")]
    unmapped: String,

    /// Deduplicate per gene (requires --gene-tag or --per-contig)
    #[arg(long = "per-gene")]
    per_gene: bool,

    /// BAM tag containing gene assignment (default: XF)
    #[arg(long = "gene-tag")]
    gene_tag: Option<String>,

    /// Skip reads with gene tag matching this regex
    #[arg(long = "skip-tags-regex")]
    skip_tags_regex: Option<String>,

    /// Use contig name as gene (requires --per-gene)
    #[arg(long = "per-contig")]
    per_contig: bool,

    #[command(flatten)]
    common: CommonArgs,
}

#[derive(clap::Args)]
#[allow(clippy::struct_excessive_bools)]
struct DedupArgs {
    /// Input BAM file
    #[arg(short = 'I', long = "stdin")]
    input: Option<String>,

    /// Dedup method: unique, percentile, cluster, adjacency, directional
    #[arg(long = "method", default_value = "directional")]
    method: String,

    /// Ignore UMI — deduplicate by position only
    #[arg(long = "ignore-umi")]
    ignore_umi: bool,

    /// Output file (default: stdout)
    #[arg(short = 'S', long = "stdout")]
    output: Option<String>,

    #[command(flatten)]
    input_format: InputFormatArgs,

    #[command(flatten)]
    output_format: OutputFormatArgs,

    /// Random seed for reproducible tie-breaking
    #[arg(long = "random-seed", default_value = "0")]
    random_seed: u64,

    /// UMI separator in read name
    #[arg(long = "umi-separator", default_value = "_")]
    umi_separator: String,

    /// Only process reads on this chromosome
    #[arg(long = "chrom")]
    chrom: Option<String>,

    /// Edit distance threshold for UMI clustering
    #[arg(long = "edit-distance-threshold", default_value = "1")]
    edit_distance_threshold: u32,

    #[command(flatten)]
    position: PositionArgs,

    /// Minimum mapping quality for a read to be retained
    #[arg(long = "mapping-quality", default_value = "0")]
    mapping_quality: u8,

    /// Tag that records multimapping (NH, X0 or XT); among duplicates with equal MAPQ the read with fewer hits is kept
    #[arg(long = "multimapping-detection-method", value_parser = ["NH", "X0", "XT"])]
    multimapping_detection_method: Option<String>,

    /// Buffer a whole contig before grouping instead of a 1 kb window; uses more memory
    #[arg(long = "buffer-whole-contig", alias = "whole-contig", action = ArgAction::SetTrue, overrides_with = "buffer_whole_contig")]
    buffer_whole_contig: bool,

    /// Random subset of reads to process (0.0-1.0)
    #[arg(long = "subset")]
    subset: Option<f32>,

    /// UMI extraction method: `read_id` or `tag`
    #[arg(long = "extract-umi-method", default_value = "read_id")]
    extract_umi_method: String,

    /// BAM tag to extract UMI from (when extract-umi-method=tag)
    #[arg(long = "umi-tag")]
    umi_tag: Option<String>,

    /// Deduplicate per gene (requires --gene-tag)
    #[arg(long = "per-gene")]
    per_gene: bool,

    /// BAM tag containing gene assignment
    #[arg(long = "gene-tag")]
    gene_tag: Option<String>,

    /// Skip reads with gene tag matching this regex
    #[arg(long = "skip-tags-regex")]
    skip_tags_regex: Option<String>,

    /// Output stats file prefix
    #[arg(long = "output-stats")]
    output_stats: Option<String>,

    /// Enable paired-end deduplication
    #[arg(long = "paired")]
    paired: bool,

    /// Ignore template length when grouping reads (paired mode)
    #[arg(long = "ignore-tlen")]
    ignore_tlen: bool,

    /// Filter UMIs against whitelist
    #[arg(long = "filter-umi")]
    filter_umi: bool,

    /// UMI whitelist file (or read1 whitelist for paired UMIs)
    #[arg(long = "umi-whitelist")]
    umi_whitelist: Option<String>,

    /// Read2 UMI whitelist file (paired UMI mode, Cartesian product with read1)
    #[arg(long = "umi-whitelist-paired")]
    umi_whitelist_paired: Option<String>,

    #[command(flatten)]
    common: CommonArgs,
}

#[derive(clap::Args)]
struct CountArgs {
    /// Input BAM file
    #[arg(short = 'I', long = "stdin")]
    input: Option<String>,

    /// Output file (default: stdout)
    #[arg(short = 'S', long = "stdout")]
    output: Option<String>,

    #[command(flatten)]
    input_format: InputFormatArgs,

    /// Minimum mapping quality for a read to be retained
    #[arg(long = "mapping-quality", default_value = "0")]
    mapping_quality: u8,

    /// Dedup method: unique, percentile, cluster, adjacency, directional
    #[arg(long = "method", default_value = "directional")]
    method: String,

    /// BAM tag containing gene assignment
    #[arg(long = "gene-tag", default_value = "XF")]
    gene_tag: String,

    /// Skip reads with gene tag matching this regex
    #[arg(long = "skip-tags-regex")]
    skip_tags_regex: Option<String>,

    /// UMI extraction method: umis, `read_id`, tag
    #[arg(long = "extract-umi-method", default_value = "read_id")]
    extract_umi_method: String,

    /// Count per cell barcode
    #[arg(long = "per-cell")]
    per_cell: bool,

    /// Output wide-format cell counts (requires --per-cell)
    #[arg(long = "wide-format-cell-counts")]
    wide_format: bool,

    /// Edit distance threshold for UMI clustering
    #[arg(long = "edit-distance-threshold", default_value = "1")]
    edit_distance_threshold: u32,

    /// No randomness in this command; accepted for umi-tools compatibility
    #[arg(long = "random-seed")]
    _random_seed: Option<u64>,

    #[command(flatten)]
    common: CommonArgs,
}

#[derive(clap::Args)]
struct CountTabArgs {
    /// Input TSV file (default: stdin)
    #[arg(short = 'I', long = "stdin")]
    input: Option<String>,

    /// Output file (default: stdout)
    #[arg(short = 'S', long = "stdout")]
    output: Option<String>,

    /// Count per cell barcode
    #[arg(long = "per-cell")]
    per_cell: bool,

    /// Barcode separator in read name
    #[arg(long = "barcode-separator", default_value = "_")]
    separator: String,

    /// Dedup method: unique, percentile, cluster, adjacency, directional
    #[arg(long = "method", default_value = "directional")]
    method: String,

    /// Edit distance threshold for UMI clustering
    #[arg(long = "edit-distance-threshold", default_value = "1")]
    edit_distance_threshold: u32,

    /// No randomness in this command; accepted for umi-tools compatibility
    #[arg(long = "random-seed")]
    _random_seed: Option<u64>,

    #[command(flatten)]
    common: CommonArgs,
}

/// Options shared by the commands that read alignments.
#[derive(clap::Args)]
struct InputFormatArgs {
    /// Input format: sam, bam or cram. Detected from the file content, so this has no effect.
    #[arg(long = "in-format", value_parser = ["sam", "bam", "cram"])]
    _in_format: Option<String>,

    /// Input is SAM. Detected from the file content, so this has no effect.
    #[arg(short = 'i', long = "in-sam", action = ArgAction::SetTrue, overrides_with = "_in_sam")]
    _in_sam: bool,

    /// FASTA reference for reading and writing CRAM. Local path only; defaults to the UR field of the input header.
    #[arg(long = "reference-filename")]
    reference_filename: Option<String>,

    /// htslib format options for reading. Accepted for umi-tools compatibility; has no effect.
    #[arg(long = "input-options")]
    input_options: Option<String>,
}

impl InputFormatArgs {
    fn note_ignored_flags(&self) {
        if self.input_options.is_some() {
            note_ignored("--input-options");
        }
    }
}

/// Grouping-key options shared by dedup and group.
#[derive(clap::Args, Clone, Copy)]
struct PositionArgs {
    /// Treat a spliced read as different from an unspliced one at the same position
    #[arg(long = "spliced-is-unique", action = ArgAction::SetTrue, overrides_with = "spliced_is_unique")]
    spliced_is_unique: bool,

    /// Bases soft-clipped from the 5' end before a read counts as spliced
    #[arg(long = "soft-clip-threshold", default_value = "4")]
    soft_clip_threshold: f64,

    /// Use read length as well as position and UMI to identify duplicates
    #[arg(long = "read-length", action = ArgAction::SetTrue, overrides_with = "read_length")]
    read_length: bool,
}

impl PositionArgs {
    const fn options(self) -> PositionOptions {
        PositionOptions {
            spliced_is_unique: self.spliced_is_unique,
            soft_clip_threshold: self.soft_clip_threshold,
            read_length: self.read_length,
        }
    }
}

/// Options shared by the commands that write alignments.
#[derive(clap::Args)]
struct OutputFormatArgs {
    /// Output format: sam, bam or cram (default: from the --stdout extension, else bam)
    #[arg(long = "out-format", value_parser = ["sam", "bam", "cram"])]
    out_format: Option<String>,

    /// Output SAM (same as --out-format=sam)
    #[arg(short = 'o', long = "out-sam", action = ArgAction::SetTrue, overrides_with = "out_sam")]
    out_sam: bool,

    /// htslib format options for writing. Accepted for umi-tools compatibility; has no effect.
    #[arg(long = "output-options")]
    output_options: Option<String>,
}

impl OutputFormatArgs {
    fn note_ignored_flags(&self) {
        if self.output_options.is_some() {
            note_ignored("--output-options");
        }
    }

    fn resolve(&self, output_path: Option<&str>) -> Result<AlignmentFormat> {
        let explicit = self
            .out_format
            .as_deref()
            .map(AlignmentFormat::parse)
            .transpose()
            .map_err(|name| anyhow::anyhow!("unknown output format '{name}'"))?;
        Ok(determine_format(output_path, self.out_sam, explicit))
    }
}

/// Logging options shared with `umi_tools`.
#[derive(clap::Args, Clone)]
struct CommonArgs {
    /// Append the run summary to this file instead of printing it to stderr
    #[arg(short = 'L', long = "log")]
    log: Option<String>,

    /// Print the run summary to stderr (the default)
    #[arg(long = "log2stderr", action = ArgAction::SetTrue, overrides_with = "_log2stderr")]
    _log2stderr: bool,

    /// Write umi-tools-rs notes and errors to this file instead of stderr. htslib messages still go to stderr.
    #[arg(short = 'E', long = "error")]
    error: Option<String>,

    /// Verbosity: 0 silences the run summary and notes; higher levels have no further effect
    #[arg(short = 'v', long = "verbose", default_value = "1")]
    verbose: u8,

    /// gzip level for .gz outputs, 1-9. umi-tools defaults to 6.
    #[arg(long = "compresslevel", default_value = "3", value_parser = clap::value_parser!(u32).range(1..=9))]
    compresslevel: u32,

    /// Same as --help
    #[arg(long = "help-extended", action = ArgAction::Help)]
    _help_extended: Option<bool>,

    /// Directory for temporary files. Accepted for umi-tools compatibility; has no effect.
    #[arg(long = "temp-dir")]
    temp_dir: Option<String>,

    /// Timing output file. Accepted for umi-tools compatibility; has no effect.
    #[arg(long = "timeit")]
    timeit: Option<String>,

    /// Name for the timing row. Accepted for umi-tools compatibility; has no effect.
    #[arg(long = "timeit-name")]
    timeit_name: Option<String>,

    /// Write a header to the timing file. Accepted for umi-tools compatibility; has no effect.
    #[arg(long = "timeit-header", action = ArgAction::SetTrue, overrides_with = "timeit_header")]
    timeit_header: bool,
}

impl CommonArgs {
    fn note_ignored_flags(&self) {
        let flags = [
            ("--temp-dir", self.temp_dir.is_some()),
            ("--timeit", self.timeit.is_some()),
            ("--timeit-name", self.timeit_name.is_some()),
            ("--timeit-header", self.timeit_header),
        ];
        for (flag, given) in flags {
            if given {
                note_ignored(flag);
            }
        }
    }
}

impl Commands {
    const fn common(&self) -> &CommonArgs {
        match self {
            Self::Extract(args) => &args.common,
            Self::Whitelist(args) => &args.common,
            Self::Group(args) => &args.common,
            Self::Dedup(args) => &args.common,
            Self::Count(args) => &args.common,
            Self::CountTab(args) => &args.common,
        }
    }

    fn note_ignored_flags(&self) {
        self.common().note_ignored_flags();
        if let Self::Whitelist(args) = self
            && args.plot_prefix.is_some()
        {
            note("--plot-prefix is accepted for umi-tools compatibility; plots are not generated");
        }
        if let Self::Group(args) = self
            && args.multimapping_detection_method.is_some()
        {
            note(
                "--multimapping-detection-method is accepted for umi-tools compatibility; \
                 group keeps every read, so it has no effect",
            );
        }
    }
}

/// Where the run summary goes: stderr, or the `--log` file appended to as `umi_tools` does.
struct RunLog {
    file: Option<File>,
    quiet: bool,
}

impl RunLog {
    fn open(path: Option<&str>, args: &[String], quiet: bool) -> Result<Self> {
        let file = path
            .map(|p| {
                let mut file = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(p)
                    .with_context(|| format!("failed to open log file: {p}"))?;
                if !quiet {
                    writeln!(
                        file,
                        "# umi-tools-rs version: {}",
                        env!("CARGO_PKG_VERSION")
                    )?;
                    writeln!(file, "# output generated by {}", args.join(" "))?;
                    writeln!(file, "# job started at {}", timestamp())?;
                }
                Ok::<_, anyhow::Error>(file)
            })
            .transpose()?;
        Ok(Self { file, quiet })
    }

    fn finish(mut self, summary: &str) -> Result<()> {
        if self.quiet {
            return Ok(());
        }
        match self.file.as_mut() {
            Some(file) => {
                writeln!(file, "{summary}")?;
                writeln!(file, "# job finished at {}", timestamp())?;
            }
            None => eprintln!("{summary}"),
        }
        Ok(())
    }
}

fn timestamp() -> humantime::Rfc3339Timestamp {
    humantime::format_rfc3339_seconds(std::time::SystemTime::now())
}

/// Destination for notes and the final error message: stderr, or the `--error` file.
static DIAGNOSTICS: Mutex<Option<File>> = Mutex::new(None);
static QUIET: AtomicBool = AtomicBool::new(false);

fn install_diagnostics(error_path: Option<&str>, quiet: bool) -> Result<()> {
    QUIET.store(quiet, Ordering::Relaxed);
    if let Some(path) = error_path {
        let file =
            File::create(path).with_context(|| format!("failed to create error file: {path}"))?;
        *DIAGNOSTICS.lock().unwrap_or_else(PoisonError::into_inner) = Some(file);
    }
    Ok(())
}

fn diagnostic(message: &str) {
    let mut sink = DIAGNOSTICS.lock().unwrap_or_else(PoisonError::into_inner);
    match sink.as_mut() {
        Some(file) => {
            let _ = writeln!(file, "{message}");
        }
        None => eprintln!("{message}"),
    }
}

fn note(message: &str) {
    if !QUIET.load(Ordering::Relaxed) {
        diagnostic(&format!("note: {message}"));
    }
}

fn note_ignored(flag: &str) {
    note(&format!(
        "{flag} is accepted for umi-tools compatibility and has no effect"
    ));
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let common = cli.command.common().clone();
    let args: Vec<String> = std::env::args().skip(1).collect();

    match run_logged(cli.command, &common, &args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            diagnostic(&format!("Error: {e:?}"));
            ExitCode::FAILURE
        }
    }
}

fn run_logged(command: Commands, common: &CommonArgs, args: &[String]) -> Result<()> {
    let quiet = common.verbose == 0;
    install_diagnostics(common.error.as_deref(), quiet)?;
    command.note_ignored_flags();
    let log = RunLog::open(common.log.as_deref(), args, quiet)?;
    let summary = run(command)?;
    log.finish(&summary)
}

#[allow(clippy::too_many_lines)]
fn run(command: Commands) -> Result<String> {
    match command {
        Commands::Extract(ExtractArgs {
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
            ignore_read_pair_suffixes,
            reconcile_pairs,
            error_correct_cell,
            blacklist,
            filtered_out,
            filtered_out2,
            either_read,
            _random_seed: _,
            common,
        }) => {
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
                ignore_read_pair_suffixes,
                reconcile_pairs,
                error_correct_cell,
                blacklist.as_deref(),
                filtered_out.as_deref(),
                filtered_out2.as_deref(),
                either_read,
                common.compresslevel,
            )
        }
        Commands::Whitelist(WhitelistArgs {
            bc_pattern,
            extract_method,
            input,
            output,
            prime3,
            knee_method,
            set_cell_number,
            expect_cells,
            error_correct_threshold,
            ed_above_threshold,
            plot_prefix: _,
            filtered_out,
            subset_reads,
            _random_seed: _,
            common,
        }) => run_whitelist_cmd(
            &bc_pattern,
            &extract_method,
            input.as_deref(),
            output.as_deref(),
            prime3,
            &knee_method,
            set_cell_number,
            expect_cells,
            error_correct_threshold,
            ed_above_threshold.as_deref(),
            filtered_out.as_deref(),
            subset_reads,
            common.compresslevel,
        ),
        Commands::Group(GroupArgs {
            input,
            method,
            ignore_umi,
            output,
            input_format,
            output_format,
            random_seed,
            umi_separator,
            chrom,
            group_out,
            output_bam,
            no_sort_output,
            subset,
            position,
            mapping_quality,
            multimapping_detection_method: _,
            buffer_whole_contig,
            output_unmapped,
            paired,
            chimeric_pairs,
            unmapped,
            per_gene,
            gene_tag,
            skip_tags_regex,
            per_contig,
            common: _,
        }) => run_group_cmd(
            input.as_deref(),
            &method,
            ignore_umi,
            output.as_deref(),
            &input_format,
            &output_format,
            random_seed,
            &umi_separator,
            chrom.as_deref(),
            group_out.as_deref(),
            output_bam,
            no_sort_output,
            subset,
            position.options(),
            mapping_quality,
            buffer_whole_contig,
            output_unmapped,
            paired,
            chimeric_pairs.as_deref(),
            &unmapped,
            per_gene,
            gene_tag.as_deref(),
            skip_tags_regex.as_deref(),
            per_contig,
        ),
        Commands::Dedup(DedupArgs {
            input,
            method,
            ignore_umi,
            output,
            input_format,
            output_format,
            random_seed,
            umi_separator,
            chrom,
            edit_distance_threshold,
            position,
            mapping_quality,
            multimapping_detection_method,
            buffer_whole_contig,
            subset,
            extract_umi_method,
            umi_tag,
            per_gene,
            gene_tag,
            skip_tags_regex,
            output_stats,
            paired,
            ignore_tlen,
            filter_umi,
            umi_whitelist,
            umi_whitelist_paired,
            common: _,
        }) => run_dedup_cmd(
            input.as_deref(),
            &method,
            ignore_umi,
            output.as_deref(),
            &input_format,
            &output_format,
            random_seed,
            &umi_separator,
            chrom.as_deref(),
            edit_distance_threshold,
            position.options(),
            mapping_quality,
            multimapping_detection_method.as_deref(),
            buffer_whole_contig,
            subset,
            &extract_umi_method,
            umi_tag.as_deref(),
            per_gene,
            gene_tag.as_deref(),
            skip_tags_regex.as_deref(),
            output_stats.as_deref(),
            paired,
            ignore_tlen,
            filter_umi,
            umi_whitelist.as_deref(),
            umi_whitelist_paired.as_deref(),
        ),
        Commands::Count(CountArgs {
            input,
            output,
            input_format,
            mapping_quality,
            method,
            gene_tag,
            skip_tags_regex,
            extract_umi_method: _,
            per_cell,
            wide_format,
            edit_distance_threshold,
            _random_seed: _,
            common,
        }) => run_count_cmd(
            input.as_deref(),
            output.as_deref(),
            &input_format,
            mapping_quality,
            &method,
            &gene_tag,
            skip_tags_regex.as_deref(),
            per_cell,
            wide_format,
            edit_distance_threshold,
            common.compresslevel,
        ),
        Commands::CountTab(CountTabArgs {
            input,
            output,
            per_cell,
            separator,
            method,
            edit_distance_threshold,
            _random_seed: _,
            common,
        }) => run_count_tab_cmd(
            input.as_deref(),
            output.as_deref(),
            per_cell,
            &separator,
            &method,
            edit_distance_threshold,
            common.compresslevel,
        ),
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

fn open_output(path: Option<&str>, compresslevel: u32) -> Result<Box<dyn Write>> {
    match path {
        Some(p) => {
            let file =
                File::create(p).with_context(|| format!("failed to create output file: {p}"))?;
            if is_gzipped(p) {
                Ok(Box::new(GzEncoder::new(
                    file,
                    Compression::new(compresslevel),
                )))
            } else {
                Ok(Box::new(file))
            }
        }
        None => Ok(Box::new(io::stdout().lock())),
    }
}

#[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
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
    ignore_read_pair_suffixes: bool,
    reconcile_pairs: bool,
    error_correct_cell: bool,
    blacklist_path: Option<&str>,
    filtered_out_path: Option<&str>,
    filtered_out2_path: Option<&str>,
    either_read: bool,
    compresslevel: u32,
) -> Result<String> {
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

    let (whitelist, correction_map) = if let Some(wl_path) = whitelist_path {
        let (wl, cm) = load_whitelist_with_correction(wl_path, error_correct_cell)?;
        (Some(wl), cm)
    } else {
        (None, None)
    };

    let blacklist = blacklist_path.map(load_blacklist).transpose()?;

    let config = ExtractConfig {
        pattern,
        pattern2,
        umi_separator: sep_byte,
        quality_filter_threshold,
        quality_encoding: qe,
        whitelist,
        correction_map,
        blacklist,
        ignore_read_pair_suffixes,
        reconcile_pairs,
    };

    let reader1 = open_input(input_path)?;

    let stats = if let Some(r2_path) = read2_in_path {
        let reader2 = open_input(Some(r2_path))?;
        if either_read {
            let writer1 = open_output(output_path, compresslevel)?;
            let writer2 = open_output(read2_out_path, compresslevel)
                .context("--read2-out is required when --either-read is specified")?;
            extract_reads_either_read(&config, reader1, reader2, writer1, writer2)
                .context("either-read extraction failed")?
        } else if read2_stdout {
            let writer = open_output(output_path, compresslevel)?;
            let filt1 = filtered_out_path
                .map(|p| open_output(Some(p), compresslevel))
                .transpose()
                .context("failed to open --filtered-out")?;
            let filt2 = filtered_out2_path
                .map(|p| open_output(Some(p), compresslevel))
                .transpose()
                .context("failed to open --filtered-out2")?;
            extract_reads_paired_r1_pattern(&config, reader1, reader2, writer, filt1, filt2)
                .context("paired-end extraction failed")?
        } else {
            let writer1 = open_output(output_path, compresslevel)?;
            let writer2 = open_output(read2_out_path, compresslevel)
                .context("--read2-out is required when --read2-in is specified")?;
            extract_reads_paired(&config, reader1, reader2, writer1, writer2)
                .context("paired-end extraction failed")?
        }
    } else {
        let writer1 = open_output(output_path, compresslevel)?;
        extract_reads(&config, reader1, writer1).context("extraction failed")?
    };

    Ok(format!(
        "Reads input: {}, output: {}, too short: {}, no match: {}, quality filtered: {}, whitelist filtered: {}",
        stats.input_reads,
        stats.output_reads,
        stats.too_short,
        stats.no_match,
        stats.quality_filtered,
        stats.whitelist_filtered,
    ))
}

type WhitelistWithCorrection = (HashSet<Vec<u8>>, Option<HashMap<Vec<u8>, Vec<u8>>>);

fn load_whitelist_with_correction(
    path: &str,
    error_correct: bool,
) -> Result<WhitelistWithCorrection> {
    let file =
        File::open(path).with_context(|| format!("failed to open whitelist file: {path}"))?;
    let reader = io::BufReader::new(file);
    let mut whitelist = HashSet::new();
    let mut correction_map = HashMap::new();

    for line in reader.lines() {
        let line = line.with_context(|| format!("failed to read whitelist file: {path}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut cols = trimmed.split('\t');
        let barcode = cols.next().unwrap().as_bytes().to_vec();
        whitelist.insert(barcode.clone());

        if error_correct
            && let Some(variants_col) = cols.next()
            && !variants_col.is_empty()
        {
            for variant in variants_col.split(',') {
                let v = variant.trim();
                if !v.is_empty() {
                    correction_map.insert(v.as_bytes().to_vec(), barcode.clone());
                }
            }
        }
    }

    let cm = if error_correct && !correction_map.is_empty() {
        Some(correction_map)
    } else {
        None
    };
    Ok((whitelist, cm))
}

fn load_blacklist(path: &str) -> Result<HashSet<Vec<u8>>> {
    let file =
        File::open(path).with_context(|| format!("failed to open blacklist file: {path}"))?;
    let reader = io::BufReader::new(file);
    let mut set = HashSet::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("failed to read blacklist file: {path}"))?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            let barcode = trimmed.split('\t').next().unwrap();
            set.insert(barcode.as_bytes().to_vec());
        }
    }
    Ok(set)
}

#[allow(clippy::too_many_arguments)]
fn run_whitelist_cmd(
    bc_pattern: &str,
    extract_method: &str,
    input_path: Option<&str>,
    output_path: Option<&str>,
    prime3: bool,
    knee_method: &str,
    set_cell_number: Option<usize>,
    expect_cells: Option<usize>,
    error_correct_threshold: usize,
    ed_above_threshold: Option<&str>,
    filtered_out_path: Option<&str>,
    subset_reads: usize,
    compresslevel: u32,
) -> Result<String> {
    let pattern = parse_pattern(bc_pattern, extract_method, prime3)?;

    let km = match knee_method {
        "distance" => KneeMethod::Distance,
        "density" => KneeMethod::Density,
        other => bail!("unknown knee method '{other}'; expected 'distance' or 'density'"),
    };

    let ed_above = ed_above_threshold
        .map(|s| match s {
            "discard" => Ok(EdAboveThreshold::Discard),
            "correct" => Ok(EdAboveThreshold::Correct),
            other => {
                bail!("unknown --ed-above-threshold '{other}'; expected 'discard' or 'correct'")
            }
        })
        .transpose()?;

    let config = WhitelistConfig {
        pattern,
        knee_method: km,
        cell_number: set_cell_number,
        expect_cells,
        error_correct_threshold,
        ed_above_threshold: ed_above,
        subset_reads,
    };

    let reader = open_input(input_path)?;
    let writer = open_output(output_path, compresslevel)?;
    let filt_out = filtered_out_path
        .map(|p| open_output(Some(p), compresslevel))
        .transpose()
        .context("failed to open --filtered-out")?;

    let stats =
        run_whitelist(&config, reader, writer, filt_out).context("whitelist command failed")?;

    Ok(format!(
        "Reads input: {}, no barcode match: {}",
        stats.input_reads, stats.no_match,
    ))
}

#[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
fn run_group_cmd(
    input_path: Option<&str>,
    method: &str,
    ignore_umi: bool,
    output_path: Option<&str>,
    input_format: &InputFormatArgs,
    output_format: &OutputFormatArgs,
    random_seed: u64,
    umi_separator: &str,
    chrom: Option<&str>,
    group_out: Option<&str>,
    output_bam: bool,
    no_sort_output: bool,
    subset: Option<f32>,
    position: PositionOptions,
    mapping_quality: u8,
    buffer_whole_contig: bool,
    output_unmapped: bool,
    paired: bool,
    chimeric_pairs: Option<&str>,
    unmapped: &str,
    per_gene: bool,
    gene_tag: Option<&str>,
    skip_tags_regex: Option<&str>,
    per_contig: bool,
) -> Result<String> {
    let input = input_path.context("--stdin is required for group (BAM input path)")?;
    if output_path.is_some() && !output_bam {
        bail!("--stdout requires --output-bam");
    }
    input_format.note_ignored_flags();
    output_format.note_ignored_flags();

    let dedup_method = match method {
        "unique" => DedupMethod::Unique,
        "percentile" => DedupMethod::Percentile,
        "cluster" => DedupMethod::Cluster,
        "adjacency" => DedupMethod::Adjacency,
        "directional" => DedupMethod::Directional,
        other => bail!("unknown method '{other}'"),
    };

    let sep_byte = umi_separator.as_bytes().first().copied().unwrap_or(b'_');

    let chimeric = match chimeric_pairs {
        Some("discard") => ChimericPairs::Discard,
        Some("output") => ChimericPairs::Output,
        Some("use") | None => ChimericPairs::Use,
        Some(other) => {
            bail!("unknown --chimeric-pairs '{other}'; expected 'discard', 'output', or 'use'")
        }
    };

    let unmapped_handling = if output_unmapped {
        UnmappedHandling::Output
    } else {
        match unmapped {
            "discard" => UnmappedHandling::Discard,
            "output" => UnmappedHandling::Output,
            "use" => UnmappedHandling::Use,
            other => {
                bail!("unknown --unmapped-reads '{other}'; expected 'discard', 'output', or 'use'")
            }
        }
    };

    let config = GroupConfig {
        method: dedup_method,
        ignore_umi,
        umi_separator: sep_byte,
        random_seed,
        output_path: output_path.map(String::from),
        output_format: output_format.resolve(output_path)?,
        reference: input_format.reference_filename.clone(),
        output_bam,
        no_sort_output,
        chrom: chrom.map(String::from),
        group_out: group_out.map(String::from),
        edit_distance_threshold: 1,
        position,
        subset,
        mapping_quality,
        buffer_whole_contig,
        per_gene,
        gene_tag: gene_tag.map(String::from),
        skip_tags_regex: skip_tags_regex.map(String::from),
        per_contig,
        paired,
        chimeric_pairs: chimeric,
        unmapped_handling,
    };

    let stats = run_group(&config, input).context("group failed")?;

    Ok(format!(
        "Reads input: {}, output: {}",
        stats.input_reads, stats.output_reads,
    ))
}

#[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
fn run_dedup_cmd(
    input_path: Option<&str>,
    method: &str,
    ignore_umi: bool,
    output_path: Option<&str>,
    input_format: &InputFormatArgs,
    output_format: &OutputFormatArgs,
    random_seed: u64,
    umi_separator: &str,
    chrom: Option<&str>,
    edit_distance_threshold: u32,
    position: PositionOptions,
    mapping_quality: u8,
    multimapping_detection_method: Option<&str>,
    buffer_whole_contig: bool,
    subset: Option<f32>,
    extract_umi_method: &str,
    umi_tag: Option<&str>,
    per_gene: bool,
    gene_tag: Option<&str>,
    skip_tags_regex: Option<&str>,
    output_stats: Option<&str>,
    paired: bool,
    ignore_tlen: bool,
    filter_umi: bool,
    umi_whitelist_path: Option<&str>,
    umi_whitelist_paired_path: Option<&str>,
) -> Result<String> {
    let input = input_path.context("--stdin is required for dedup (BAM input path)")?;
    input_format.note_ignored_flags();
    output_format.note_ignored_flags();

    let dedup_method = match method {
        "unique" => DedupMethod::Unique,
        "percentile" => DedupMethod::Percentile,
        "cluster" => DedupMethod::Cluster,
        "adjacency" => DedupMethod::Adjacency,
        "directional" => DedupMethod::Directional,
        other => bail!("unknown dedup method '{other}'"),
    };

    let sep_byte = umi_separator.as_bytes().first().copied().unwrap_or(b'_');

    let multimapping_detection = multimapping_detection_method
        .map(|name| {
            MultimappingDetection::parse(name)
                .ok_or_else(|| anyhow::anyhow!("unknown --multimapping-detection-method '{name}'"))
        })
        .transpose()?;

    let umi_whitelist = if filter_umi {
        let wl_path =
            umi_whitelist_path.context("--umi-whitelist is required when --filter-umi is set")?;
        Some(load_umi_whitelist(wl_path, umi_whitelist_paired_path)?)
    } else {
        None
    };

    let config = DedupConfig {
        method: dedup_method,
        ignore_umi,
        umi_separator: sep_byte,
        random_seed,
        output_path: output_path.map(String::from),
        output_format: output_format.resolve(output_path)?,
        reference: input_format.reference_filename.clone(),
        chrom: chrom.map(String::from),
        edit_distance_threshold,
        position,
        subset,
        mapping_quality,
        multimapping_detection,
        buffer_whole_contig,
        extract_umi_method: extract_umi_method.to_string(),
        umi_tag: umi_tag.map(String::from),
        per_gene,
        gene_tag: gene_tag.map(String::from),
        skip_tags_regex: skip_tags_regex.map(String::from),
        output_stats: output_stats.map(String::from),
        paired,
        ignore_tlen,
        umi_whitelist,
    };

    let stats = run_dedup(&config, input).context("dedup failed")?;

    Ok(format!(
        "Reads input: {}, output: {}, positions: {}",
        stats.input_reads, stats.output_reads, stats.positions,
    ))
}

fn load_umi_whitelist(path: &str, paired_path: Option<&str>) -> Result<HashSet<Vec<u8>>> {
    let load_barcodes = |p: &str| -> Result<Vec<Vec<u8>>> {
        let file = File::open(p).with_context(|| format!("failed to open UMI whitelist: {p}"))?;
        let reader = io::BufReader::new(file);
        let mut barcodes = Vec::new();
        for line in reader.lines() {
            let line = line.with_context(|| format!("failed to read UMI whitelist: {p}"))?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let barcode = trimmed.split('\t').next().unwrap();
            barcodes.push(barcode.as_bytes().to_vec());
        }
        Ok(barcodes)
    };

    let barcodes1 = load_barcodes(path)?;

    if let Some(p2) = paired_path {
        // Paired mode: Cartesian product of both whitelist files
        let barcodes2 = load_barcodes(p2)?;
        let product_size = barcodes1.len().saturating_mul(barcodes2.len());
        if product_size > 10_000_000 {
            bail!(
                "Cartesian product of UMI whitelists would produce {product_size} entries \
                 ({} x {}); refusing to proceed (limit: 10,000,000)",
                barcodes1.len(),
                barcodes2.len()
            );
        }
        let mut whitelist = HashSet::new();
        for b1 in &barcodes1 {
            for b2 in &barcodes2 {
                let mut combined = b1.clone();
                combined.extend_from_slice(b2);
                whitelist.insert(combined);
            }
        }
        Ok(whitelist)
    } else {
        Ok(barcodes1.into_iter().collect())
    }
}

#[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
fn run_count_cmd(
    input_path: Option<&str>,
    output_path: Option<&str>,
    input_format: &InputFormatArgs,
    mapping_quality: u8,
    method: &str,
    gene_tag: &str,
    skip_tags_regex: Option<&str>,
    per_cell: bool,
    wide_format: bool,
    edit_distance_threshold: u32,
    compresslevel: u32,
) -> Result<String> {
    let input = input_path.context("--stdin is required for count (BAM input path)")?;
    input_format.note_ignored_flags();

    let dedup_method = match method {
        "unique" => DedupMethod::Unique,
        "percentile" => DedupMethod::Percentile,
        "cluster" => DedupMethod::Cluster,
        "adjacency" => DedupMethod::Adjacency,
        "directional" => DedupMethod::Directional,
        other => bail!("unknown method '{other}'"),
    };

    let config = CountConfig {
        method: dedup_method,
        gene_tag: gene_tag.to_string(),
        skip_tags_regex: skip_tags_regex.map(String::from),
        per_cell,
        wide_format,
        edit_distance_threshold,
        reference: input_format.reference_filename.clone(),
        mapping_quality,
    };

    let mut output = open_output(output_path, compresslevel)?;
    let stats = run_count(&config, input, &mut output).context("count failed")?;

    Ok(format!(
        "Reads input: {}, counted: {}",
        stats.input_reads, stats.counted_reads,
    ))
}

fn run_count_tab_cmd(
    input_path: Option<&str>,
    output_path: Option<&str>,
    per_cell: bool,
    separator: &str,
    method: &str,
    edit_distance_threshold: u32,
    compresslevel: u32,
) -> Result<String> {
    let dedup_method = match method {
        "unique" => DedupMethod::Unique,
        "percentile" => DedupMethod::Percentile,
        "cluster" => DedupMethod::Cluster,
        "adjacency" => DedupMethod::Adjacency,
        "directional" => DedupMethod::Directional,
        other => bail!("unknown method '{other}'"),
    };

    let sep_byte = separator.as_bytes().first().copied().unwrap_or(b'_');

    let config = CountTabConfig {
        method: dedup_method,
        per_cell,
        separator: sep_byte,
        edit_distance_threshold,
    };

    let input = open_input(input_path)?;
    let mut reader = io::BufReader::new(input);
    let mut output = open_output(output_path, compresslevel)?;
    let stats = run_count_tab(&config, &mut reader, &mut output).context("count_tab failed")?;

    Ok(format!(
        "Reads input: {}, counted: {}",
        stats.input_reads, stats.counted_reads,
    ))
}

fn is_gzipped(path: &str) -> bool {
    Path::new(path)
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    use clap::error::ErrorKind;

    use super::*;

    fn parse(args: &[&str]) -> Cli {
        Cli::try_parse_from(std::iter::once("umi-tools-rs").chain(args.iter().copied()))
            .expect("arguments should parse")
    }

    #[test]
    fn long_flags_accept_unambiguous_prefixes() {
        let cli = parse(&["dedup", "--stdin=in.bam", "--reference-file=ref.fa"]);
        let Commands::Dedup(args) = cli.command else {
            panic!("expected dedup");
        };
        assert_eq!(
            args.input_format.reference_filename.as_deref(),
            Some("ref.fa")
        );
    }

    #[test]
    fn exact_flag_wins_over_longer_flag_with_same_prefix() {
        let cli = parse(&["extract", "--bc-pattern=NNN", "--filtered-out=a.fq"]);
        let Commands::Extract(args) = cli.command else {
            panic!("expected extract");
        };
        assert_eq!(args.filtered_out.as_deref(), Some("a.fq"));
        assert!(args.filtered_out2.is_none());
    }

    #[test]
    fn unmapped_prefix_selects_unmapped_reads() {
        let cli = parse(&["group", "--stdin=in.bam", "--unmapped=use"]);
        let Commands::Group(args) = cli.command else {
            panic!("expected group");
        };
        assert_eq!(args.unmapped, "use");
    }

    #[test]
    fn log_file_gets_header_summary_and_footer_and_is_appended_to() {
        let path =
            std::env::temp_dir().join(format!("umi-tools-rs-runlog-{}.log", std::process::id()));
        let _ = std::fs::remove_file(&path);
        let path_str = path.to_str().unwrap();
        let args = ["dedup".to_string(), "--stdin=in.bam".to_string()];

        for _ in 0..2 {
            let log = RunLog::open(Some(path_str), &args, false).unwrap();
            log.finish("Reads input: 1, output: 1").unwrap();
        }

        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 10, "two runs append two blocks:\n{text}");
        assert!(lines[0].starts_with("# umi-tools-rs version: "));
        assert_eq!(lines[1], "# output generated by dedup --stdin=in.bam");
        assert!(lines[2].starts_with("# job started at "));
        assert_eq!(lines[3], "Reads input: 1, output: 1");
        assert!(lines[4].starts_with("# job finished at "));
        assert_eq!(&lines[5..10], &lines[0..5]);
    }

    #[test]
    fn quiet_log_creates_an_empty_file() {
        let path =
            std::env::temp_dir().join(format!("umi-tools-rs-quietlog-{}.log", std::process::id()));
        let _ = std::fs::remove_file(&path);
        let log = RunLog::open(Some(path.to_str().unwrap()), &[], true).unwrap();
        log.finish("summary").unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(text, "");
    }

    #[test]
    fn verbose_and_error_are_shared_by_all_commands() {
        let cli = parse(&["extract", "--bc-pattern=NNN", "-v", "0", "-E", "err.txt"]);
        let common = cli.command.common();
        assert_eq!(common.verbose, 0);
        assert_eq!(common.error.as_deref(), Some("err.txt"));
    }

    #[test]
    fn compresslevel_is_range_checked() {
        assert!(Cli::try_parse_from(["umi-tools-rs", "extract", "--compresslevel=0"]).is_err());
        let cli = parse(&["extract", "--bc-pattern=NNN", "--compresslevel=9"]);
        assert_eq!(cli.command.common().compresslevel, 9);
    }

    #[test]
    fn help_extended_and_version_work_on_subcommands() {
        Cli::command().debug_assert();
        let err = Cli::try_parse_from(["umi-tools-rs", "dedup", "--help-extended"])
            .err()
            .expect("help stops parsing");
        assert_eq!(err.kind(), ErrorKind::DisplayHelp);
        let err = Cli::try_parse_from(["umi-tools-rs", "dedup", "--version"])
            .err()
            .expect("version stops parsing");
        assert_eq!(err.kind(), ErrorKind::DisplayVersion);
    }

    #[test]
    fn upstream_profiling_and_seed_flags_parse_everywhere() {
        let required = |command: &str| {
            if matches!(command, "extract" | "whitelist") {
                "--bc-pattern=NNN"
            } else {
                "--stdin=in.bam"
            }
        };
        for command in [
            "extract",
            "whitelist",
            "group",
            "dedup",
            "count",
            "count_tab",
        ] {
            let cli = parse(&[
                command,
                required(command),
                "--temp-dir=/tmp",
                "--timeit=t.tsv",
                "--timeit-name=x",
                "--timeit-header",
                "--timeit-header",
            ]);
            assert_eq!(cli.command.common().timeit.as_deref(), Some("t.tsv"));
        }
        for command in ["extract", "whitelist", "count", "count_tab"] {
            parse(&[command, required(command), "--random-seed=1"]);
        }
    }

    #[test]
    fn ambiguous_prefix_is_rejected() {
        assert!(Cli::try_parse_from(["umi-tools-rs", "dedup", "--std=x"]).is_err());
    }

    #[test]
    fn out_sam_may_be_repeated() {
        let cli = parse(&["dedup", "--stdin=in.bam", "--out-sam", "--out-sam"]);
        let Commands::Dedup(args) = cli.command else {
            panic!("expected dedup");
        };
        assert!(args.output_format.out_sam);
    }
}
