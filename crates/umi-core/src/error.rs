use thiserror::Error;

#[derive(Debug, Error)]
pub enum ExtractError {
    #[error("invalid barcode pattern: {0}")]
    InvalidPattern(String),

    #[error("read too short ({read_len} bp) for pattern ({pattern_len} bp)")]
    ReadTooShort { read_len: usize, pattern_len: usize },

    #[error("regex did not match read sequence")]
    RegexNoMatch,

    #[error(
        "--set-cell-number ({requested}) must be smaller than the number of observed cell barcodes ({observed})"
    )]
    InvalidCellNumber { requested: usize, observed: usize },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTQ parse error: {0}")]
    FastqParse(String),

    #[error(
        "No local minima was accepted. Recommend checking the plot output and counts per local minima (requires `--plot-prefix`option) and then re-running with manually selected threshold (`--set-cell-number` option)"
    )]
    NoThreshold,
}
