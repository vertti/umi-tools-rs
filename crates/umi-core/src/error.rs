use thiserror::Error;

#[derive(Debug, Error)]
pub enum ExtractError {
    #[error("invalid barcode pattern: {0}")]
    InvalidPattern(String),

    #[error("read too short ({read_len} bp) for pattern ({pattern_len} bp)")]
    ReadTooShort { read_len: usize, pattern_len: usize },

    #[error("regex did not match read sequence")]
    RegexNoMatch,

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTQ parse error: {0}")]
    FastqParse(String),
}
