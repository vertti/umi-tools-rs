//! SAM/BAM/CRAM input and output shared by dedup, group and count.

use std::path::Path;

use rust_htslib::bam::{self, HeaderView};
use rust_htslib::errors::Error as HtsError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlignmentFormat {
    Sam,
    #[default]
    Bam,
    Cram,
}

impl AlignmentFormat {
    /// # Errors
    ///
    /// Returns the unrecognised name.
    pub fn parse(name: &str) -> Result<Self, String> {
        match name.to_ascii_lowercase().as_str() {
            "sam" => Ok(Self::Sam),
            "bam" => Ok(Self::Bam),
            "cram" => Ok(Self::Cram),
            _ => Err(name.to_string()),
        }
    }

    fn from_extension(path: &str) -> Option<Self> {
        let ext = Path::new(path).extension()?;
        if ext.eq_ignore_ascii_case("sam") {
            Some(Self::Sam)
        } else if ext.eq_ignore_ascii_case("cram") {
            Some(Self::Cram)
        } else {
            None
        }
    }

    const fn htslib(self) -> bam::Format {
        match self {
            Self::Sam => bam::Format::Sam,
            Self::Bam => bam::Format::Bam,
            Self::Cram => bam::Format::Cram,
        }
    }
}

/// Resolves the output format the way `umi_tools` does: `--out-sam` wins, then
/// `--out-format`, then a `.sam`/`.cram` file extension, otherwise BAM.
#[must_use]
pub fn determine_format(
    path: Option<&str>,
    sam_flag: bool,
    explicit: Option<AlignmentFormat>,
) -> AlignmentFormat {
    if sam_flag {
        return AlignmentFormat::Sam;
    }
    if let Some(format) = explicit {
        return format;
    }
    path.and_then(AlignmentFormat::from_extension)
        .unwrap_or_default()
}

/// Where and how alignments are written. `path` of `None` means stdout.
#[derive(Debug, Clone, Copy)]
pub struct AlignmentOutput<'a> {
    pub path: Option<&'a str>,
    pub format: AlignmentFormat,
    pub reference: Option<&'a str>,
}

/// Opens an alignment file. The format is detected from the content, and
/// `reference` is the FASTA used to decode CRAM.
///
/// # Errors
///
/// Returns htslib errors from opening the file or the reference.
pub fn open_reader(path: &str, reference: Option<&str>) -> Result<bam::Reader, HtsError> {
    let mut reader = bam::Reader::from_path(path)?;
    if let Some(reference) = reference {
        reader.set_reference(reference)?;
    }
    Ok(reader)
}

/// # Errors
///
/// Returns htslib errors from creating the file or loading the CRAM reference.
pub fn open_writer(
    header: &bam::Header,
    output: AlignmentOutput<'_>,
) -> Result<bam::Writer, HtsError> {
    let format = output.format.htslib();
    let mut writer = match output.path {
        Some(path) => bam::Writer::from_path(path, header, format)?,
        None => bam::Writer::from_stdout(header, format)?,
    };
    if output.format == AlignmentFormat::Cram
        && let Some(reference) = output.reference
    {
        writer.set_reference(reference)?;
    }
    Ok(writer)
}

/// Header for coordinate-sorted output, rewritten the way `samtools sort` does it.
#[must_use]
pub fn coordinate_sorted_header(template: &HeaderView) -> bam::Header {
    let text = coordinate_sorted_header_text(template.as_bytes());
    bam::Header::from_template(&HeaderView::from_bytes(&text))
}

/// Sets `SO:coordinate` on the `@HD` line, keeping an existing `VN` (else 1.6)
/// and dropping any other `@HD` fields, as `samtools sort` does.
fn coordinate_sorted_header_text(text: &[u8]) -> Vec<u8> {
    let mut lines: Vec<&[u8]> = text.split(|&b| b == b'\n').collect();
    if lines.last().is_some_and(|l| l.is_empty()) {
        lines.pop();
    }

    let hd_index = lines.iter().position(|l| l.starts_with(b"@HD"));
    let version = hd_index
        .and_then(|i| {
            lines[i]
                .split(|&b| b == b'\t')
                .find_map(|field| field.strip_prefix(b"VN:"))
        })
        .unwrap_or(b"1.6");

    let mut hd = b"@HD\tVN:".to_vec();
    hd.extend_from_slice(version);
    hd.extend_from_slice(b"\tSO:coordinate");

    let mut out = Vec::with_capacity(text.len() + hd.len() + 1);
    match hd_index {
        Some(i) => lines[i] = &hd,
        None => lines.insert(0, &hd),
    }
    for line in lines {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_is_case_insensitive() {
        assert_eq!(AlignmentFormat::parse("CRAM"), Ok(AlignmentFormat::Cram));
        assert_eq!(AlignmentFormat::parse("sam"), Ok(AlignmentFormat::Sam));
        assert_eq!(AlignmentFormat::parse("vcf"), Err("vcf".to_string()));
    }

    #[test]
    fn sam_flag_beats_explicit_format_and_extension() {
        assert_eq!(
            determine_format(Some("out.cram"), true, Some(AlignmentFormat::Cram)),
            AlignmentFormat::Sam
        );
    }

    #[test]
    fn explicit_format_beats_extension() {
        assert_eq!(
            determine_format(Some("out.sam"), false, Some(AlignmentFormat::Cram)),
            AlignmentFormat::Cram
        );
    }

    #[test]
    fn extension_decides_when_nothing_explicit() {
        assert_eq!(
            determine_format(Some("out.SAM"), false, None),
            AlignmentFormat::Sam
        );
        assert_eq!(
            determine_format(Some("out.cram"), false, None),
            AlignmentFormat::Cram
        );
        assert_eq!(
            determine_format(Some("out.bam"), false, None),
            AlignmentFormat::Bam
        );
        assert_eq!(
            determine_format(Some("out.txt"), false, None),
            AlignmentFormat::Bam
        );
        assert_eq!(determine_format(None, false, None), AlignmentFormat::Bam);
    }

    #[test]
    fn missing_hd_line_is_inserted_first() {
        let text = b"@SQ\tSN:chr1\tLN:1000\n@RG\tID:g\n@CO\thello\n";
        assert_eq!(
            coordinate_sorted_header_text(text),
            b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:g\n@CO\thello\n"
        );
    }

    #[test]
    fn existing_hd_keeps_version_and_drops_other_fields() {
        let text = b"@HD\tVN:1.0\tSO:unsorted\tGO:none\n@SQ\tSN:chr1\tLN:1000\n";
        assert_eq!(
            coordinate_sorted_header_text(text),
            b"@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        );
    }

    #[test]
    fn already_sorted_header_is_unchanged() {
        let text = b"@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@PG\tID:x\tPN:x\n";
        assert_eq!(coordinate_sorted_header_text(text), text);
    }

    #[test]
    fn header_without_trailing_newline_gets_one() {
        let text = b"@SQ\tSN:chr1\tLN:1000";
        assert_eq!(
            coordinate_sorted_header_text(text),
            b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
        );
    }
}
