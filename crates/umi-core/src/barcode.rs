//! UMI and cell barcode extraction, mirroring `umi_tools.sam_methods`.

use rust_htslib::bam::Record;
use rust_htslib::bam::record::Aux;

/// Where a read's UMI and cell barcode are encoded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BarcodeSource {
    /// The last `separator`-delimited field of the read name is the UMI and
    /// the one before it the cell barcode.
    ReadId { separator: Vec<u8> },
    /// Aux tags, optionally cut at `split` or stripped of `delimiter`.
    Tag {
        umi_tag: Vec<u8>,
        umi_split: Option<Vec<u8>>,
        umi_delimiter: Option<Vec<u8>>,
        cell_tag: Option<Vec<u8>>,
        cell_split: Option<Vec<u8>>,
        cell_delimiter: Option<Vec<u8>>,
    },
    /// `UMI_` and `CELL_` fields of a colon-separated read name, as written by umis.
    Umis,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BarcodeExtractor {
    pub source: BarcodeSource,
    pub per_cell: bool,
}

impl Default for BarcodeExtractor {
    fn default() -> Self {
        Self {
            source: BarcodeSource::ReadId {
                separator: b"_".to_vec(),
            },
            per_cell: false,
        }
    }
}

/// A read's UMI and, with `per_cell`, its cell barcode. `cell` is empty otherwise.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Barcode {
    pub umi: Vec<u8>,
    pub cell: Vec<u8>,
}

#[derive(Debug, thiserror::Error)]
pub enum BarcodeError {
    /// `umi_tools` skips the read.
    #[error("read {0} lacks the UMI or cell tag")]
    MissingTag(String),
    /// `umi_tools` aborts.
    #[error("could not extract the UMI or cell barcode from read {0}")]
    Malformed(String),
}

impl BarcodeExtractor {
    /// # Errors
    ///
    /// `MissingTag` when a tag is absent, `Malformed` when the name or tag
    /// cannot hold a barcode.
    pub fn extract(&self, record: &Record) -> Result<Barcode, BarcodeError> {
        match &self.source {
            BarcodeSource::ReadId { separator } => {
                let fields = split_all(record.qname(), separator);
                let umi = fields[fields.len() - 1].to_vec();
                let cell = if self.per_cell {
                    fields
                        .len()
                        .checked_sub(2)
                        .map(|i| fields[i].to_vec())
                        .ok_or_else(|| malformed(record))?
                } else {
                    Vec::new()
                };
                Ok(Barcode { umi, cell })
            }
            BarcodeSource::Tag {
                umi_tag,
                umi_split,
                umi_delimiter,
                cell_tag,
                cell_split,
                cell_delimiter,
            } => {
                let mut umi = string_tag(record, umi_tag)?;
                if let Some(split) = umi_split {
                    umi = first_field(&umi, split);
                }
                if let Some(delimiter) = umi_delimiter {
                    umi = remove_all(&umi, delimiter);
                }
                let mut cell = Vec::new();
                if self.per_cell {
                    let tag = cell_tag.as_ref().ok_or_else(|| malformed(record))?;
                    cell = string_tag(record, tag)?;
                    if !cell.is_empty() {
                        if let Some(split) = cell_split {
                            cell = first_field(&cell, split);
                        }
                        if let Some(delimiter) = cell_delimiter {
                            cell = remove_all(&cell, delimiter);
                        }
                    }
                }
                Ok(Barcode { umi, cell })
            }
            BarcodeSource::Umis => {
                let mut umi = None;
                let mut cell = None;
                for element in record.qname().split(|&b| b == b':') {
                    if let Some(value) = element.strip_prefix(b"UMI_") {
                        umi = Some(value.to_vec());
                    } else if let Some(value) = element.strip_prefix(b"CELL_")
                        && self.per_cell
                    {
                        cell = Some(value.to_vec());
                    }
                }
                Ok(Barcode {
                    umi: umi.ok_or_else(|| malformed(record))?,
                    cell: cell.unwrap_or_default(),
                })
            }
        }
    }

    /// Barcode used for grouping, or `None` when `umi_tools` would skip the read.
    ///
    /// # Errors
    ///
    /// Returns the extraction error `umi_tools` would abort on.
    pub fn for_grouping(
        &self,
        record: &Record,
        ignore_umi: bool,
    ) -> Result<Option<Barcode>, BarcodeError> {
        if ignore_umi && !self.per_cell {
            return Ok(Some(Barcode::default()));
        }
        match self.extract(record) {
            Ok(mut barcode) => {
                if ignore_umi {
                    barcode.umi.clear();
                }
                Ok(Some(barcode))
            }
            Err(BarcodeError::MissingTag(_)) if !ignore_umi => Ok(None),
            Err(e) => Err(e),
        }
    }
}

fn malformed(record: &Record) -> BarcodeError {
    BarcodeError::Malformed(String::from_utf8_lossy(record.qname()).into_owned())
}

fn string_tag(record: &Record, tag: &[u8]) -> Result<Vec<u8>, BarcodeError> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(value.as_bytes().to_vec()),
        Ok(Aux::Char(value)) => Ok(vec![value]),
        Ok(_) => Err(malformed(record)),
        Err(_) => Err(BarcodeError::MissingTag(
            String::from_utf8_lossy(record.qname()).into_owned(),
        )),
    }
}

/// Python's `str.split(sep)`: every field, including empty ones.
fn split_all<'a>(haystack: &'a [u8], needle: &[u8]) -> Vec<&'a [u8]> {
    if needle.is_empty() {
        return vec![haystack];
    }
    let mut fields = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i + needle.len() <= haystack.len() {
        if &haystack[i..i + needle.len()] == needle {
            fields.push(&haystack[start..i]);
            i += needle.len();
            start = i;
        } else {
            i += 1;
        }
    }
    fields.push(&haystack[start..]);
    fields
}

fn first_field(value: &[u8], separator: &[u8]) -> Vec<u8> {
    split_all(value, separator)[0].to_vec()
}

fn remove_all(value: &[u8], delimiter: &[u8]) -> Vec<u8> {
    split_all(value, delimiter).concat()
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::HeaderView;

    use super::*;

    fn read(name: &str, tags: &str) -> Record {
        let header = HeaderView::from_bytes(b"@SQ\tSN:chr1\tLN:100000\n");
        let line = format!("{name}\t0\tchr1\t101\t60\t10M\t*\t0\t0\t*\t*{tags}");
        Record::from_sam(&header, line.as_bytes()).unwrap()
    }

    fn tag_source(umi_split: &str, umi_delimiter: &str, cell_split: &str) -> BarcodeSource {
        let opt = |s: &str| (!s.is_empty()).then(|| s.as_bytes().to_vec());
        BarcodeSource::Tag {
            umi_tag: b"RX".to_vec(),
            umi_split: opt(umi_split),
            umi_delimiter: opt(umi_delimiter),
            cell_tag: Some(b"CB".to_vec()),
            cell_split: opt(cell_split),
            cell_delimiter: None,
        }
    }

    #[test]
    fn read_id_takes_last_fields_with_a_multibyte_separator() {
        let extractor = BarcodeExtractor {
            source: BarcodeSource::ReadId {
                separator: b"::".to_vec(),
            },
            per_cell: true,
        };
        let barcode = extractor.extract(&read("id::CELL::UMI", "")).unwrap();
        assert_eq!(barcode.umi, b"UMI");
        assert_eq!(barcode.cell, b"CELL");
    }

    #[test]
    fn read_id_without_separator_is_the_whole_name_and_fails_per_cell() {
        let single = BarcodeExtractor::default();
        assert_eq!(single.extract(&read("plain", "")).unwrap().umi, b"plain");
        let per_cell = BarcodeExtractor {
            per_cell: true,
            ..BarcodeExtractor::default()
        };
        assert!(matches!(
            per_cell.extract(&read("plain", "")),
            Err(BarcodeError::Malformed(_))
        ));
    }

    #[test]
    fn tag_split_and_delimiter_follow_umi_tools() {
        let with_delimiter = BarcodeExtractor {
            source: tag_source("", "-", "-"),
            per_cell: true,
        };
        let barcode = with_delimiter
            .extract(&read("r", "\tRX:Z:GT-GACC\tCB:Z:ACAAGG-1"))
            .unwrap();
        assert_eq!(barcode.umi, b"GTGACC");
        assert_eq!(
            barcode.cell, b"ACAAGG",
            "GEM suffix dropped by the default cell split"
        );

        let with_split = BarcodeExtractor {
            source: tag_source("-", "", ""),
            per_cell: false,
        };
        let barcode = with_split
            .extract(&read("r", "\tRX:Z:GT-GACC\tCB:Z:ACAAGG-1"))
            .unwrap();
        assert_eq!(barcode.umi, b"GT");
        assert!(barcode.cell.is_empty());
    }

    #[test]
    fn missing_tag_skips_and_non_string_tag_aborts() {
        let extractor = BarcodeExtractor {
            source: tag_source("", "", ""),
            per_cell: false,
        };
        assert!(matches!(
            extractor.extract(&read("r", "\tNH:i:1")),
            Err(BarcodeError::MissingTag(_))
        ));
        assert!(matches!(
            extractor.extract(&read("r", "\tRX:i:7")),
            Err(BarcodeError::Malformed(_))
        ));
        assert!(
            extractor
                .for_grouping(&read("r", "\tNH:i:1"), false)
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn umis_names_yield_umi_and_optionally_cell() {
        let name = "M0:1:2:3:4:5:6:CELL_ACAAGG:UMI_GTGACC:SAMPLE_X";
        let single = BarcodeExtractor {
            source: BarcodeSource::Umis,
            per_cell: false,
        };
        let barcode = single.extract(&read(name, "")).unwrap();
        assert_eq!(barcode.umi, b"GTGACC");
        assert!(barcode.cell.is_empty());
        let per_cell = BarcodeExtractor {
            source: BarcodeSource::Umis,
            per_cell: true,
        };
        assert_eq!(per_cell.extract(&read(name, "")).unwrap().cell, b"ACAAGG");
        assert!(matches!(
            single.extract(&read("no:umi:here", "")),
            Err(BarcodeError::Malformed(_))
        ));
    }

    #[test]
    fn ignore_umi_keeps_the_cell() {
        let extractor = BarcodeExtractor {
            source: BarcodeSource::Umis,
            per_cell: true,
        };
        let barcode = extractor
            .for_grouping(&read("a:b:c:d:e:f:g:CELL_AC:UMI_GT", ""), true)
            .unwrap()
            .unwrap();
        assert!(barcode.umi.is_empty());
        assert_eq!(barcode.cell, b"AC");
    }
}
