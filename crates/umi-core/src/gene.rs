//! Per-gene assignment of reads, mirroring the `per_gene` branch of `umi_tools`' `get_bundles`.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use regex::bytes::Regex;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{HeaderView, Record};

pub const DEFAULT_SKIP_REGEX: &str = "^(__|Unassigned)";

/// `--per-gene` and its companions.
#[derive(Debug, Clone)]
pub struct GeneOptions {
    pub per_gene: bool,
    pub per_contig: bool,
    pub gene_tag: Option<Vec<u8>>,
    /// Tag holding the assignment status; `None` means the gene tag itself.
    pub assigned_tag: Option<Vec<u8>>,
    pub skip_regex: String,
    pub transcript_map: Option<String>,
}

impl Default for GeneOptions {
    fn default() -> Self {
        Self {
            per_gene: false,
            per_contig: false,
            gene_tag: None,
            assigned_tag: None,
            skip_regex: DEFAULT_SKIP_REGEX.to_string(),
            transcript_map: None,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum GeneError {
    #[error("need to use --per-gene with --per-contig")]
    PerContigRequiresPerGene,
    #[error("need to use --per-gene with --gene-tag")]
    GeneTagRequiresPerGene,
    #[error("for per-gene applications, must supply --per-contig or --gene-tag")]
    PerGeneNeedsSource,
    #[error("need to use either --per-contig OR --gene-tag, please do not provide both")]
    BothSources,
    #[error("need to use --per-contig and --per-gene with --gene-transcript-map")]
    TranscriptMapRequiresPerContig,
    #[error("skip-regex '{0}' is not a valid regex")]
    InvalidRegex(String),
    #[error("failed to read --gene-transcript-map {0}: {1}")]
    TranscriptMapRead(String, String),
    #[error("--gene-transcript-map line does not have two tab-separated fields: {0}")]
    TranscriptMapLine(String),
}

impl GeneOptions {
    /// `umi_tools.Utilities.validateSamOptions`.
    ///
    /// # Errors
    ///
    /// Returns the same complaints `umi_tools` raises for inconsistent options.
    pub fn validate(&self) -> Result<(), GeneError> {
        if self.per_gene {
            if self.gene_tag.is_some() && self.per_contig {
                return Err(GeneError::BothSources);
            }
            if !self.per_contig && self.gene_tag.is_none() {
                return Err(GeneError::PerGeneNeedsSource);
            }
        }
        if self.per_contig && !self.per_gene {
            return Err(GeneError::PerContigRequiresPerGene);
        }
        if self.gene_tag.is_some() && !self.per_gene {
            return Err(GeneError::GeneTagRequiresPerGene);
        }
        if self.transcript_map.is_some() && !self.per_contig {
            return Err(GeneError::TranscriptMapRequiresPerContig);
        }
        Regex::new(&self.skip_regex)
            .map_err(|_| GeneError::InvalidRegex(self.skip_regex.clone()))?;
        Ok(())
    }
}

/// Genes and their transcripts from `--gene-transcript-map`, in file order and
/// restricted to contigs present in the BAM header.
///
/// `umi_tools` keeps each gene's transcripts in a `set`, so its contig order (and
/// therefore which read represents a UMI group) changes with Python's hash seed.
/// Here transcripts keep file order.
#[derive(Debug, Clone, Default)]
pub struct TranscriptMap {
    pub genes: Vec<(Vec<u8>, Vec<u32>)>,
    gene_index: HashMap<Vec<u8>, usize>,
    gene_of_tid: HashMap<u32, usize>,
}

impl TranscriptMap {
    /// # Errors
    ///
    /// Returns an error when the file cannot be read or a line lacks two fields.
    pub fn load(path: &str, header: &HeaderView) -> Result<Self, GeneError> {
        let read_error =
            |e: std::io::Error| GeneError::TranscriptMapRead(path.to_string(), e.to_string());
        let file = File::open(path).map_err(read_error)?;
        let mut map = Self::default();
        for line in BufReader::new(file).lines() {
            let line = line.map_err(read_error)?;
            if line.starts_with('#') {
                continue;
            }
            // umi_tools stops at the first blank line.
            if line.trim().is_empty() {
                break;
            }
            let mut fields = line.trim().split('\t');
            let (Some(gene), Some(transcript), None) =
                (fields.next(), fields.next(), fields.next())
            else {
                return Err(GeneError::TranscriptMapLine(line));
            };
            let Some(tid) = header.tid(transcript.as_bytes()) else {
                continue;
            };
            map.add(gene.as_bytes(), tid);
        }
        Ok(map)
    }

    fn add(&mut self, gene: &[u8], tid: u32) {
        let index = *self.gene_index.entry(gene.to_vec()).or_insert_with(|| {
            self.genes.push((gene.to_vec(), Vec::new()));
            self.genes.len() - 1
        });
        if !self.genes[index].1.contains(&tid) {
            self.genes[index].1.push(tid);
        }
        self.gene_of_tid.insert(tid, index);
    }

    #[must_use]
    pub fn gene_of_tid(&self, tid: u32) -> Option<&[u8]> {
        self.gene_of_tid
            .get(&tid)
            .map(|&index| self.genes[index].0.as_slice())
    }

    fn transcript_count(&self, gene: &[u8]) -> Option<usize> {
        self.gene_index
            .get(gene)
            .map(|&index| self.genes[index].1.len())
    }
}

/// Decides which gene each read is grouped under.
pub struct GeneAssigner {
    per_contig: bool,
    gene_tag: Option<Vec<u8>>,
    assigned_tag: Option<Vec<u8>>,
    skip_regex: Regex,
    map: Option<TranscriptMap>,
    header: HeaderView,
}

impl GeneAssigner {
    /// `None` unless `--per-gene` is set.
    ///
    /// # Errors
    ///
    /// Returns option validation and transcript-map errors.
    pub fn new(options: &GeneOptions, header: &HeaderView) -> Result<Option<Self>, GeneError> {
        options.validate()?;
        if !options.per_gene {
            return Ok(None);
        }
        let map = options
            .transcript_map
            .as_deref()
            .map(|path| TranscriptMap::load(path, header))
            .transpose()?;
        Ok(Some(Self {
            per_contig: options.per_contig,
            gene_tag: options.gene_tag.clone(),
            assigned_tag: options.assigned_tag.clone(),
            skip_regex: Regex::new(&options.skip_regex)
                .map_err(|_| GeneError::InvalidRegex(options.skip_regex.clone()))?,
            map,
            header: header.clone(),
        }))
    }

    #[must_use]
    pub const fn transcript_map(&self) -> Option<&TranscriptMap> {
        self.map.as_ref()
    }

    #[must_use]
    pub const fn per_contig(&self) -> bool {
        self.per_contig
    }

    /// Gene to group the read under, or `None` when `umi_tools` skips the read:
    /// a missing gene or status tag, an empty gene, or a status matching the skip regex.
    #[must_use]
    pub fn gene(&self, record: &Record) -> Option<Vec<u8>> {
        if self.per_contig {
            let tid = u32::try_from(record.tid()).ok()?;
            return Some(match &self.map {
                Some(map) => map.gene_of_tid(tid)?.to_vec(),
                None => self.header.tid2name(tid).to_vec(),
            });
        }
        let gene_tag = self.gene_tag.as_deref()?;
        let assigned_tag = self.assigned_tag.as_deref().unwrap_or(gene_tag);
        let assigned = string_tag(record, assigned_tag)?;
        let gene = string_tag(record, gene_tag)?;
        if gene.is_empty() || self.skip_regex.is_match(&assigned) {
            return None;
        }
        Some(gene)
    }

    #[must_use]
    pub fn flusher(&self) -> GeneFlusher {
        GeneFlusher::new(self.map.as_ref())
    }
}

fn string_tag(record: &Record, tag: &[u8]) -> Option<Vec<u8>> {
    match record.aux(tag).ok()? {
        Aux::String(value) => Some(value.as_bytes().to_vec()),
        Aux::Char(value) => Some(vec![value]),
        _ => None,
    }
}

/// Which buffered genes to flush before a read is added.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Flush {
    None,
    /// Every buffered gene, in name order.
    All,
    /// One gene whose transcripts have all been seen.
    Gene(Vec<u8>),
}

/// `get_bundles.check_output` for per-gene mode.
///
/// Flush everything when the contig changes, or with a transcript map flush the
/// previous gene once the contig changes and all of that gene's transcripts have
/// contributed reads.
pub struct GeneFlusher {
    transcripts_per_gene: Option<HashMap<Vec<u8>, usize>>,
    observed: HashMap<Vec<u8>, HashSet<i32>>,
    last: Option<(i32, Vec<u8>)>,
}

impl GeneFlusher {
    fn new(map: Option<&TranscriptMap>) -> Self {
        Self {
            transcripts_per_gene: map.map(|map| {
                map.genes
                    .iter()
                    .map(|(gene, _)| (gene.clone(), map.transcript_count(gene).unwrap_or(0)))
                    .collect()
            }),
            observed: HashMap::new(),
            last: None,
        }
    }

    #[must_use]
    pub fn before_read(&self, tid: i32) -> Flush {
        let Some((last_tid, last_gene)) = &self.last else {
            return Flush::None;
        };
        if tid == *last_tid {
            return Flush::None;
        }
        let Some(counts) = &self.transcripts_per_gene else {
            return Flush::All;
        };
        let observed = self.observed.get(last_gene).map_or(0, HashSet::len);
        if counts.get(last_gene) == Some(&observed) {
            Flush::Gene(last_gene.clone())
        } else {
            Flush::None
        }
    }

    pub fn after_read(&mut self, tid: i32, gene: &[u8]) {
        if self.transcripts_per_gene.is_some() {
            self.observed.entry(gene.to_vec()).or_default().insert(tid);
        }
        self.last = Some((tid, gene.to_vec()));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> HeaderView {
        HeaderView::from_bytes(
            b"@SQ\tSN:tx1\tLN:1000\n@SQ\tSN:tx2\tLN:1000\n@SQ\tSN:tx3\tLN:1000\n",
        )
    }

    fn read(tid: i32, tags: &str) -> Record {
        let contig = ["tx1", "tx2", "tx3"][usize::try_from(tid).unwrap()];
        let line = format!("r\t0\t{contig}\t101\t60\t10M\t*\t0\t0\t*\t*{tags}");
        Record::from_sam(&header(), line.as_bytes()).unwrap()
    }

    fn tagged_options() -> GeneOptions {
        GeneOptions {
            per_gene: true,
            gene_tag: Some(b"XT".to_vec()),
            ..GeneOptions::default()
        }
    }

    #[test]
    fn validation_matches_umi_tools() {
        assert!(GeneOptions::default().validate().is_ok());
        let per_gene_only = GeneOptions {
            per_gene: true,
            ..GeneOptions::default()
        };
        assert!(matches!(
            per_gene_only.validate(),
            Err(GeneError::PerGeneNeedsSource)
        ));
        let both = GeneOptions {
            per_contig: true,
            ..tagged_options()
        };
        assert!(matches!(both.validate(), Err(GeneError::BothSources)));
        let contig_only = GeneOptions {
            per_contig: true,
            ..GeneOptions::default()
        };
        assert!(matches!(
            contig_only.validate(),
            Err(GeneError::PerContigRequiresPerGene)
        ));
        let map_without_contig = GeneOptions {
            transcript_map: Some("map.tsv".into()),
            ..tagged_options()
        };
        assert!(matches!(
            map_without_contig.validate(),
            Err(GeneError::TranscriptMapRequiresPerContig)
        ));
    }

    #[test]
    fn gene_tag_skips_by_status_and_empty_gene() {
        let assigner = GeneAssigner::new(&tagged_options(), &header())
            .unwrap()
            .unwrap();
        assert_eq!(
            assigner.gene(&read(0, "\tXT:Z:ENSG1")),
            Some(b"ENSG1".to_vec())
        );
        assert_eq!(
            assigner.gene(&read(0, "\tXT:Z:Unassigned_NoFeatures")),
            None
        );
        assert_eq!(assigner.gene(&read(0, "\tXT:Z:__ambiguous")), None);
        assert_eq!(assigner.gene(&read(0, "\tXT:Z:")), None, "empty gene");
        assert_eq!(assigner.gene(&read(0, "\tNH:i:1")), None, "missing tag");

        let with_status = GeneAssigner::new(
            &GeneOptions {
                assigned_tag: Some(b"XS".to_vec()),
                ..tagged_options()
            },
            &header(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(
            with_status.gene(&read(0, "\tXT:Z:ENSG1\tXS:Z:Assigned")),
            Some(b"ENSG1".to_vec())
        );
        assert_eq!(
            with_status.gene(&read(0, "\tXT:Z:ENSG1\tXS:Z:Unassigned_MultiMapping")),
            None
        );
    }

    #[test]
    fn per_contig_uses_contig_or_mapped_gene() {
        let contig = GeneAssigner::new(
            &GeneOptions {
                per_gene: true,
                per_contig: true,
                ..GeneOptions::default()
            },
            &header(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(contig.gene(&read(1, "")), Some(b"tx2".to_vec()));

        let mut map = TranscriptMap::default();
        map.add(b"geneA", 0);
        map.add(b"geneA", 1);
        map.add(b"geneB", 2);
        assert_eq!(map.gene_of_tid(1), Some(&b"geneA"[..]));
        assert_eq!(map.transcript_count(b"geneA"), Some(2));
        assert_eq!(map.genes.len(), 2);
    }

    #[test]
    fn flusher_follows_check_output() {
        let mut plain = GeneFlusher::new(None);
        assert_eq!(
            plain.before_read(0),
            Flush::None,
            "first read never flushes"
        );
        plain.after_read(0, b"g1");
        assert_eq!(plain.before_read(0), Flush::None);
        assert_eq!(plain.before_read(1), Flush::All);

        let mut map = TranscriptMap::default();
        map.add(b"geneA", 0);
        map.add(b"geneA", 1);
        map.add(b"geneB", 2);
        let mut meta = GeneFlusher::new(Some(&map));
        meta.after_read(0, b"geneA");
        assert_eq!(
            meta.before_read(1),
            Flush::None,
            "geneA still has an unseen transcript"
        );
        meta.after_read(1, b"geneA");
        assert_eq!(meta.before_read(2), Flush::Gene(b"geneA".to_vec()));
    }
}
