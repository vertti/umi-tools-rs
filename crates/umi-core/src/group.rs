use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};

use rust_htslib::bam::{self, Record};

use crate::alignment_io::{self, AlignmentFormat, AlignmentOutput};
use crate::alignment_sort::RecordOutput;
use crate::barcode::{Barcode, BarcodeError, BarcodeExtractor};
use crate::clustering::cluster_umis;
use crate::dedup::{
    DedupMethod, GroupKey, PositionOptions, PythonRandom, TieBreakRng, five_prime_position,
    get_read_position,
};
use crate::gene::{Flush, GeneAssigner, GeneError, GeneOptions};
use crate::pairing::PairingOptions;

#[allow(clippy::struct_excessive_bools)]
pub struct GroupConfig {
    pub method: DedupMethod,
    pub ignore_umi: bool,
    pub barcode: BarcodeExtractor,
    pub umi_group_tag: Vec<u8>,
    pub random_seed: u64,
    pub output_path: Option<String>,
    pub output_format: AlignmentFormat,
    pub reference: Option<String>,
    pub output_bam: bool,
    pub no_sort_output: bool,
    pub chrom: Option<String>,
    pub group_out: Option<String>,
    pub edit_distance_threshold: u32,
    pub position: PositionOptions,
    pub subset: Option<f32>,
    pub mapping_quality: u8,
    pub buffer_whole_contig: bool,
    pub gene: GeneOptions,
    pub pairing: PairingOptions,
    pub ignore_tlen: bool,
}

pub struct GroupStats {
    pub input_reads: u64,
    pub output_reads: u64,
}

struct GroupSlot {
    records: Vec<Record>,
    count: u32,
    insertion_order: u32,
}

type Drained<K> = BTreeMap<K, BTreeMap<GroupKey, HashMap<Vec<u8>, GroupSlot>>>;

struct GroupBuffer<K: Ord = i64> {
    groups: Drained<K>,
    insertion_counters: BTreeMap<K, BTreeMap<GroupKey, u32>>,
}

impl GroupBuffer<i64> {
    fn drain_up_to(&mut self, threshold: i64) -> Drained<i64> {
        let rest = self.groups.split_off(&(threshold + 1));
        let drained = std::mem::replace(&mut self.groups, rest);
        let rest_counters = self.insertion_counters.split_off(&(threshold + 1));
        let _ = std::mem::replace(&mut self.insertion_counters, rest_counters);
        drained
    }
}

impl<K: Ord + Clone> GroupBuffer<K> {
    const fn new() -> Self {
        Self {
            groups: BTreeMap::new(),
            insertion_counters: BTreeMap::new(),
        }
    }

    fn add(&mut self, record: Record, pos: K, key: GroupKey, umi: Vec<u8>) {
        let umi_map = self
            .groups
            .entry(pos.clone())
            .or_default()
            .entry(key.clone())
            .or_default();

        if let Some(slot) = umi_map.get_mut(&umi) {
            slot.count += 1;
            slot.records.push(record);
            return;
        }

        let counter = self
            .insertion_counters
            .entry(pos)
            .or_default()
            .entry(key)
            .or_default();
        let order = *counter;
        *counter += 1;

        umi_map.insert(
            umi,
            GroupSlot {
                records: vec![record],
                count: 1,
                insertion_order: order,
            },
        );
    }

    fn drain_key(&mut self, key: &K) -> Drained<K> {
        self.insertion_counters.remove(key);
        self.groups
            .remove(key)
            .map(|key_map| BTreeMap::from([(key.clone(), key_map)]))
            .unwrap_or_default()
    }

    fn drain_all(&mut self) -> Drained<K> {
        let drained = std::mem::take(&mut self.groups);
        self.insertion_counters.clear();
        drained
    }
}

fn assign_groups(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, GroupSlot>,
    edit_threshold: u32,
) -> Vec<Vec<Vec<u8>>> {
    cluster_umis(
        method,
        umi_map
            .iter()
            .map(|(umi, slot)| (umi.as_slice(), slot.count, slot.insertion_order)),
        edit_threshold,
    )
    .into_iter()
    .map(|group| group.into_iter().map(<[u8]>::to_vec).collect())
    .collect()
}

/// Process drained position groups: assign UMI groups, annotate records, write TSV rows.
#[allow(clippy::cast_sign_loss)]
fn process_drained<K: Ord>(
    drained: Drained<K>,
    config: &GroupConfig,
    unique_id: &mut u32,
    tsv_writer: &mut Option<BufWriter<File>>,
    header_view: &bam::HeaderView,
    assigner: Option<&GeneAssigner>,
    output: &mut RecordOutput,
) -> Result<(), GroupError> {
    for key_map in drained.into_values() {
        for (_, mut umi_map) in key_map {
            let groups = assign_groups(config.method, &umi_map, config.edit_distance_threshold);

            for group in &groups {
                let top_umi = &group[0];
                let group_count: u32 = group.iter().map(|u| umi_map[u].count).sum();
                let top_umi_str = std::str::from_utf8(top_umi).unwrap_or("");

                for umi in group {
                    let slot = umi_map.remove(umi).expect("UMI must exist in umi_map");

                    for record in slot.records {
                        if let Some(w) = tsv_writer.as_mut() {
                            let read_name = std::str::from_utf8(record.qname()).unwrap_or("");
                            let contig_name =
                                std::str::from_utf8(header_view.tid2name(record.tid() as u32))
                                    .unwrap_or("");
                            let umi_str = std::str::from_utf8(umi).unwrap_or("");
                            let read_pos = five_prime_position(&record);
                            // umi_tools prints the contig, not the mapped gene, for --per-contig.
                            let gene_label = match assigner {
                                Some(assigner) if assigner.per_contig() => contig_name.to_string(),
                                Some(assigner) => assigner
                                    .gene(&record)
                                    .map(|gene| String::from_utf8_lossy(&gene).into_owned())
                                    .unwrap_or_default(),
                                None => "NA".to_string(),
                            };

                            writeln!(
                                w,
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                read_name,
                                contig_name,
                                read_pos,
                                gene_label,
                                umi_str,
                                slot.count,
                                top_umi_str,
                                group_count,
                                *unique_id,
                            )
                            .map_err(|e| GroupError::TsvWrite(e.to_string()))?;
                        }

                        let mut tagged = record;
                        if config.output_bam {
                            alignment_io::set_aux(
                                &mut tagged,
                                b"UG",
                                rust_htslib::bam::record::Aux::U32(*unique_id),
                            )
                            .map_err(|e| GroupError::BamWrite(e.to_string()))?;
                            alignment_io::set_aux(
                                &mut tagged,
                                &config.umi_group_tag,
                                rust_htslib::bam::record::Aux::String(top_umi_str),
                            )
                            .map_err(|e| GroupError::BamWrite(e.to_string()))?;
                        }
                        output
                            .push(tagged)
                            .map_err(|e| GroupError::BamWrite(e.to_string()))?;
                    }
                }

                *unique_id += 1;
            }
        }
    }

    Ok(())
}

/// # Errors
///
/// Returns `GroupError` on BAM I/O failures or unknown chromosome filter.
#[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
pub fn run_group(config: &GroupConfig, input_path: &str) -> Result<GroupStats, GroupError> {
    let mut source = alignment_io::RecordSource::whole(input_path, config.reference.as_deref())
        .map_err(|e| GroupError::BamOpen(e.to_string()))?;
    let assigner = GeneAssigner::new(&config.gene, source.header())?;
    if let Some(map) = assigner.as_ref().and_then(GeneAssigner::transcript_map) {
        source = alignment_io::RecordSource::by_contig(
            input_path,
            config.reference.as_deref(),
            map.genes.clone(),
        )
        .map_err(|e| GroupError::BamOpen(e.to_string()))?;
    }
    let mut flusher = assigner.as_ref().map(GeneAssigner::flusher);
    let header_view = source.header().clone();

    let writer = if config.output_bam {
        let header = if config.no_sort_output {
            bam::Header::from_template(&header_view)
        } else {
            alignment_io::coordinate_sorted_header(&header_view)
        };
        Some(
            alignment_io::open_writer(
                &header,
                AlignmentOutput {
                    path: config.output_path.as_deref(),
                    format: config.output_format,
                    reference: config.reference.as_deref(),
                },
            )
            .map_err(|e| GroupError::BamWrite(e.to_string()))?,
        )
    } else {
        None
    };

    // Optional chromosome filter
    let chrom_filter: Option<i32> = config
        .chrom
        .as_ref()
        .map(|c| {
            let tid = source
                .header()
                .tid(c.as_bytes())
                .ok_or_else(|| GroupError::UnknownChrom(c.clone()))?;
            #[allow(clippy::cast_possible_wrap)]
            Ok::<i32, GroupError>(tid as i32)
        })
        .transpose()?;

    // Open TSV writer
    let mut tsv_writer: Option<BufWriter<File>> = config
        .group_out
        .as_ref()
        .map(|path| {
            let file =
                File::create(path).map_err(|e| GroupError::TsvWrite(e.to_string()))?;
            let mut w = BufWriter::new(file);
            writeln!(
                w,
                "read_id\tcontig\tposition\tgene\tumi\tumi_count\tfinal_umi\tfinal_umi_count\tunique_id"
            )
            .map_err(|e| GroupError::TsvWrite(e.to_string()))?;
            Ok::<_, GroupError>(w)
        })
        .transpose()?;

    let mut buffer = GroupBuffer::<i64>::new();
    let mut gene_buffer = GroupBuffer::<Vec<u8>>::new();
    let mut stats = GroupStats {
        input_reads: 0,
        output_reads: 0,
    };

    #[allow(clippy::cast_possible_truncation)]
    let mut rng = PythonRandom::new(config.random_seed as u32);

    let mut output = RecordOutput::new(
        writer,
        (!config.no_sort_output).then(|| alignment_io::coordinate_sorted_header(&header_view)),
    );
    let mut unique_id: u32 = 0;

    let mut last_start: i64 = 0;
    let mut last_chrom: i32 = -1;

    while let Some(record) = source
        .read_next()
        .map_err(|e| GroupError::BamRead(e.to_string()))?
    {
        let tid = record.tid();
        if chrom_filter.is_some_and(|filter_tid| tid != filter_tid) {
            continue;
        }

        let triage = config.pairing.triage(&record, true);
        for _ in 0..triage.copies {
            output
                .push(record.clone())
                .map_err(|e| GroupError::BamWrite(e.to_string()))?;
        }
        if triage.is_read2 {
            continue;
        }
        stats.input_reads += 1;
        if !triage.grouped {
            continue;
        }

        // Subset check consumes one RNG call per grouped read (before buffer.add)
        if config.subset.is_some_and(|s| rng.random() >= f64::from(s)) {
            continue;
        }

        if record.mapq() < config.mapping_quality {
            continue;
        }

        let Some(Barcode { umi, cell }) =
            config.barcode.for_grouping(&record, config.ignore_umi)?
        else {
            continue;
        };

        if let (Some(assigner), Some(flusher)) = (&assigner, flusher.as_mut()) {
            let Some(gene) = assigner.gene(&record) else {
                continue;
            };
            let done = match flusher.before_read(tid) {
                Flush::None => None,
                Flush::All => Some(gene_buffer.drain_all()),
                Flush::Gene(done) => Some(gene_buffer.drain_key(&done)),
            };
            if let Some(drained) = done {
                process_drained(
                    drained,
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    Some(assigner),
                    &mut output,
                )?;
            }
            flusher.after_read(tid, &gene);

            let key: GroupKey = GroupKey::for_gene(cell);
            gene_buffer.add(record, gene, key, umi);
        } else {
            // Standard coordinate mode
            let position = get_read_position(&record, config.position.soft_clip_threshold);
            let start = position.start;

            if tid != last_chrom {
                process_drained(
                    buffer.drain_all(),
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    assigner.as_ref(),
                    &mut output,
                )?;
            } else if !config.buffer_whole_contig && start > last_start + 1000 {
                let threshold = start - 1000;
                process_drained(
                    buffer.drain_up_to(threshold),
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    assigner.as_ref(),
                    &mut output,
                )?;
            }

            last_start = start;
            last_chrom = tid;

            let key = GroupKey::for_position(
                &record,
                &position,
                &config.position,
                config.pairing.paired && !config.ignore_tlen,
                cell,
            );

            buffer.add(record, position.pos, key, umi);
        }
    }

    process_drained(
        buffer.drain_all(),
        config,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
        assigner.as_ref(),
        &mut output,
    )?;
    process_drained(
        gene_buffer.drain_all(),
        config,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
        assigner.as_ref(),
        &mut output,
    )?;

    // Flush TSV
    if let Some(w) = tsv_writer.as_mut() {
        w.flush().map_err(|e| GroupError::TsvWrite(e.to_string()))?;
    }

    stats.output_reads = output
        .finish()
        .map_err(|e| GroupError::BamWrite(e.to_string()))?;

    Ok(stats)
}

#[derive(Debug, thiserror::Error)]
pub enum GroupError {
    #[error("failed to open BAM: {0}")]
    BamOpen(String),
    #[error("failed to read BAM record: {0}")]
    BamRead(String),
    #[error("failed to write BAM/SAM: {0}")]
    BamWrite(String),
    #[error("failed to write TSV: {0}")]
    TsvWrite(String),
    #[error("unknown chromosome: {0}")]
    UnknownChrom(String),
    #[error("invalid regex: {0}")]
    InvalidRegex(String),
    #[error(transparent)]
    Barcode(#[from] BarcodeError),
    #[error(transparent)]
    Gene(#[from] GeneError),
}
