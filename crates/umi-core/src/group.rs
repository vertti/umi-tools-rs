use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

use rust_htslib::bam::{self, Record};

use crate::alignment_io::{self, AlignmentFormat, AlignmentOutput};
use crate::barcode::{Barcode, BarcodeError, BarcodeExtractor};
use crate::dedup::{
    DedupMethod, GroupKey, PositionOptions, PythonRandom, TieBreakRng, build_adjacency_list,
    build_directional_adjacency_list, connected_components, five_prime_position, get_read_position,
    median, min_set_cover,
};
use crate::gene::{Flush, GeneAssigner, GeneError, GeneOptions};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ChimericPairs {
    Discard,
    Output,
    Use,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum UnmappedHandling {
    Discard,
    Output,
    Use,
}

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
    pub paired: bool,
    pub chimeric_pairs: ChimericPairs,
    pub unmapped_handling: UnmappedHandling,
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

/// Assign UMIs to groups. Returns groups where each group is a list of UMIs
/// sorted by count descending, lex ascending. First UMI is the representative.
#[allow(clippy::too_many_lines)]
fn assign_groups(
    method: DedupMethod,
    umi_map: &HashMap<Vec<u8>, GroupSlot>,
    edit_threshold: u32,
) -> Vec<Vec<Vec<u8>>> {
    let counts: HashMap<&[u8], u32> = umi_map
        .iter()
        .map(|(k, v)| (k.as_slice(), v.count))
        .collect();
    let orders: HashMap<&[u8], u32> = umi_map
        .iter()
        .map(|(k, v)| (k.as_slice(), v.insertion_order))
        .collect();

    let lex_sort = |a: &[u8], b: &[u8]| -> std::cmp::Ordering {
        counts[b].cmp(&counts[a]).then_with(|| a.cmp(b))
    };

    match method {
        DedupMethod::Unique => {
            let mut umis: Vec<Vec<u8>> = umi_map.keys().cloned().collect();
            umis.sort_by(|a, b| orders[a.as_slice()].cmp(&orders[b.as_slice()]));
            umis.into_iter().map(|u| vec![u]).collect()
        }

        DedupMethod::Percentile => {
            if counts.len() <= 1 {
                return umi_map.keys().cloned().map(|u| vec![u]).collect();
            }
            let all_counts: Vec<u32> = counts.values().copied().collect();
            let threshold = median(&all_counts) / 100.0;
            let mut umis: Vec<Vec<u8>> = umi_map
                .iter()
                .filter(|(_, slot)| f64::from(slot.count) > threshold)
                .map(|(umi, _)| umi.clone())
                .collect();
            umis.sort_by(|a, b| orders[a.as_slice()].cmp(&orders[b.as_slice()]));
            umis.into_iter().map(|u| vec![u]).collect()
        }

        DedupMethod::Cluster => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            components
                .into_iter()
                .map(|mut comp| {
                    comp.sort_by(|a, b| lex_sort(a, b));
                    comp.into_iter().map(<[u8]>::to_vec).collect()
                })
                .collect()
        }

        DedupMethod::Adjacency => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_adjacency_list(&umis, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            // Adjacency splits components via min_set_cover, grouping
            // connected nodes around each lead UMI.
            let mut groups = Vec::new();
            for component in components {
                if component.len() == 1 {
                    groups.push(component.into_iter().map(<[u8]>::to_vec).collect());
                } else {
                    let lead_umis = min_set_cover(&component, &adj_list, &counts);
                    let mut observed: HashSet<&[u8]> = lead_umis.iter().copied().collect();
                    for &lead in &lead_umis {
                        let connected: HashSet<&[u8]> = adj_list
                            .get(lead)
                            .map_or_else(HashSet::new, |ns| ns.iter().copied().collect());
                        let mut group = vec![lead.to_vec()];
                        for node in connected {
                            if observed.insert(node) {
                                group.push(node.to_vec());
                            }
                        }
                        groups.push(group);
                    }
                }
            }
            groups
        }

        DedupMethod::Directional => {
            let umis: Vec<&[u8]> = umi_map.keys().map(Vec::as_slice).collect();
            let adj_list = build_directional_adjacency_list(&umis, &counts, edit_threshold);
            let components = connected_components(&umis, &counts, &orders, &adj_list);
            // Directed BFS can produce overlapping components. Filter already-
            // observed UMIs so each UMI is assigned to exactly one group,
            // matching Python's _group_directional logic.
            let mut observed: HashSet<&[u8]> = HashSet::new();
            let mut groups = Vec::new();
            for mut comp in components {
                comp.sort_by(|a, b| lex_sort(a, b));
                if comp.len() == 1 {
                    observed.insert(comp[0]);
                    groups.push(comp.into_iter().map(<[u8]>::to_vec).collect());
                } else {
                    let mut filtered: Vec<Vec<u8>> = Vec::new();
                    for node in comp {
                        if observed.insert(node) {
                            filtered.push(node.to_vec());
                        }
                    }
                    if !filtered.is_empty() {
                        groups.push(filtered);
                    }
                }
            }
            groups
        }
    }
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
) -> Result<Vec<Record>, GroupError> {
    let mut output_records = Vec::new();

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
                        tagged
                            .push_aux(
                                b"UG",
                                #[allow(clippy::cast_possible_wrap)]
                                rust_htslib::bam::record::Aux::I32(*unique_id as i32),
                            )
                            .ok();
                        tagged
                            .push_aux(
                                &config.umi_group_tag,
                                rust_htslib::bam::record::Aux::String(top_umi_str),
                            )
                            .ok();

                        output_records.push(tagged);
                    }
                }

                *unique_id += 1;
            }
        }
    }

    Ok(output_records)
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

    let mut writer = if config.output_bam {
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

    let output_unmapped = config.unmapped_handling == UnmappedHandling::Output
        || config.unmapped_handling == UnmappedHandling::Use;

    let mut buffer = GroupBuffer::<i64>::new();
    let mut gene_buffer = GroupBuffer::<Vec<u8>>::new();
    let mut stats = GroupStats {
        input_reads: 0,
        output_reads: 0,
    };

    #[allow(clippy::cast_possible_truncation)]
    let mut rng = PythonRandom::new(config.random_seed as u32);

    let mut output_records: Vec<Record> = Vec::new();
    let mut unique_id: u32 = 0;

    let mut last_start: i64 = 0;
    let mut last_chrom: i32 = -1;

    while let Some(record) = source
        .read_next()
        .map_err(|e| GroupError::BamRead(e.to_string()))?
    {
        // R2 reads are passthrough (no grouping, no tags).
        if record.is_last_in_template() {
            if record.is_unmapped() {
                if output_unmapped {
                    output_records.push(record);
                }
            } else {
                output_records.push(record);
            }
            continue;
        }

        // Handle unmapped reads (R1 in paired mode, or any read in single-end)
        if record.is_unmapped() {
            if output_unmapped {
                output_records.push(record);
            }
            continue;
        }

        let tid = record.tid();

        if chrom_filter.is_some_and(|filter_tid| tid != filter_tid) {
            continue;
        }

        stats.input_reads += 1;

        // Subset check consumes one RNG call per mapped read (before buffer.add)
        if config.subset.is_some_and(|s| rng.random() >= f64::from(s)) {
            continue;
        }

        // Paired-mode filtering for R1 reads
        if config.paired {
            let is_chimeric =
                !record.is_mate_unmapped() && record.tid() != record.mtid() && record.mtid() >= 0;

            if is_chimeric {
                match config.chimeric_pairs {
                    ChimericPairs::Discard => continue,
                    ChimericPairs::Output => {
                        output_records.push(record);
                        continue;
                    }
                    ChimericPairs::Use => {} // fall through to grouping with TLEN=0
                }
            }

            if record.is_mate_unmapped() {
                match config.unmapped_handling {
                    UnmappedHandling::Discard => continue,
                    UnmappedHandling::Output => {
                        output_records.push(record);
                        continue;
                    }
                    UnmappedHandling::Use => {} // fall through to grouping with TLEN=0
                }
            }
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
                output_records.extend(process_drained(
                    drained,
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    Some(assigner),
                )?);
            }
            flusher.after_read(tid, &gene);

            let key: GroupKey = (false, 0, 0, 0, cell);
            gene_buffer.add(record, gene, key, umi);
        } else {
            // Standard coordinate mode
            let position = get_read_position(&record, config.position.soft_clip_threshold);
            let start = position.start;

            if tid != last_chrom {
                output_records.extend(process_drained(
                    buffer.drain_all(),
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    assigner.as_ref(),
                )?);
            } else if !config.buffer_whole_contig && start > last_start + 1000 {
                let threshold = start - 1000;
                output_records.extend(process_drained(
                    buffer.drain_up_to(threshold),
                    config,
                    &mut unique_id,
                    &mut tsv_writer,
                    &header_view,
                    assigner.as_ref(),
                )?);
            }

            last_start = start;
            last_chrom = tid;

            // For paired non-chimeric reads, include signed TLEN in the group key.
            // Python sorts GroupKeys as tuples: (is_reverse, is_spliced, tlen, r_length).
            // We place signed tlen in position 2 (i64) to match Python's sorted() ordering.
            let tlen =
                if config.paired && !record.is_mate_unmapped() && record.tid() == record.mtid() {
                    record.insert_size()
                } else {
                    0
                };
            let (splice, length) = config.position.key_parts(&position, &record);
            let key: GroupKey = (record.is_reverse(), splice, tlen, length, cell);

            buffer.add(record, position.pos, key, umi);
        }
    }

    output_records.extend(process_drained(
        buffer.drain_all(),
        config,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
        assigner.as_ref(),
    )?);
    output_records.extend(process_drained(
        gene_buffer.drain_all(),
        config,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
        assigner.as_ref(),
    )?);

    // Flush TSV
    if let Some(w) = tsv_writer.as_mut() {
        w.flush().map_err(|e| GroupError::TsvWrite(e.to_string()))?;
    }

    // Sort by coordinate unless --no-sort-output.
    // Unmapped reads are placed after all mapped reads (matching Python).
    if !config.no_sort_output {
        let (mut mapped, unmapped): (Vec<_>, Vec<_>) =
            output_records.into_iter().partition(|r| !r.is_unmapped());
        mapped.sort_by(|a, b| a.tid().cmp(&b.tid()).then_with(|| a.pos().cmp(&b.pos())));
        mapped.extend(unmapped);
        output_records = mapped;
    }

    stats.output_reads = output_records.len() as u64;

    if let Some(writer) = writer.as_mut() {
        for r in &output_records {
            writer
                .write(r)
                .map_err(|e| GroupError::BamWrite(e.to_string()))?;
        }
    }

    drop(writer);

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
