use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

use rust_htslib::bam::{self, Read as BamRead, Record};

use crate::dedup::{
    DedupMethod, GroupKey, PythonRandom, TieBreakRng, build_adjacency_list,
    build_directional_adjacency_list, connected_components, extract_umi_from_name,
    get_read_position, median, min_set_cover,
};

#[allow(clippy::struct_excessive_bools)]
pub struct GroupConfig {
    pub method: DedupMethod,
    pub ignore_umi: bool,
    pub umi_separator: u8,
    pub random_seed: u64,
    pub out_sam: bool,
    pub output_bam: bool,
    pub no_sort_output: bool,
    pub chrom: Option<String>,
    pub group_out: Option<String>,
    pub edit_distance_threshold: u32,
    pub subset: Option<f32>,
    pub output_unmapped: bool,
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

struct GroupBuffer {
    groups: BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, GroupSlot>>>,
    insertion_counters: BTreeMap<i64, BTreeMap<GroupKey, u32>>,
}

impl GroupBuffer {
    const fn new() -> Self {
        Self {
            groups: BTreeMap::new(),
            insertion_counters: BTreeMap::new(),
        }
    }

    fn add(&mut self, record: Record, pos: i64, key: GroupKey, umi: Vec<u8>) {
        let umi_map = self.groups.entry(pos).or_default().entry(key).or_default();

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

    fn drain_up_to(
        &mut self,
        threshold: i64,
    ) -> BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, GroupSlot>>> {
        let rest = self.groups.split_off(&(threshold + 1));
        let drained = std::mem::replace(&mut self.groups, rest);
        let rest_counters = self.insertion_counters.split_off(&(threshold + 1));
        let _ = std::mem::replace(&mut self.insertion_counters, rest_counters);
        drained
    }

    fn drain_all(&mut self) -> BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, GroupSlot>>> {
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
                    comp
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
                    groups.push(component);
                } else {
                    let lead_umis = min_set_cover(&component, &adj_list, &counts);
                    let mut observed: HashSet<Vec<u8>> = lead_umis.iter().cloned().collect();
                    for lead in &lead_umis {
                        let connected: HashSet<&Vec<u8>> = adj_list
                            .get(lead)
                            .map_or_else(HashSet::new, |ns| ns.iter().collect());
                        let mut group = vec![lead.clone()];
                        for node in connected {
                            if !observed.contains(node) {
                                group.push(node.clone());
                                observed.insert(node.clone());
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
            let mut observed: HashSet<Vec<u8>> = HashSet::new();
            let mut groups = Vec::new();
            for mut comp in components {
                comp.sort_by(|a, b| lex_sort(a, b));
                if comp.len() == 1 {
                    observed.insert(comp[0].clone());
                    groups.push(comp);
                } else {
                    let filtered: Vec<Vec<u8>> = comp
                        .into_iter()
                        .filter(|node| {
                            if observed.contains(node) {
                                false
                            } else {
                                observed.insert(node.clone());
                                true
                            }
                        })
                        .collect();
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
fn process_drained(
    drained: BTreeMap<i64, BTreeMap<GroupKey, HashMap<Vec<u8>, GroupSlot>>>,
    method: DedupMethod,
    edit_threshold: u32,
    unique_id: &mut u32,
    tsv_writer: &mut Option<BufWriter<File>>,
    header_view: &bam::HeaderView,
) -> Vec<Record> {
    let mut output_records = Vec::new();

    for (pos, key_map) in drained {
        for (_, mut umi_map) in key_map {
            let groups = assign_groups(method, &umi_map, edit_threshold);

            for group in &groups {
                let top_umi = &group[0];
                let group_count: u32 = group.iter().map(|u| umi_map[u].count).sum();
                let top_umi_str = std::str::from_utf8(top_umi).unwrap_or("");

                for umi in group {
                    let slot = umi_map.remove(umi).expect("UMI must exist in umi_map");

                    for record in slot.records {
                        if let Some(w) = tsv_writer.as_mut() {
                            let read_name = std::str::from_utf8(record.qname()).unwrap_or("");
                            let contig =
                                std::str::from_utf8(header_view.tid2name(record.tid() as u32))
                                    .unwrap_or("");
                            let umi_str = std::str::from_utf8(umi).unwrap_or("");

                            let _ = writeln!(
                                w,
                                "{}\t{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}",
                                read_name,
                                contig,
                                pos,
                                umi_str,
                                slot.count,
                                top_umi_str,
                                group_count,
                                *unique_id,
                            );
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
                            .push_aux(b"BX", rust_htslib::bam::record::Aux::String(top_umi_str))
                            .ok();

                        output_records.push(tagged);
                    }
                }

                *unique_id += 1;
            }
        }
    }

    output_records
}

/// # Errors
///
/// Returns `GroupError` on BAM I/O failures or unknown chromosome filter.
#[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
pub fn run_group(config: &GroupConfig, input_path: &str) -> Result<GroupStats, GroupError> {
    let mut reader =
        bam::Reader::from_path(input_path).map_err(|e| GroupError::BamOpen(e.to_string()))?;
    let header = bam::Header::from_template(reader.header());
    let header_view = reader.header().clone();

    let format = if config.out_sam {
        bam::Format::Sam
    } else {
        bam::Format::Bam
    };

    let mut writer = bam::Writer::from_stdout(&header, format)
        .map_err(|e| GroupError::BamWrite(e.to_string()))?;

    // Optional chromosome filter
    let chrom_filter: Option<i32> = config
        .chrom
        .as_ref()
        .map(|c| {
            let tid = reader
                .header()
                .tid(c.as_bytes())
                .ok_or_else(|| GroupError::UnknownChrom(c.clone()))?;
            #[allow(clippy::cast_possible_wrap)]
            Ok(tid as i32)
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
            Ok(w)
        })
        .transpose()?;

    let mut buffer = GroupBuffer::new();
    let mut stats = GroupStats {
        input_reads: 0,
        output_reads: 0,
    };

    #[allow(clippy::cast_possible_truncation)]
    let mut rng = PythonRandom::new(config.random_seed as u32);

    let mut output_records: Vec<Record> = Vec::new();
    let mut unmapped_records: Vec<Record> = Vec::new();
    let mut unique_id: u32 = 0;

    let mut last_start: i64 = 0;
    let mut last_chrom: i32 = -1;

    for result in reader.records() {
        let record = result.map_err(|e| GroupError::BamRead(e.to_string()))?;

        if record.is_unmapped() {
            if config.output_unmapped {
                unmapped_records.push(record);
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

        let (start, pos) = get_read_position(&record);

        if tid != last_chrom {
            output_records.extend(process_drained(
                buffer.drain_all(),
                config.method,
                config.edit_distance_threshold,
                &mut unique_id,
                &mut tsv_writer,
                &header_view,
            ));
        } else if start > last_start + 1000 {
            let threshold = start - 1000;
            output_records.extend(process_drained(
                buffer.drain_up_to(threshold),
                config.method,
                config.edit_distance_threshold,
                &mut unique_id,
                &mut tsv_writer,
                &header_view,
            ));
        }

        last_start = start;
        last_chrom = tid;

        let key: GroupKey = (record.is_reverse(), false, 0, 0);

        let umi = if config.ignore_umi {
            Vec::new()
        } else {
            extract_umi_from_name(&record, config.umi_separator)
        };

        buffer.add(record, pos, key, umi);
    }

    output_records.extend(process_drained(
        buffer.drain_all(),
        config.method,
        config.edit_distance_threshold,
        &mut unique_id,
        &mut tsv_writer,
        &header_view,
    ));

    // Flush TSV
    if let Some(w) = tsv_writer.as_mut() {
        w.flush().map_err(|e| GroupError::TsvWrite(e.to_string()))?;
    }

    // Sort by coordinate unless --no-sort-output
    if !config.no_sort_output {
        output_records.sort_by(|a, b| a.tid().cmp(&b.tid()).then_with(|| a.pos().cmp(&b.pos())));
    }

    // Append unmapped reads after sorted mapped reads (no tags added)
    output_records.extend(unmapped_records);

    stats.output_reads = output_records.len() as u64;

    if config.output_bam {
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
}
