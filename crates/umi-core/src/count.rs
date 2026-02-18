use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{self, BufRead, Write as IoWrite};

use rust_htslib::bam::{self, Read as BamRead, record::Aux};
use thiserror::Error;

use crate::dedup::{DedupMethod, count_umis, extract_umi_umis};

#[derive(Error, Debug)]
pub enum CountError {
    #[error("BAM open error: {0}")]
    BamOpen(String),
    #[error("BAM read error: {0}")]
    BamRead(String),
    #[error("invalid regex: {0}")]
    InvalidRegex(String),
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
}

pub struct CountConfig {
    pub method: DedupMethod,
    pub gene_tag: String,
    pub skip_tags_regex: Option<String>,
    pub per_cell: bool,
    pub wide_format: bool,
    pub edit_distance_threshold: u32,
}

pub struct CountStats {
    pub input_reads: u64,
    pub counted_reads: u64,
}

pub struct CountTabConfig {
    pub method: DedupMethod,
    pub per_cell: bool,
    pub separator: u8,
    pub edit_distance_threshold: u32,
}

/// UMI count map: `umi -> (count, insertion_order)`.
type UmiCountMap = HashMap<Vec<u8>, (u32, u32)>;

#[allow(clippy::missing_errors_doc)]
pub fn run_count(
    config: &CountConfig,
    bam_path: &str,
    output: &mut dyn IoWrite,
) -> Result<CountStats, CountError> {
    let mut reader =
        bam::Reader::from_path(bam_path).map_err(|e| CountError::BamOpen(e.to_string()))?;

    let skip_regex = config
        .skip_tags_regex
        .as_ref()
        .map(|s| regex::Regex::new(s).map_err(|e| CountError::InvalidRegex(e.to_string())))
        .transpose()?;

    // BTreeMap for genes so output is sorted
    let mut data: BTreeMap<String, CellUmiMap> = BTreeMap::new();
    let mut stats = CountStats {
        input_reads: 0,
        counted_reads: 0,
    };

    for result in reader.records() {
        let record = result.map_err(|e| CountError::BamRead(e.to_string()))?;

        if record.is_unmapped() {
            continue;
        }
        if record.is_paired() && record.is_last_in_template() {
            continue;
        }

        stats.input_reads += 1;

        let gene = match record.aux(config.gene_tag.as_bytes()) {
            Ok(Aux::String(s)) => s.to_string(),
            _ => continue,
        };

        if skip_regex.as_ref().is_some_and(|re| re.is_match(&gene)) {
            continue;
        }

        let (umi, cell) = extract_umi_umis(record.qname());

        let cell_key = if config.per_cell {
            cell.map(|c| String::from_utf8_lossy(&c).into_owned())
        } else {
            None
        };

        stats.counted_reads += 1;

        let cell_map = data.entry(gene).or_default();
        cell_map.add(cell_key, umi);
    }

    if config.per_cell && config.wide_format {
        write_wide_format(&data, config, output)?;
    } else if config.per_cell {
        write_long_format(&data, config, output)?;
    } else {
        write_gene_counts(&data, config, output)?;
    }

    Ok(stats)
}

#[allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]
pub fn run_count_tab(
    config: &CountTabConfig,
    input: &mut dyn BufRead,
    output: &mut dyn IoWrite,
) -> Result<CountStats, CountError> {
    let mut stats = CountStats {
        input_reads: 0,
        counted_reads: 0,
    };

    if config.per_cell {
        writeln!(output, "cell\tgene\tcount")?;
    } else {
        writeln!(output, "gene\tcount")?;
    }

    let mut current_gene: Option<String> = None;
    let mut cell_umis = CellUmiMap::default();

    let mut line_buf = String::new();
    loop {
        line_buf.clear();
        let n = input.read_line(&mut line_buf)?;
        if n == 0 {
            break;
        }
        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }

        let mut cols = line.splitn(2, '\t');
        let Some(read_name) = cols.next() else {
            continue;
        };
        let Some(gene) = cols.next() else {
            continue;
        };
        let gene = gene.to_string();

        stats.input_reads += 1;

        // When gene changes, flush previous gene
        if current_gene.as_ref().is_some_and(|g| *g != gene) {
            flush_count_tab_gene(
                current_gene.as_deref().expect("checked above"),
                &cell_umis,
                config,
                output,
            )?;
            cell_umis = CellUmiMap::default();
        }
        current_gene = Some(gene);

        let sep = config.separator;
        let parts: Vec<&str> = read_name.split(|c: char| c as u8 == sep).collect();
        let umi = parts
            .last()
            .map_or_else(Vec::new, |s| s.as_bytes().to_vec());

        let cell_key = if config.per_cell && parts.len() >= 2 {
            Some(parts[parts.len() - 2].to_string())
        } else {
            None
        };

        stats.counted_reads += 1;
        cell_umis.add(cell_key, umi);
    }

    if let Some(ref gene) = current_gene {
        flush_count_tab_gene(gene, &cell_umis, config, output)?;
    }

    Ok(stats)
}

#[derive(Default)]
struct CellUmiMap {
    cells: Vec<(Option<String>, UmiCountMap)>,
    cell_index: HashMap<Option<String>, usize>,
    next_order: u32,
}

impl CellUmiMap {
    fn add(&mut self, cell: Option<String>, umi: Vec<u8>) {
        let idx = if let Some(&i) = self.cell_index.get(&cell) {
            i
        } else {
            let i = self.cells.len();
            self.cell_index.insert(cell.clone(), i);
            self.cells.push((cell, HashMap::new()));
            i
        };
        let entry = self.cells[idx].1.entry(umi).or_insert_with(|| {
            let order = self.next_order;
            self.next_order += 1;
            (0, order)
        });
        entry.0 += 1;
    }

    fn dedup_count(
        &self,
        method: DedupMethod,
        edit_threshold: u32,
    ) -> Vec<(&Option<String>, usize)> {
        self.cells
            .iter()
            .map(|(cell, umi_map)| {
                let counts: HashMap<Vec<u8>, u32> =
                    umi_map.iter().map(|(k, &(c, _))| (k.clone(), c)).collect();
                let orders: HashMap<Vec<u8>, u32> =
                    umi_map.iter().map(|(k, &(_, o))| (k.clone(), o)).collect();
                let n = count_umis(method, &counts, &orders, edit_threshold);
                (cell, n)
            })
            .collect()
    }
}

fn write_gene_counts(
    data: &BTreeMap<String, CellUmiMap>,
    config: &CountConfig,
    output: &mut dyn IoWrite,
) -> Result<(), CountError> {
    writeln!(output, "gene\tcount")?;
    for (gene, cell_map) in data {
        let results = cell_map.dedup_count(config.method, config.edit_distance_threshold);
        let total: usize = results.iter().map(|(_, n)| n).sum();
        writeln!(output, "{gene}\t{total}")?;
    }
    Ok(())
}

fn write_long_format(
    data: &BTreeMap<String, CellUmiMap>,
    config: &CountConfig,
    output: &mut dyn IoWrite,
) -> Result<(), CountError> {
    writeln!(output, "gene\tcell\tcount")?;
    for (gene, cell_map) in data {
        let results = cell_map.dedup_count(config.method, config.edit_distance_threshold);
        let mut sorted: Vec<_> = results
            .into_iter()
            .filter_map(|(cell, n)| cell.as_ref().map(|c| (c.clone(), n)))
            .collect();
        sorted.sort_by(|a, b| a.0.cmp(&b.0));
        for (cell, count) in sorted {
            writeln!(output, "{gene}\t{cell}\t{count}")?;
        }
    }
    Ok(())
}

fn write_wide_format(
    data: &BTreeMap<String, CellUmiMap>,
    config: &CountConfig,
    output: &mut dyn IoWrite,
) -> Result<(), CountError> {
    let mut all_cells: BTreeSet<String> = BTreeSet::new();
    for cell_map in data.values() {
        for (cell, _) in &cell_map.cells {
            if let Some(c) = cell {
                all_cells.insert(c.clone());
            }
        }
    }
    let cell_list: Vec<&String> = all_cells.iter().collect();

    write!(output, "gene")?;
    for cell in &cell_list {
        write!(output, "\t{cell}")?;
    }
    writeln!(output)?;

    for (gene, cell_map) in data {
        let results = cell_map.dedup_count(config.method, config.edit_distance_threshold);
        let cell_counts: HashMap<&str, usize> = results
            .into_iter()
            .filter_map(|(cell, n)| cell.as_ref().map(|c| (c.as_str(), n)))
            .collect();

        write!(output, "{gene}")?;
        for cell in &cell_list {
            let count = cell_counts.get(cell.as_str()).copied().unwrap_or(0);
            write!(output, "\t{count}")?;
        }
        writeln!(output)?;
    }
    Ok(())
}

fn flush_count_tab_gene(
    gene: &str,
    cell_umis: &CellUmiMap,
    config: &CountTabConfig,
    output: &mut dyn IoWrite,
) -> Result<(), CountError> {
    let results = cell_umis.dedup_count(config.method, config.edit_distance_threshold);

    if config.per_cell {
        for (cell, count) in results {
            let cell_str = cell.as_deref().unwrap_or("");
            writeln!(output, "{cell_str}\t{gene}\t{count}")?;
        }
    } else {
        let total: usize = results.iter().map(|(_, n)| n).sum();
        writeln!(output, "{gene}\t{total}")?;
    }
    Ok(())
}
