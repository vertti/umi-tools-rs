use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::{self, BufRead, Write as IoWrite};

use thiserror::Error;

use crate::alignment_io::RecordSource;
use crate::barcode::{Barcode, BarcodeError, BarcodeExtractor};
use crate::dedup::{DedupMethod, PythonRandom, TieBreakRng, count_umis};
use crate::gene::{Flush, GeneAssigner, GeneError, GeneOptions};
use crate::pairing::{PairingError, PairingOptions};

#[derive(Error, Debug)]
pub enum CountError {
    #[error("BAM open error: {0}")]
    BamOpen(String),
    #[error("BAM read error: {0}")]
    BamRead(String),
    #[error("count needs --gene-tag or --per-contig")]
    NotPerGene,
    #[error("unknown chromosome: {0}")]
    UnknownChrom(String),
    #[error(transparent)]
    Pairing(#[from] PairingError),
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error(transparent)]
    Barcode(#[from] BarcodeError),
    #[error(transparent)]
    Gene(#[from] GeneError),
}

pub struct CountConfig {
    pub method: DedupMethod,
    pub gene: GeneOptions,
    pub barcode: BarcodeExtractor,
    pub ignore_umi: bool,
    pub wide_format: bool,
    pub edit_distance_threshold: u32,
    pub reference: Option<String>,
    pub mapping_quality: u8,
    pub pairing: PairingOptions,
    pub chrom: Option<String>,
    pub subset: Option<f32>,
    pub random_seed: u64,
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

/// Counts unique molecules per gene (and cell), one row per bundle in the
/// order `get_bundles` yields them: genes flush when the contig changes.
///
/// # Errors
///
/// Returns BAM, option and I/O errors.
pub fn run_count(
    config: &CountConfig,
    bam_path: &str,
    output: &mut dyn IoWrite,
) -> Result<CountStats, CountError> {
    config.pairing.validate(false)?;
    let mut source = RecordSource::whole(bam_path, config.reference.as_deref())
        .map_err(|e| CountError::BamOpen(e.to_string()))?;
    let assigner =
        GeneAssigner::new(&config.gene, source.header())?.ok_or(CountError::NotPerGene)?;
    let chrom_filter: Option<i32> = config
        .chrom
        .as_ref()
        .map(|c| {
            let tid = source
                .header()
                .tid(c.as_bytes())
                .ok_or_else(|| CountError::UnknownChrom(c.clone()))?;
            #[allow(clippy::cast_possible_wrap)]
            Ok::<i32, CountError>(tid as i32)
        })
        .transpose()?;
    #[allow(clippy::cast_possible_truncation)]
    let mut rng = PythonRandom::new(config.random_seed as u32);
    if let Some(map) = assigner.transcript_map() {
        source = RecordSource::by_contig(bam_path, config.reference.as_deref(), map.genes.clone())
            .map_err(|e| CountError::BamOpen(e.to_string()))?;
    }
    let mut flusher = assigner.flusher();

    let mut buffer: GeneBuffer = BTreeMap::new();
    let mut rows: Vec<CountRow> = Vec::new();
    let mut stats = CountStats {
        input_reads: 0,
        counted_reads: 0,
    };

    while let Some(record) = source
        .read_next()
        .map_err(|e| CountError::BamRead(e.to_string()))?
    {
        let tid = record.tid();
        if chrom_filter.is_some_and(|filter_tid| tid != filter_tid) {
            continue;
        }

        let triage = config.pairing.triage(&record, false);
        if triage.is_read2 {
            continue;
        }
        stats.input_reads += 1;
        if !triage.grouped {
            continue;
        }

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

        let Some(gene) = assigner.gene(&record) else {
            continue;
        };

        match flusher.before_read(tid) {
            Flush::None => {}
            Flush::All => flush_genes(std::mem::take(&mut buffer), config, &mut rows),
            Flush::Gene(done) => {
                if let Some(cells) = buffer.remove(&done) {
                    flush_genes(BTreeMap::from([(done, cells)]), config, &mut rows);
                }
            }
        }
        flusher.after_read(tid, &gene);

        stats.counted_reads += 1;
        buffer
            .entry(gene)
            .or_default()
            .entry(cell)
            .or_default()
            .add(umi);
    }
    flush_genes(buffer, config, &mut rows);

    if config.barcode.per_cell {
        // A later bundle for the same gene and cell overwrites the earlier one, as umi_tools' dict does.
        let mut table: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
        for row in rows {
            table
                .entry(row.gene)
                .or_default()
                .insert(row.cell, row.count);
        }
        if config.wide_format {
            write_wide_format(&table, output)?;
        } else {
            write_long_format(&table, output)?;
        }
    } else {
        writeln!(output, "gene\tcount")?;
        for row in rows {
            writeln!(output, "{}\t{}", row.gene, row.count)?;
        }
    }

    Ok(stats)
}

type GeneBuffer = BTreeMap<Vec<u8>, BTreeMap<Vec<u8>, UmiCounts>>;

struct CountRow {
    gene: String,
    cell: String,
    count: usize,
}

/// One row per (gene, cell) bundle: genes in name order, cells in name order.
fn flush_genes(buffer: GeneBuffer, config: &CountConfig, rows: &mut Vec<CountRow>) {
    for (gene, cells) in buffer {
        for (cell, umis) in cells {
            rows.push(CountRow {
                gene: String::from_utf8_lossy(&gene).into_owned(),
                cell: String::from_utf8_lossy(&cell).into_owned(),
                count: umis.dedup_count(config.method, config.edit_distance_threshold),
            });
        }
    }
}

#[derive(Default)]
struct UmiCounts {
    counts: UmiCountMap,
    next_order: u32,
}

impl UmiCounts {
    fn add(&mut self, umi: Vec<u8>) {
        let order = self.next_order;
        let entry = self.counts.entry(umi).or_insert((0, order));
        if entry.0 == 0 {
            self.next_order += 1;
        }
        entry.0 += 1;
    }

    fn dedup_count(&self, method: DedupMethod, edit_threshold: u32) -> usize {
        let counts: HashMap<Vec<u8>, u32> = self
            .counts
            .iter()
            .map(|(k, &(c, _))| (k.clone(), c))
            .collect();
        let orders: HashMap<Vec<u8>, u32> = self
            .counts
            .iter()
            .map(|(k, &(_, o))| (k.clone(), o))
            .collect();
        count_umis(method, &counts, &orders, edit_threshold)
    }
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

fn write_long_format(
    table: &BTreeMap<String, BTreeMap<String, usize>>,
    output: &mut dyn IoWrite,
) -> Result<(), CountError> {
    writeln!(output, "gene\tcell\tcount")?;
    for (gene, cells) in table {
        for (cell, count) in cells {
            writeln!(output, "{gene}\t{cell}\t{count}")?;
        }
    }
    Ok(())
}

fn write_wide_format(
    table: &BTreeMap<String, BTreeMap<String, usize>>,
    output: &mut dyn IoWrite,
) -> Result<(), CountError> {
    let all_cells: BTreeSet<&String> = table.values().flat_map(BTreeMap::keys).collect();

    write!(output, "gene")?;
    for cell in &all_cells {
        write!(output, "\t{cell}")?;
    }
    writeln!(output)?;

    for (gene, cells) in table {
        write!(output, "{gene}")?;
        for cell in &all_cells {
            let count = cells.get(*cell).copied().unwrap_or(0);
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
