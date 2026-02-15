# Roadmap

Implementation roadmap derived from the umi-tools compatibility test suite (`tests/tests.yaml` in the Python umi-tools repo). The test suite has **66 tests** total; we use it to track feature parity.

## Current state: 23 / 66 tests passing

| Phase | Tests | Status |
|-------|-------|--------|
| Baseline | `extract_single_string`, `extract_single`, `extract_3prime`, `extract_quality` | Done |
| Phase 1 | `extract_read2_only_string`, `extract_read2_only_regex` | Done |
| Phase 2 | `extract_scrb_seq`, `extract_scrb_seq_string` | Done |
| Phase 3 | `extract_scrb_seq_suffix`, `extract_scrb_seq_prefiltered` | Done |
| Phase 4 | `extract_indrop_fuzzy`, `extract_indrop_blacklist`, `extract_indrop_output_filtered` | Done |
| Phase 5 | `extract_either_read` | Done |
| Phase 6 | `whitelist_scrb_seq`, `whitelist_indrop`, `whitelist_indrop_set_cell`, `whitelist_indrop_3_errors`, `whitelist_indrop_density`, `whitelist_indrop_expect_cells_density`, `whitelist_indrop_filtered_out`, `whitelist_indrop_ed_above_threshold_discard`, `whitelist_indrop_ed_above_threshold_correct` | Done |

---

## Phase 7: BAM infrastructure + dedup (16 tests)

New subcommand: `dedup`. Requires BAM I/O (rust-htslib or noodles).

| Test | Key flags |
|------|-----------|
| `dedup_single_ignore` | `--ignore-umi` |
| `dedup_single_chrom` | `--chrom=chr19` |
| `dedup_single_unique` | `--method=unique` |
| `dedup_single_perc` | `--method=percentile` |
| `dedup_single_cluster` | `--method=cluster` |
| `dedup_single_adj` | `--method=adjacency` |
| `dedup_single_dir` | `--method=directional` |
| `dedup_single_stats` | `--method=cluster`, `--output-stats` |
| `dedup_single_dir_edit_dist` | `--method=directional`, `--edit-distance-threshold=2` |
| `dedup_single_subset` | `--method=directional`, `--subset=0.1` |
| `dedup_single_sep` | `--umi-separator=:` |
| `dedup_single_tag` | `--umi-tag=RX`, `--extract-umi-method=tag` |
| `dedup_single_tag_missing` | `--umi-tag=RX`, `--extract-umi-method=tag`, `--output-stats` |
| `dedup_single_gene_tag` | `--per-gene`, `--gene-tag=XF`, `--skip-tags-regex` |
| `dedup_paired_umi_whitelist` | `--paired`, `--filter-umi`, `--umi-whitelist`, `--umi-whitelist-paired` |
| `dedup_paired_ignore_tlen_tag` | `--paired`, `--ignore-tlen` |

Key features: 5 dedup methods (unique, percentile, cluster, adjacency, directional), BAM read/write, UMI extraction from read name or tags, per-gene mode, paired-end support, output stats.

## Phase 8: Group subcommand (15 tests)

New subcommand: `group`. Similar BAM infrastructure to dedup but outputs group assignments.

| Test | Key flags |
|------|-----------|
| `group_gene_tag` | `--per-gene`, `--gene-tag=XF`, `--group-out`, `--output-bam` |
| `group_unique` | `--method=unique`, `--group-out` |
| `group_cluster` | `--method=cluster`, `--group-out` |
| `group_adjacency` | `--method=adjacency`, `--group-out` |
| `group_directional` | `--method=directional`, `--group-out` |
| `group_directional_subset` | `--method=directional`, `--subset=0.1`, `--group-out` |
| `group_directional_unmapped` | `--method=directional`, `--output-unmapped`, `--group-out` |
| `group_unsorted` | `--method=directional`, `--no-sort-output` |
| `group_paired_discard_chimeric` | `--paired`, `--chimeric-pairs=discard` |
| `group_paired_output_chimeric` | `--paired`, `--chimeric-pairs=output` |
| `group_paired_use_chimeric` | `--paired`, `--chimeric-pairs=use` |
| `group_paired_discard_unmapped` | `--paired`, `--unmapped=discard` |
| `group_paired_output_unmapped` | `--unmapped=output` |
| `group_paired_use_unmapped` | `--paired`, `--unmapped=use` |
| `group_contig_no_gene_tag` | `--per-contig`, `--per-gene`, `--group-out` |

Key features: group TSV output, output-bam with BX/UG tags, chimeric/unmapped pair handling, per-contig mode.

## Phase 9: Count subcommands (5 tests)

New subcommands: `count`, `count_tab`.

| Test | Key flags |
|------|-----------|
| `count_single_gene_tag` | `count`, `--gene-tag=XF`, `--skip-tags-regex`, `--extract-umi-method=umis` |
| `count_single_cells_gene_tag` | `count`, `--per-cell` |
| `count_single_cells_wide_gene_tag` | `count`, `--per-cell`, `--wide-format-cell-counts` |
| `count_tab_single` | `count_tab` (TSV input) |
| `count_tab_single_per_cell` | `count_tab`, `--per-cell` |

Key features: gene-level UMI counting from BAM, per-cell counts, wide-format output, tab-delimited input mode.

---

## Not tracked: help tests (7 tests)

The test suite includes `--help` output tests for each subcommand. These will pass naturally as subcommands are implemented and are not tracked as separate milestones.
