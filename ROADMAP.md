# Roadmap

Implementation roadmap derived from the umi-tools compatibility test suite (`tests/tests.yaml` in the Python umi-tools repo). The test suite has **66 tests** total; we use it to track feature parity.

## Current state: 4 / 66 tests passing

| Test | Status |
|------|--------|
| `extract_single_string` | Done |
| `extract_single` | Done |
| `extract_3prime` | Done |
| `extract_quality` | Done |

All four are single-end extract with stdin/stdout, covering string method, regex method, 3' patterns, and quality filtering.

---

## Phase 1: Paired-end extract basics (2 tests)

New flags: `--read2-in`, `--read2-out`, `--bc-pattern2`

| Test | Method | Key flags |
|------|--------|-----------|
| `extract_read2_only_string` | string | `--read2-in`, `--bc-pattern2=NNNNNNNNNN`, `--read2-out` |
| `extract_read2_only_regex` | regex | `--read2-in`, `--bc-pattern2="^(?P<umi_1>.{10}).*"`, `--read2-out` |

These extract UMI from read2 only (pattern on read2, not read1). Stdout gets processed read1; `--read2-out` gets processed read2.

## Phase 2: Paired-end + whitelist (2 tests)

New flags: `--read2-stdout`, `--whitelist`, cell barcode support (`cell_N` capture groups)

| Test | Method | Key flags |
|------|--------|-----------|
| `extract_scrb_seq` | regex | `--read2-in`, `--read2-stdout`, `--whitelist`, `--bc-pattern="^(?P<cell_1>.{6})(?P<umi_1>.{10})"` |
| `extract_scrb_seq_string` | string | `--read2-in`, `--read2-stdout`, `--whitelist`, `--bc-pattern=CCCCCCNNNNNNNNNN` |

Pattern is on read1 (cell + UMI); `--read2-stdout` outputs processed read2. Whitelist filters cell barcodes.

## Phase 3: Paired-end edge cases (2 tests)

New flags: `--ignore-read-pair-suffixes`, `--reconcile-pairs`

| Test | Method | Key flags |
|------|--------|-----------|
| `extract_scrb_seq_suffix` | string | Phase 2 flags + `--ignore-read-pair-suffixes` |
| `extract_scrb_seq_prefiltered` | string | `--read2-in`, `--read2-stdout`, `--reconcile-pairs` (no whitelist) |

`--ignore-read-pair-suffixes` strips `/1` `/2` suffixes before matching read names. `--reconcile-pairs` handles pre-filtered inputs where read1/read2 sets may not be in sync.

## Phase 4: Fuzzy regex + cell error correction (3 tests)

New features: `{s<=N}` fuzzy regex syntax, `--error-correct-cell`, `--blacklist`, `--filtered-out`/`--filtered-out2`

| Test | Key flags |
|------|-----------|
| `extract_indrop_fuzzy` | fuzzy regex (`{s<=2}`), `--error-correct-cell`, `--whitelist` |
| `extract_indrop_blacklist` | exact regex, `--error-correct-cell`, `--whitelist`, `--blacklist` |
| `extract_indrop_output_filtered` | fuzzy regex, `--error-correct-cell`, `--whitelist`, `--filtered-out`, `--filtered-out2` |

The indrop pattern: `(?P<cell_1>.{8,12})(?P<discard_2>GAGTGATTGCTTGTGACGCCTT{s<=2})(?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*`

## Phase 5: Either-read mode (1 test)

New flags: `--either-read`, `--bc-pattern2` with fuzzy regex

| Test | Key flags |
|------|-----------|
| `extract_either_read` | `--either-read`, `--bc-pattern2` (same fuzzy pattern on both reads), `--read2-out` |

Tries pattern on both reads; uses whichever matches.

## Phase 6: Whitelist subcommand (9 tests)

New subcommand: `whitelist`. Identifies valid cell barcodes from FASTQ using knee-point detection.

| Test | Key flags |
|------|-----------|
| `whitelist_scrb_seq` | `--bc-pattern=CCCCCCNNNNNNNNNN`, `--plot-prefix` |
| `whitelist_indrop` | regex pattern, basic |
| `whitelist_indrop_filtered_out` | `--filtered-out` |
| `whitelist_indrop_density` | `--knee-method=density`, `--plot-prefix` |
| `whitelist_indrop_set_cell` | `--set-cell-number=1000` |
| `whitelist_indrop_expect_cells_density` | `--knee-method=density`, `--expect-cells=6000` |
| `whitelist_indrop_3_errors` | `--error-correct-threshold=3` |
| `whitelist_indrop_ed_above_threshold_discard` | `--ed-above-threshold=discard`, `--error-correct-threshold=3` |
| `whitelist_indrop_ed_above_threshold_correct` | `--ed-above-threshold=correct`, `--error-correct-threshold=3` |

Key features: knee-point methods (distance, density), cell number override, error correction thresholds, edit-distance above-threshold handling.

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
