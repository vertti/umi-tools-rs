# umi-tools-rs

[![CI](https://github.com/vertti/umi-tools-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/vertti/umi-tools-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Crates.io](https://img.shields.io/crates/v/umi-tools-rs)](https://crates.io/crates/umi-tools-rs)

A drop-in replacement for [UMI-tools](https://github.com/CGATOxford/UMI-tools), written in Rust. Same flags, same output — just faster.

## Performance

Measured for v2.0.0 with [hyperfine](https://github.com/sharkdp/hyperfine) against UMI-tools 1.1.6 on the bundled benchmark inputs:

| Command | Speedup |
|:--------|--------:|
| `extract` | **13.2–33.7x** |
| `whitelist` | **31.6x** |
| `dedup` | **34.6x** |
| `group` | **14.1x** |
| `count` | **80.1x** |
| `count_tab` | **76.2x** |

These are end-to-end command timings, including process startup. Speedups depend on input size and workload; the short count benchmarks have higher timing variability.

Run benchmarks yourself: `mise run bench`

## Installation

### From crates.io

```sh
cargo install umi-tools-rs
```

### Bioconda

Community-maintained [Bioconda packages](https://anaconda.org/bioconda/umi-tools-rs/files) are available for Linux x86_64 and Intel macOS:

```sh
conda create -n umi-tools-rs --override-channels \
  -c conda-forge -c bioconda --strict-channel-priority umi-tools-rs
conda activate umi-tools-rs
```

Bioconda versions may lag behind GitHub releases. For native Apple Silicon builds, use the prebuilt binaries below or install from crates.io.

### Prebuilt binaries

Download from [GitHub Releases](https://github.com/vertti/umi-tools-rs/releases/latest):

```sh
# Linux x86_64
curl -fsSL https://github.com/vertti/umi-tools-rs/releases/latest/download/umi-tools-rs-x86_64-unknown-linux-gnu.tar.gz \
  | tar xz -C /usr/local/bin

# macOS Apple Silicon
curl -fsSL https://github.com/vertti/umi-tools-rs/releases/latest/download/umi-tools-rs-aarch64-apple-darwin.tar.gz \
  | tar xz -C /usr/local/bin
```

### From source

```sh
cargo install --path .
```

## Usage

Replace `umi_tools` with `umi-tools-rs` in your existing commands:

```sh
# Extract UMIs
umi-tools-rs extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --stdin=R1.fastq.gz --read2-in=R2.fastq.gz \
  --stdout=R1_extracted.fastq.gz --read2-out=R2_extracted.fastq.gz

# Whitelist cell barcodes
umi-tools-rs whitelist --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --stdin=R1.fastq.gz --stdout=whitelist.tsv

# Deduplicate BAM reads
umi-tools-rs dedup --method=directional --stdin=aligned.bam --stdout=deduped.bam

# Count unique molecules per gene
umi-tools-rs count --gene-tag=XF --per-cell --stdin=aligned.bam > counts.tsv
```

### CRAM

`dedup`, `group` and `count` read CRAM, and `dedup` and `group` write it. The output format follows the `--stdout` file extension, or set it with `--out-format=cram`. Pass `--reference-filename` when the `UR` field of the input header does not point at a local FASTA:

```sh
umi-tools-rs dedup --stdin=aligned.cram --reference-filename=genome.fa --stdout=deduped.cram
```

## Testing

`mise run compat` runs the upstream UMI-tools test suite against the Rust binary (it expects a checkout of UMI-tools next to this repository, or `UMI_TOOLS_TESTS_DIR`). Flags the upstream suite does not exercise are covered by `tests/golden/cases.yaml`, whose references come from Python umi_tools; `mise run golden` regenerates them.

Synthetic extraction tests also cover combinations of patterns, output destinations, filtering, and pair reconciliation, with expected outputs checked against Python umi_tools.

`mise run duplicates` runs pinned [jscpd](https://github.com/kucherenko/jscpd) and fails on new duplicated Rust blocks of at least 10 lines and 100 tokens, including copies with renamed identifiers. CI runs this alongside the Rust checks. Inline test modules are excluded; existing production duplicates are recorded in `.jscpd-baseline.json` so they can be removed incrementally. This catches copied blocks, not every repeated algorithm, and does not replace compatibility tests.

After reviewing an intentional duplication or removing existing copies, refresh the baseline with `uvx --from jscpd==5.2.1 jscpd --config .jscpd.json --baseline .jscpd-baseline.json --update-baseline src crates` and review its diff. Do not refresh it just to make a failure disappear.

## Compatibility notes

Optparse-style abbreviations work as in UMI-tools, for example `--unmapped` for `--unmapped-reads`. Flags UMI-tools accepts that umi-tools-rs does not implement are rejected with an error rather than silently ignored. Known differences:

- The run summary goes to stderr, or to the `--log` file when given. UMI-tools writes its log to stdout by default.
- `--log` appends a shorter header than UMI-tools writes.
- `--error` captures umi-tools-rs notes and errors; htslib messages still go to stderr.
- `--compresslevel` defaults to 3 rather than 6.
- `--reference-filename` takes a local path; URL references are not fetched.
- `--input-options`, `--output-options`, `--temp-dir`, `--timeit`, `--timeit-name` and `--timeit-header` print a note and have no effect. `--plot-prefix` prints a note and no plots are generated. `group --multimapping-detection-method` prints a note because group keeps every read. `whitelist --ignore-read-pair-suffixes` and the alignment and position flags `count_tab` inherits from UMI-tools print a note because they cannot apply to a table.
- `--in-format` and `--in-sam` have no effect because the input format is detected from the file content.
- `--gene-transcript-map` reads a gene's transcripts in file order. UMI-tools iterates a set, so which read represents a UMI group there varies with the Python hash seed; counts agree.
- `count` needs `--gene-tag` or `--per-contig`, as UMI-tools does; earlier releases defaulted to `XF`.
- `group --unpaired-reads=output` writes each unpaired read twice, once ungrouped and once grouped, because UMI-tools yields it and then keeps processing it.
- Help output is rendered by clap and does not match the UMI-tools text.
