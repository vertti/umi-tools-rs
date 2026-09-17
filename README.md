# umi-tools-rs

[![CI](https://github.com/vertti/umi-tools-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/vertti/umi-tools-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Crates.io](https://img.shields.io/crates/v/umi-tools-rs)](https://crates.io/crates/umi-tools-rs)

A drop-in replacement for [UMI-tools](https://github.com/CGATOxford/UMI-tools), written in Rust. Same flags, same output — just faster.

## Performance

Measured with [hyperfine](https://github.com/sharkdp/hyperfine) against UMI-tools 1.1.6:

| Command | Speedup |
|:--------|--------:|
| `extract` | **14-43x** |
| `whitelist` | **31x** |
| `dedup` | **50x** |
| `group` | **17x** |
| `count` | **66x** |
| `count_tab` | **103x** |

Run benchmarks yourself: `mise run bench`

## Installation

### From crates.io

```sh
cargo install umi-tools-rs
```

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

## Compatibility notes

Optparse-style abbreviations work as in UMI-tools, for example `--unmapped` for `--unmapped-reads`. Flags UMI-tools accepts that umi-tools-rs does not implement are rejected with an error rather than silently ignored. Known differences:

- The run summary goes to stderr, or to the `--log` file when given. UMI-tools writes its log to stdout by default.
- `--log` appends a shorter header than UMI-tools writes.
- `--error` captures umi-tools-rs notes and errors; htslib messages still go to stderr.
- `--compresslevel` defaults to 3 rather than 6.
- `--reference-filename` takes a local path; URL references are not fetched.
- `--input-options`, `--output-options`, `--temp-dir`, `--timeit`, `--timeit-name` and `--timeit-header` print a note and have no effect. `--plot-prefix` prints a note and no plots are generated. `group --multimapping-detection-method` prints a note because group keeps every read.
- `--in-format` and `--in-sam` have no effect because the input format is detected from the file content.
- Help output is rendered by clap and does not match the UMI-tools text.
