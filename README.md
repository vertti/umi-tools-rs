# umi-tools-rs

[![CI](https://github.com/vertti/umi-tools-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/vertti/umi-tools-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Crates.io](https://img.shields.io/crates/v/umi-tools-rs)](https://crates.io/crates/umi-tools-rs)

A drop-in replacement for [UMI-tools](https://github.com/CGATOxford/UMI-tools), written in Rust. Same flags, same output — just faster.

## Performance

Measured with [hyperfine](https://github.com/sharkdp/hyperfine) against UMI-tools 1.1.6:

| Command | Speedup |
|:--------|--------:|
| `extract` | **14-36x** |
| `whitelist` | **32x** |
| `dedup` | **50x** |
| `group` | **17x** |
| `count` | **80x** |
| `count_tab` | **91x** |

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

# Linux aarch64
curl -fsSL https://github.com/vertti/umi-tools-rs/releases/latest/download/umi-tools-rs-aarch64-unknown-linux-gnu.tar.gz \
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
