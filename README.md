# umi-tools-rs

A fast drop-in replacement for [UMI-tools](https://github.com/CGATOxford/UMI-tools), written in Rust.

**14-33x faster** than Python UMI-tools on real 10x Genomics data.

## Benchmarks

Measured with [hyperfine](https://github.com/sharkdp/hyperfine) on 10x Genomics hgmm_100 data (`extract`, pattern `CCCCCCCCCCCCCCCCNNNNNNNNNN`):

| Dataset | umi-tools-rs | UMI-tools (Python) | Speedup |
|:--------|-------------:|-------------------:|--------:|
| 100K reads (1.8 MB gz) | 26 ms | 853 ms | **33x** |
| 912K reads (19 MB gz) | 207 ms | 2,964 ms | **14x** |

Run benchmarks yourself: `mise run bench`

## Feature support

| Command | Description | Status |
|:--------|:------------|:------:|
| `extract` | Extract UMIs from FASTQ reads | Done |
| `whitelist` | Build cell barcode whitelist | Done |
| `dedup` | Deduplicate aligned BAM reads | Done |
| `group` | Group PCR duplicates in BAM | Planned |
| `count` | Count unique molecules per gene | Planned |
| `count_tab` | Count from flatfile input | Planned |

## Installation

```sh
cargo install --path .
```

Or build from source:

```sh
cargo build --release
# Binary at target/release/umi-tools-rs
```

## Usage

```sh
# Extract UMIs (single-end)
umi-tools-rs extract --bc-pattern=NNNNNNNN --stdin=reads.fastq.gz --stdout=extracted.fastq.gz

# Extract UMIs (paired-end)
umi-tools-rs extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --stdin=R1.fastq.gz --read2-in=R2.fastq.gz \
  --stdout=R1_extracted.fastq.gz --read2-out=R2_extracted.fastq.gz

# Whitelist cell barcodes
umi-tools-rs whitelist --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --stdin=R1.fastq.gz --stdout=whitelist.tsv

# Extract with cell barcode filtering
umi-tools-rs extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
  --stdin=R1.fastq.gz --stdout=extracted.fastq.gz \
  --whitelist=whitelist.tsv --error-correct-cell

# Deduplicate BAM reads
umi-tools-rs dedup --stdin=aligned.bam --stdout=deduped.bam

# Dedup with specific method and chromosome filter
umi-tools-rs dedup --method=directional --chrom=chr19 \
  --stdin=aligned.bam --stdout=deduped.bam
```

Aims to be CLI-compatible with Python UMI-tools — same flags, same output format.
