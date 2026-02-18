#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="$(dirname "$0")/data"
PATTERN="CCCCCCCCCCCCCCCCNNNNNNNNNN"

# Build release binary
echo "Building release binary..."
cargo build --release --quiet

RUST_BIN="$(cargo metadata --format-version 1 --no-deps | python3 -c 'import sys,json; print(json.load(sys.stdin)["target_directory"])')/release/umi-tools-rs"

if [[ ! -x "$RUST_BIN" ]]; then
    echo "ERROR: rust binary not found at $RUST_BIN"
    exit 1
fi

UMI_TOOLS="uv run --group bench umi_tools"
if ! $UMI_TOOLS --version &>/dev/null; then
    echo "ERROR: umi_tools not available via uv"
    exit 1
fi

if ! command -v hyperfine &>/dev/null; then
    echo "ERROR: hyperfine not found on PATH"
    echo "Install with: mise install"
    exit 1
fi

echo ""
echo "Versions:"
echo "  umi-tools-rs: $("$RUST_BIN" --version 2>&1 || echo 'unknown')"
echo "  umi_tools:    $($UMI_TOOLS --version 2>&1 || echo 'unknown')"
echo "  hyperfine:    $(hyperfine --version)"
echo ""

run_benchmark() {
    local label=$1
    local input=$2

    echo "=== Benchmark: $label ==="
    echo "Input: $input ($(du -h "$input" | cut -f1))"
    echo ""

    hyperfine \
        --warmup 2 \
        --min-runs 5 \
        --export-markdown "benchmarks/results_${label}.md" \
        --command-name "umi-tools-rs" \
        "$RUST_BIN extract --bc-pattern=$PATTERN --stdin=$input --stdout=/dev/null" \
        --command-name "umi_tools (python)" \
        "$UMI_TOOLS extract --bc-pattern=$PATTERN --stdin=$input --stdout=/dev/null --log=/dev/null"

    echo ""
}

run_benchmark "100k" "$DATA_DIR/hgmm_100_R1_100k.fastq.gz"
run_benchmark "full" "$DATA_DIR/hgmm_100_R1.fastq.gz"

# --- Dedup benchmarks ---

SCRIPT_DIR="$(dirname "$0")"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
UMI_TOOLS_TESTS_DIR="${UMI_TOOLS_TESTS_DIR:-$REPO_ROOT/../umi-tools/tests}"
CHR19_BAM="$UMI_TOOLS_TESTS_DIR/chr19.bam"

if [[ -f "$CHR19_BAM" ]]; then
    run_dedup_benchmark() {
        local label=$1
        local method=$2

        echo "=== Benchmark: dedup_${label} ==="
        echo "Input: $CHR19_BAM ($(du -h "$CHR19_BAM" | cut -f1))"
        echo ""

        hyperfine \
            --warmup 2 \
            --min-runs 5 \
            --export-markdown "benchmarks/results_dedup_${label}.md" \
            --command-name "umi-tools-rs" \
            "$RUST_BIN dedup --method=$method --out-sam --random-seed=123456789 --stdin=$CHR19_BAM > /dev/null" \
            --command-name "umi_tools (python)" \
            "$UMI_TOOLS dedup --method=$method --out-sam --random-seed=123456789 --stdin=$CHR19_BAM --stdout=/dev/null --log=/dev/null"

        echo ""
    }

    run_dedup_benchmark "directional" "directional"

    # --- Group benchmarks ---

    echo "=== Benchmark: group_directional ==="
    echo "Input: $CHR19_BAM ($(du -h "$CHR19_BAM" | cut -f1))"
    echo ""

    hyperfine \
        --warmup 2 \
        --min-runs 5 \
        --export-markdown "benchmarks/results_group_directional.md" \
        --command-name "umi-tools-rs" \
        "$RUST_BIN group --method=directional --output-bam --stdin=$CHR19_BAM > /dev/null" \
        --command-name "umi_tools (python)" \
        "$UMI_TOOLS group --method=directional --output-bam --stdin=$CHR19_BAM --stdout=/dev/null --log=/dev/null"

    echo ""

    # --- Count benchmarks ---

    CHR19_GENE_BAM="$UMI_TOOLS_TESTS_DIR/chr19_gene_tags.bam"
    CHR19_GENE_TSV="$UMI_TOOLS_TESTS_DIR/chr19_gene_assigned.tsv"

    if [[ -f "$CHR19_GENE_BAM" ]]; then
        echo "=== Benchmark: count ==="
        echo "Input: $CHR19_GENE_BAM ($(du -h "$CHR19_GENE_BAM" | cut -f1))"
        echo ""

        hyperfine \
            --warmup 2 \
            --min-runs 5 \
            --export-markdown "benchmarks/results_count.md" \
            --command-name "umi-tools-rs" \
            "$RUST_BIN count --method=directional --gene-tag=XF --skip-tags-regex='^[__|Unassigned]' --extract-umi-method=umis --stdin=$CHR19_GENE_BAM > /dev/null" \
            --command-name "umi_tools (python)" \
            "$UMI_TOOLS count --method=directional --gene-tag=XF --skip-tags-regex='^[__|Unassigned]' --extract-umi-method=umis --stdin=$CHR19_GENE_BAM --stdout=/dev/null --log=/dev/null"

        echo ""
    else
        echo "Skipping count benchmark: $CHR19_GENE_BAM not found"
    fi

    if [[ -f "$CHR19_GENE_TSV" ]]; then
        echo "=== Benchmark: count_tab ==="
        echo "Input: $CHR19_GENE_TSV ($(du -h "$CHR19_GENE_TSV" | cut -f1))"
        echo ""

        hyperfine \
            --warmup 2 \
            --min-runs 5 \
            --export-markdown "benchmarks/results_count_tab.md" \
            --command-name "umi-tools-rs" \
            "$RUST_BIN count_tab --stdin=$CHR19_GENE_TSV > /dev/null" \
            --command-name "umi_tools (python)" \
            "$UMI_TOOLS count_tab --stdin=$CHR19_GENE_TSV --stdout=/dev/null --log=/dev/null"

        echo ""
    else
        echo "Skipping count_tab benchmark: $CHR19_GENE_TSV not found"
    fi

else
    echo "Skipping BAM-based benchmarks: $CHR19_BAM not found"
    echo "Set UMI_TOOLS_TESTS_DIR or clone umi-tools next to this repo"
fi

# --- Whitelist benchmarks ---

FASTQ_100K="$DATA_DIR/hgmm_100_R1_100k.fastq.gz"
if [[ -f "$FASTQ_100K" ]]; then
    echo "=== Benchmark: whitelist ==="
    echo "Input: $FASTQ_100K ($(du -h "$FASTQ_100K" | cut -f1))"
    echo ""

    hyperfine \
        --warmup 2 \
        --min-runs 5 \
        --export-markdown "benchmarks/results_whitelist.md" \
        --command-name "umi-tools-rs" \
        "$RUST_BIN whitelist --bc-pattern=$PATTERN --stdin=$FASTQ_100K --stdout=/dev/null" \
        --command-name "umi_tools (python)" \
        "$UMI_TOOLS whitelist --bc-pattern=$PATTERN --stdin=$FASTQ_100K --stdout=/dev/null --log=/dev/null"

    echo ""
fi

echo "Results written to benchmarks/results_*.md"
