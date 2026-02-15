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

echo "Results written to benchmarks/results_*.md"
