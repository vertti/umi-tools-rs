#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="$(dirname "$0")/data"
TAR_URL="https://cf.10xgenomics.com/samples/cell-exp/1.3.0/hgmm_100/hgmm_100_fastqs.tar"
R1_PATH_IN_TAR="fastqs/hgmm_100_S1_L001_R1_001.fastq.gz"
FULL_GZ="$DATA_DIR/hgmm_100_R1.fastq.gz"

mkdir -p "$DATA_DIR"

# Download and extract R1 if not already present
if [[ ! -f "$FULL_GZ" ]]; then
    echo "Downloading hgmm_100 tar (~798 MB)..."
    TAR_FILE="$DATA_DIR/hgmm_100_fastqs.tar"
    curl -L -o "$TAR_FILE" "$TAR_URL"

    echo "Extracting R1 FASTQ..."
    tar -xf "$TAR_FILE" -C "$DATA_DIR" "$R1_PATH_IN_TAR"
    mv "$DATA_DIR/$R1_PATH_IN_TAR" "$FULL_GZ"
    rmdir "$DATA_DIR/fastqs"
    rm "$TAR_FILE"

    echo "Verifying R1..."
    gzip -t "$FULL_GZ"
    echo "R1 OK: $FULL_GZ"
else
    echo "R1 already exists: $FULL_GZ"
fi

# Create subsets
create_subset() {
    local reads=$1 label=$2
    local lines=$((reads * 4))
    local out="$DATA_DIR/hgmm_100_R1_${label}.fastq.gz"

    if [[ -f "$out" ]]; then
        echo "Subset already exists: $out"
        return
    fi

    echo "Creating ${label} subset (${reads} reads)..."
    # head causes SIGPIPE on the gzip writer; tolerate it
    (gzip -dc "$FULL_GZ" || true) | head -n "$lines" | gzip > "$out"
    gzip -t "$out"
    echo "Subset OK: $out"
}

create_subset 100000 "100k"

echo "Setup complete. Files in $DATA_DIR:"
ls -lh "$DATA_DIR"
