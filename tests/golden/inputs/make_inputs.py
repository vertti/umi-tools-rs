"""Recreate the golden input BAMs from the upstream umi-tools test data.

    uv run --group bench python tests/golden/inputs/make_inputs.py

Subsampled inputs keep the golden references small; tags_sub.bam adds the
cell and UMI encodings the upstream data lacks (a separator-delimited read
name, RX/CB tags with delimiters and a 10x-style GEM suffix).
"""

import os
import subprocess
import sys
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
UPSTREAM = Path(
    os.environ.get("UMI_TOOLS_TESTS_DIR", HERE.parent.parent.parent.parent / "umi-tools" / "tests")
)

SUBSAMPLES = {
    "paired_sub.bam": ("paired.bam", "7.15"),
    "whitelist_umi_sub.bam": ("whitelist_umi_input.bam", "7.15"),
    "chr19_sub.bam": ("chr19.bam", "7.08"),
}


def subsample(target, source, fraction):
    with open(HERE / target, "wb") as out:
        subprocess.run(
            ["samtools", "view", "-b", "-s", fraction, str(UPSTREAM / source)],
            stdout=out,
            check=True,
        )


def make_tags(target, source):
    """Re-encode umis-style names (…:CELL_x:UMI_y:…) as name#CELL#UMI plus RX/CB tags.

    The RX delimiter sits at a read-dependent position so that removing it with
    --umi-tag-delimiter changes the grouping; RY carries a read-dependent suffix
    so that cutting it off with --umi-tag-split changes the grouping.
    """
    with pysam.AlignmentFile(str(UPSTREAM / source)) as inp:
        with pysam.AlignmentFile(str(HERE / target), "wb", template=inp) as out:
            for read in inp:
                fields = dict(f.split("_", 1) for f in read.query_name.split(":")[7:])
                cell, umi = fields["CELL"], fields["UMI"]
                prefix = ":".join(read.query_name.split(":")[:7])
                read.query_name = f"{prefix}#{cell}#{umi}"
                digest = sum(read.query_name.encode())
                cut = 2 + digest % 3
                read.set_tag("RX", f"{umi[:cut]}-{umi[cut:]}")
                read.set_tag("RY", f"{umi}-{1 + digest % 2}")
                read.set_tag("CB", f"{cell}-1")
                out.write(read)


def main():
    if not UPSTREAM.exists():
        sys.exit(f"upstream tests directory not found at {UPSTREAM}")
    for target, (source, fraction) in SUBSAMPLES.items():
        subsample(target, source, fraction)
    make_tags("tags_sub.bam", "chr19_gene_tags.bam")
    for bam in sorted(HERE.glob("*.bam")):
        pysam.index(str(bam))
        print(f"{bam.name}: {pysam.AlignmentFile(str(bam)).count(until_eof=True)} reads")


if __name__ == "__main__":
    main()
