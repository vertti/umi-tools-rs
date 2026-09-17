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
    so that cutting it off with --umi-tag-split changes the grouping. XS holds
    the featureCounts-style assignment status and XT a gene for every read, so
    --assigned-status-tag=XS changes which reads are skipped.
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
                assigned = not read.get_tag("XF").startswith(("Unassigned", "__"))
                read.set_tag("XT", read.get_tag("XF") if assigned else "ENSG_UNASSIGNED")
                read.set_tag("XS", "Assigned" if assigned else read.get_tag("XF"))
                out.write(read)


def make_transcripts(target, source, map_target, bins=12):
    """Split single-contig reads into position bins and present each bin as a transcript.

    Gives --per-contig several contigs to count over, and the gene map groups
    transcripts into genes for --gene-transcript-map. The map also lists a
    transcript absent from the BAM, which umi_tools ignores.
    """
    with pysam.AlignmentFile(str(HERE / source)) as inp:
        reads = sorted(inp.fetch(until_eof=True), key=lambda r: r.reference_start)
    per_bin = -(-len(reads) // bins)
    chunks = [reads[i : i + per_bin] for i in range(0, len(reads), per_bin)]
    names = [f"ENST{i + 1:011d}" for i in range(len(chunks))]
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": name, "LN": chunk[-1].reference_end - chunk[0].reference_start + 1000}
            for name, chunk in zip(names, chunks)
        ],
    }
    with pysam.AlignmentFile(str(HERE / target), "wb", header=header) as out:
        for tid, chunk in enumerate(chunks):
            offset = chunk[0].reference_start
            for read in chunk:
                read.reference_id = tid
                read.reference_start -= offset
                out.write(read)
    genes = ["ENSG_A", "ENSG_A", "ENSG_A", "ENSG_B", "ENSG_C", "ENSG_C",
             "ENSG_D", "ENSG_D", "ENSG_D", "ENSG_D", "ENSG_E", "ENSG_F"]
    with open(HERE / map_target, "w") as out:
        out.write("# gene\ttranscript\n")
        for gene, name in zip(genes, names):
            out.write(f"{gene}\t{name}\n")
        out.write("ENSG_F\tENST99999999999\n")
    # One transcript per gene: umi_tools' output is then independent of its
    # hash-seeded set order, so dedup and group can have golden references.
    with open(HERE / map_target.replace(".tsv", "_single.tsv"), "w") as out:
        for i, name in enumerate(names):
            out.write(f"ENSG_S{i + 1:02d}\t{name}\n")


def make_unpaired(target, source):
    """Strip the pairing flags from roughly a tenth of the templates.

    Gives --paired runs some read1s without the paired flag, which exercise
    --unpaired-reads; upstream data has none.
    """
    with pysam.AlignmentFile(str(HERE / source)) as inp:
        with pysam.AlignmentFile(str(HERE / target), "wb", template=inp) as out:
            for read in inp:
                if sum(read.query_name.encode()) % 10 == 0:
                    read.flag &= ~(0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80)
                    read.next_reference_id = -1
                    read.next_reference_start = -1
                    read.template_length = 0
                out.write(read)


def main():
    if not UPSTREAM.exists():
        sys.exit(f"upstream tests directory not found at {UPSTREAM}")
    for target, (source, fraction) in SUBSAMPLES.items():
        subsample(target, source, fraction)
    make_tags("tags_sub.bam", "chr19_gene_tags.bam")
    make_transcripts("transcripts_sub.bam", "chr19_sub.bam", "gene_transcript_map.tsv")
    make_unpaired("paired_mixed_sub.bam", "paired_sub.bam")
    for bam in sorted(HERE.glob("*.bam")):
        pysam.index(str(bam))
        print(f"{bam.name}: {pysam.AlignmentFile(str(bam)).count(until_eof=True)} reads")


if __name__ == "__main__":
    main()
