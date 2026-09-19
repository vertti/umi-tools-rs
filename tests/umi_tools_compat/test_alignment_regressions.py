"""Small alignment fixtures for regressions outside the upstream suite."""

import subprocess

import pysam
import pytest


def write_bam(path, tags):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as stream:
        read = pysam.AlignedSegment(stream.header)
        read.query_name = "r_AAAA"
        read.query_sequence = "ACGT"
        read.query_qualities = pysam.qualitystring_to_array("IIII")
        read.reference_id = 0
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigarstring = "4M"
        read.set_tags(tags)
        stream.write(read)
    pysam.index(str(path))


@pytest.mark.parametrize("tag", ["BX", "RX"])
def test_group_replaces_existing_annotations(rust_binary, tmp_path, tag):
    source, output, tsv = [tmp_path / name for name in ("in.bam", "out.bam", "groups.tsv")]
    write_bam(source, [("UG", 999), (tag, "OLD")])
    result = subprocess.run(
        [str(rust_binary), "group", "--stdin", str(source), "--output-bam",
         "--stdout", str(output), "--group-out", str(tsv), f"--umi-group-tag={tag}"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    with pysam.AlignmentFile(output) as stream:
        read, = list(stream)
    fields = tsv.read_text().splitlines()[1].split("\t")
    assert read.get_tag("UG") == int(fields[8]) == 0
    assert read.get_tag(tag) == fields[6] == "AAAA"
    assert [name for name, _ in read.tags].count("UG") == 1
    assert [name for name, _ in read.tags].count(tag) == 1


@pytest.mark.parametrize("command", ["count", "dedup", "group"])
def test_transcript_map_replaces_an_existing_mc_tag(rust_binary, tmp_path, command):
    source, output = tmp_path / "in.bam", tmp_path / "out.bam"
    write_bam(source, [("MC", "4M")])
    mapping = tmp_path / "genes.tsv"
    mapping.write_text("GENE1\tchr1\n")
    args = [str(rust_binary), command, "--stdin", str(source), "--per-gene",
            "--per-contig", "--gene-transcript-map", str(mapping)]
    if command != "count":
        args += ["--stdout", str(output)]
    if command == "group":
        args += ["--output-bam"]
    result = subprocess.run(args, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    if command == "count":
        assert result.stdout == "gene\tcount\nGENE1\t1\n"
    else:
        with pysam.AlignmentFile(output) as stream:
            read, = list(stream)
        assert read.get_tag("MC") == "GENE1"
        assert [name for name, _ in read.tags].count("MC") == 1
