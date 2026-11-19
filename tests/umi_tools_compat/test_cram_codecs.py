"""CRAM blocks may be bzip2 or lzma compressed; both must decode."""

import subprocess

import pysam
import pytest


@pytest.mark.parametrize("codec", ["use_bzip2", "use_lzma"])
def test_dedup_reads_cram_codec(rust_binary, umi_tools_tests_dir, tmp_path, codec):
    reference = umi_tools_tests_dir / "chr19_1mb.fa"
    source = umi_tools_tests_dir / "chr19_1mb.bam"
    cram = tmp_path / "in.cram"
    pysam.view("-C", "-T", str(reference), "--output-fmt-option", f"{codec}=1",
               "-o", str(cram), str(source), catch_stdout=False)
    pysam.index(str(cram))

    def dedup(path):
        result = subprocess.run(
            [str(rust_binary), "dedup", "--stdin", str(path), "--out-sam", "--random-seed=1",
             f"--reference-filename={reference}"], capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        return [line for line in result.stdout.splitlines() if not line.startswith("@")]

    assert dedup(cram) == dedup(source)
