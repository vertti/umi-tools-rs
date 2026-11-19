"""Read handling in extract and whitelist, checked against umi_tools master."""

import subprocess

import pytest

READS = "@a\nACGTACGTAC\n+\nIIIIIIIIII\n@b\nTTGGCCAATT\n+\nIIIIIIIIII\n"
SHORT = "@a\nACGTACGTAC\n+\nIIIIIIIIII\n@b\nACG\n+\nIII\n"
CELLS = "".join(f"@r{i}\n{cell}ACGTAC\n+\nIIIIIIIIII\n"
                for i, cell in enumerate(["AAAA"] * 3 + ["CCCC"] * 2 + ["GGGG"]))


def run(binary, tmp_path, options, **files):
    for name, content in files.items():
        (tmp_path / name).write_text(content)
    return subprocess.run([str(binary), *options], capture_output=True, text=True, cwd=tmp_path)


@pytest.mark.parametrize("options", [
    ["extract", "--bc-pattern=NNNNNN", "-I", "short"],
    ["extract", "--bc-pattern=NNNNNN", "--3prime", "-I", "short"],
    ["extract", "--bc-pattern=NNN", "--bc-pattern2=NNNNNN", "-I", "reads", "--read2-in=short"],
    ["whitelist", "--bc-pattern=CCCCNN", "-I", "short"],
    ["whitelist", "--bc-pattern=CCCNN", "--bc-pattern2=CCCCCC", "-I", "reads", "--read2-in=short"],
])
def test_read_shorter_than_string_pattern_aborts(rust_binary, tmp_path, options):
    result = run(rust_binary, tmp_path, options, reads=READS, short=SHORT)
    assert result.returncode == 1
    assert "panicked" not in result.stderr


def test_extract_subset_of_zero_reads_everything(rust_binary, tmp_path):
    result = run(rust_binary, tmp_path,
                 ["extract", "--bc-pattern=NNN", "-I", "reads", "--subset-reads=0"], reads=READS)
    assert result.returncode == 0, result.stderr
    assert result.stdout.count("@") == 2


def test_whitelist_subset_of_zero_reads_everything(rust_binary, tmp_path):
    options = ["whitelist", "--bc-pattern=CCCCNN", "-I", "cells", "--set-cell-number=2"]
    unlimited = run(rust_binary, tmp_path, [*options, "--subset-reads=0"], cells=CELLS)
    assert unlimited.returncode == 0, unlimited.stderr
    assert unlimited.stdout == "AAAA\t\t3\t\nCCCC\t\t2\t\n"
