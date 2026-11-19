"""Seed and subset edge cases for dedup, group and count, checked against umi_tools master."""

import subprocess

import pytest

UMIS = ["AAAA", "CCCC", "GGGG", "TTTT", "ACGT", "TGCA", "AACC", "GGTT"]


@pytest.fixture
def source(tmp_path):
    header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
    rows = [f"r{i}_{UMIS[i % 8]}\t0\tchr1\t{101 + i * 10}\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            for i in range(40)]
    path = tmp_path / "in.sam"
    path.write_text(header + "".join(rows))
    return path


def run(binary, source, command, *options):
    return subprocess.run([str(binary), command, "-I", str(source), *options],
                          capture_output=True, text=True, cwd=source.parent)


def records(sam):
    return [line for line in sam.splitlines() if not line.startswith("@")]


def test_group_seeds_python_random_with_all_bits_of_a_large_seed(rust_binary, source):
    result = run(rust_binary, source, "group", "--group-out=groups.tsv", "--subset=0.5",
                 "--random-seed=4294967297")
    assert result.returncode == 0, result.stderr
    rows = (source.parent / "groups.tsv").read_text().splitlines()[1:]
    assert " ".join(row.split("\t")[0] for row in rows) == (
        "r0_AAAA r2_GGGG r3_TTTT r7_GGTT r9_CCCC r10_GGGG r17_CCCC r19_TTTT r20_ACGT r22_AACC r25_CCCC r26_GGGG r30_AACC r33_CCCC r36_ACGT r37_TGCA r38_AACC r39_GGTT"
    )


@pytest.mark.parametrize("command", [["dedup", "--out-sam"], ["count", "--per-contig", "--per-gene"]])
def test_numpy_seed_range_is_enforced(rust_binary, source, command):
    result = run(rust_binary, source, *command, "--random-seed=4294967296")
    assert result.returncode == 1
    assert "Seed must be between 0 and 2**32 - 1" in result.stderr
    assert run(rust_binary, source, *command, "--random-seed=4294967295").returncode == 0


def test_subset_of_zero_keeps_every_read(rust_binary, source):
    result = run(rust_binary, source, "dedup", "--out-sam", "--subset=0")
    assert result.returncode == 0, result.stderr
    assert len(records(result.stdout)) == 40
    grouped = run(rust_binary, source, "group", "--group-out=groups.tsv", "--subset=0")
    assert grouped.returncode == 0, grouped.stderr
    assert len((source.parent / "groups.tsv").read_text().splitlines()) == 41
    counted = run(rust_binary, source, "count", "--per-contig", "--per-gene", "--subset=0")
    assert counted.stdout == "gene\tcount\nchr1\t8\n"


def test_empty_umi_whitelist_filters_nothing(rust_binary, source):
    (source.parent / "empty.txt").write_text("")
    result = run(rust_binary, source, "dedup", "--out-sam", "--filter-umi",
                 "--umi-whitelist=empty.txt")
    assert result.returncode == 0, result.stderr
    assert len(records(result.stdout)) == 40
