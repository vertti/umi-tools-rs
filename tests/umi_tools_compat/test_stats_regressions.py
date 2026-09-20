"""Statistics must conserve reads and handle empty eligible input."""

import csv
import subprocess

import pytest


def write_counts(path, counts):
    header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
    rows = [f"r{i}_{umi}\t0\tchr1\t101\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            for umi, count in counts.items() for i in range(count)]
    path.write_text(header + "".join(rows))


def run_stats(binary, tmp_path, counts, method="adjacency"):
    source = tmp_path / "in.sam"
    write_counts(source, counts)
    result = subprocess.run(
        [str(binary), "dedup", "--stdin", str(source), "--stdout", str(tmp_path / "out.bam"),
         "--output-stats", str(tmp_path / "stats"), f"--method={method}"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    with (tmp_path / "stats_per_umi.tsv").open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_adjacency_cluster_totals_conserve_reads(rust_binary, tmp_path):
    rows = run_stats(rust_binary, tmp_path, {"AAAA": 10, "AAAT": 9, "AATT": 1})
    totals = {row["UMI"]: int(row["total_counts_post"]) for row in rows}
    assert totals == {"AAAA": 10, "AAAT": 10, "AATT": 0}
    assert sum(totals.values()) == 20


def test_empty_input_produces_empty_statistics(rust_binary, tmp_path):
    assert run_stats(rust_binary, tmp_path, {}) == []
    for suffix in ("per_umi_per_position", "edit_distance"):
        assert (tmp_path / f"stats_{suffix}.tsv").read_text().strip()


@pytest.mark.parametrize("method,expected", [
    ("unique", 3), ("percentile", 3), ("cluster", 1), ("adjacency", 2), ("directional", 2),
])
def test_clustering_agrees_across_commands(rust_binary, tmp_path, method, expected):
    rows = run_stats(rust_binary, tmp_path, {"AAAA": 10, "AAAT": 9, "AATT": 1}, method)
    assert sum(int(row["times_observed_post"]) for row in rows) == expected
    assert sum(int(row["total_counts_post"]) for row in rows) == 20
    source = tmp_path / "in.sam"
    tsv = tmp_path / "groups.tsv"
    grouped = subprocess.run(
        [str(rust_binary), "group", "--stdin", str(source), "--group-out", str(tsv),
         f"--method={method}"], capture_output=True, text=True,
    )
    assert grouped.returncode == 0, grouped.stderr
    with tsv.open() as stream:
        groups = list(csv.DictReader(stream, delimiter="\t"))
    assert len({row["unique_id"] for row in groups}) == expected
    assert len(groups) == 20
    counted = subprocess.run(
        [str(rust_binary), "count", "--stdin", str(source), "--per-contig", f"--method={method}"],
        capture_output=True, text=True,
    )
    assert counted.returncode == 0, counted.stderr
    assert counted.stdout == f"gene\tcount\nchr1\t{expected}\n"
