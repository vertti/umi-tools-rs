"""Statistics must conserve reads and handle empty eligible input."""

import csv
import subprocess


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
