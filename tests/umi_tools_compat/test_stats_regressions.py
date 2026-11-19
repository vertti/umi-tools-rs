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


def test_single_umi_positions_add_no_edit_distance_row(rust_binary, tmp_path):
    run_stats(rust_binary, tmp_path, {"AAAA": 2})
    assert (tmp_path / "stats_edit_distance.tsv").read_text() == (
        "unique\tunique_null\tadjacency\tadjacency_null\tedit_distance\n"
        "1\t1\t1\t1\tSingle_UMI\n"
        "0\t0\t0\t0\t0\n"
    )


MIXED_LENGTHS = "not all umis are the same length(!):  4 - 5"


@pytest.mark.parametrize("method", ["unique", "percentile", "cluster", "adjacency", "directional"])
@pytest.mark.parametrize("command", [
    ["dedup"], ["dedup", "--output-stats=stats"], ["group", "--group-out=groups.tsv", "--output-bam"],
    ["count", "--per-contig"],
])
def test_mixed_umi_lengths_in_a_bundle_abort(rust_binary, tmp_path, command, method):
    source = tmp_path / "in.sam"
    write_counts(source, {"ACGT": 1, "ACGTA": 1})
    result = subprocess.run(
        [str(rust_binary), *command, "--stdin", str(source), "--stdout", "out", f"--method={method}"],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert result.returncode == 1
    assert MIXED_LENGTHS in result.stderr


def test_count_tab_rejects_mixed_umi_lengths(rust_binary, tmp_path):
    source = tmp_path / "in.tsv"
    source.write_text("r1_ACGT\tg1\nr2_ACGTA\tg1\n")
    result = subprocess.run(
        [str(rust_binary), "count_tab", "--stdin", str(source)], capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert MIXED_LENGTHS in result.stderr


def test_null_umis_of_mixed_lengths_fail_statistics(rust_binary, tmp_path):
    source = tmp_path / "in.sam"
    write_counts(source, {"AAAA": 1, "AAAT": 1})
    longer = "".join(f"s{i}_{umi}\t0\tchr1\t301\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                     for i, umi in enumerate(["CCCCC", "CCCCG"]))
    source.write_text(source.read_text() + longer)
    result = subprocess.run(
        [str(rust_binary), "dedup", "--stdin", str(source), "--stdout", "out.bam",
         "--output-stats=stats", "--random-seed=1"], capture_output=True, text=True, cwd=tmp_path,
    )
    assert result.returncode == 1
    assert "panicked" not in result.stderr


def test_dedup_reports_processed_positions(rust_binary, tmp_path):
    source = tmp_path / "positions.sam"
    write_counts(source, {"AAAA": 2, "AAAT": 1})
    source.write_text(source.read_text().replace("r0_AAAT\t0\tchr1\t101", "r0_AAAT\t0\tchr1\t201"))
    result = subprocess.run(
        [str(rust_binary), "dedup", "--stdin", str(source), "--stdout", str(tmp_path / "out.bam")],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "positions: 2" in result.stderr


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


BUNDLES = [(101, {"AAAA": 5, "AAAT": 2, "CCCC": 1}), (301, {"GGGG": 3, "GGGT": 1}),
           (601, {"TTTT": 2}), (901, {"ACGT": 1, "TGCA": 1, "AAAA": 1})]


def bundle_stats(binary, tmp_path, *options):
    """Statistics files for BUNDLES; expected contents come from umi_tools master."""
    names = (f"q{i}" for i in range(100))
    rows = [f"{next(names)}_{umi}\t0\tchr1\t{position}\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            for position, counts in BUNDLES for umi, count in counts.items() for _ in range(count)]
    source = tmp_path / "in.sam"
    source.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n" + "".join(rows))
    (tmp_path / "umis.txt").write_text("AAAA\nCCCC\n")
    result = subprocess.run(
        [str(binary), "dedup", "--stdin", str(source), "--stdout", "out.bam",
         "--output-stats=stats", "--random-seed=1", *options],
        capture_output=True, text=True, cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    return {suffix: (tmp_path / f"stats_{suffix}.tsv").read_text()
            for suffix in ("edit_distance", "per_umi", "per_umi_per_position")}


def test_unique_method_reports_its_own_tallies_in_the_unique_columns(rust_binary, tmp_path):
    stats = bundle_stats(rust_binary, tmp_path, "--method=unique")
    assert stats["edit_distance"] == "unique\tunique_null\tunique\tunique_null\tedit_distance\n1\t1\t1\t1\tSingle_UMI\n0\t1\t0\t1\t0\n1\t0\t1\t0\t1\n0\t0\t0\t0\t2\n1\t1\t1\t1\t3\n1\t1\t1\t1\t4\n0\t0\t0\t0\t5\n"


def test_bundles_emptied_by_the_umi_whitelist_leave_no_statistics(rust_binary, tmp_path):
    stats = bundle_stats(rust_binary, tmp_path, "--filter-umi", "--umi-whitelist=umis.txt")
    assert stats["edit_distance"] == "unique\tunique_null\tdirectional\tdirectional_null\tedit_distance\n0\t0\t1\t1\tSingle_UMI\n0\t1\t0\t1\t0\n0\t0\t0\t0\t1\n0\t0\t0\t0\t2\n1\t1\t0\t0\t3\n1\t0\t1\t0\t4\n0\t0\t0\t0\t5\n"
    assert stats["per_umi"] == "UMI\tmedian_counts_pre\ttimes_observed_pre\ttotal_counts_pre\tmedian_counts_post\ttimes_observed_post\ttotal_counts_post\nAAAA\t3\t2\t6\t4\t2\t8\nAAAT\t2\t1\t2\t0\t0\t0\nACGT\t1\t1\t1\t0\t0\t0\nCCCC\t1\t1\t1\t1\t1\t1\nTGCA\t1\t1\t1\t0\t0\t0\n"
    assert stats["per_umi_per_position"] == "counts\tinstances_pre\tinstances_post\n1\t4\t2\n2\t1\t0\n5\t1\t0\n7\t0\t1\n"
