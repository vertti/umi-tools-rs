"""Compatibility tests: run umi-tools-rs dedup against the real umi-tools test suite."""

import pytest

from .conftest import run_case, upstream_params

ENABLED_TESTS = [
    "dedup_single_ignore",
    "dedup_single_chrom",
    "dedup_single_unique",
    "dedup_single_perc",
    "dedup_single_cluster",
    "dedup_single_adj",
    "dedup_single_dir",
    "dedup_single_sep",
    "dedup_single_dir_edit_dist",
    "dedup_single_tag",
    "dedup_single_subset",
    "dedup_single_gene_tag",
    "dedup_single_tag_missing",
    "dedup_single_stats",
    "dedup_paired_ignore_tlen_tag",
    "dedup_paired_umi_whitelist",
    "dedup_from_cram",
    "dedup_to_cram",
    "dedup_bam_to_cram",
]


@pytest.mark.parametrize("case", upstream_params(ENABLED_TESTS))
def test_dedup(rust_binary, umi_tools_tests_dir, case):
    run_case(rust_binary, case)
