"""Compatibility tests: run umi-tools-rs group against the real umi-tools test suite."""

import pytest

from .conftest import run_case, upstream_params

ENABLED_TESTS = [
    "group_unique",
    "group_cluster",
    "group_adjacency",
    "group_directional",
    "group_unsorted",
    "group_directional_subset",
    "group_directional_unmapped",
    "group_gene_tag",
    "group_contig_no_gene_tag",
    "group_paired_discard_chimeric",
    "group_paired_output_chimeric",
    "group_paired_use_chimeric",
    "group_paired_discard_unmapped",
    "group_paired_output_unmapped",
    "group_paired_use_unmapped",
]


@pytest.mark.parametrize("case", upstream_params(ENABLED_TESTS))
def test_group(rust_binary, umi_tools_tests_dir, case):
    run_case(rust_binary, case)
