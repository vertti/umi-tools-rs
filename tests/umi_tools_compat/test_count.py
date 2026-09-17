"""Compatibility tests: run umi-tools-rs count/count_tab against the real umi-tools test suite."""

import pytest

from .conftest import run_case, upstream_params

ENABLED_TESTS = [
    "count_single_gene_tag",
    "count_single_cells_gene_tag",
    "count_single_cells_wide_gene_tag",
    "count_tab_single",
    "count_tab_single_per_cell",
]


@pytest.mark.parametrize("case", upstream_params(ENABLED_TESTS))
def test_count(rust_binary, umi_tools_tests_dir, case):
    run_case(rust_binary, case)
