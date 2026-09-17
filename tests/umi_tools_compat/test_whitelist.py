"""Compatibility tests: run umi-tools-rs whitelist against the real umi-tools test suite."""

import pytest

from .conftest import run_case, upstream_params

ENABLED_TESTS = [
    "whitelist_scrb_seq",
    "whitelist_indrop",
    "whitelist_indrop_set_cell",
    "whitelist_indrop_3_errors",
    "whitelist_indrop_density",
    "whitelist_indrop_expect_cells_density",
    "whitelist_indrop_filtered_out",
    "whitelist_indrop_ed_above_threshold_discard",
    "whitelist_indrop_ed_above_threshold_correct",
]


@pytest.mark.parametrize("case", upstream_params(ENABLED_TESTS))
def test_whitelist(rust_binary, umi_tools_tests_dir, case):
    run_case(rust_binary, case)
