"""Compatibility tests: run umi-tools-rs against the real umi-tools test suite."""

import pytest

from .conftest import run_case, upstream_params

ENABLED_TESTS = [
    "extract_single_string",
    "extract_single",
    "extract_3prime",
    "extract_quality",
    "extract_read2_only_string",
    "extract_read2_only_regex",
    "extract_scrb_seq",
    "extract_scrb_seq_string",
    "extract_scrb_seq_suffix",
    "extract_scrb_seq_prefiltered",
    "extract_indrop_blacklist",
    "extract_indrop_fuzzy",
    "extract_indrop_output_filtered",
    "extract_either_read",
]


@pytest.mark.parametrize("case", upstream_params(ENABLED_TESTS))
def test_extract(rust_binary, umi_tools_tests_dir, case):
    run_case(rust_binary, case)
