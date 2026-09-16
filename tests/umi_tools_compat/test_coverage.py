"""Every upstream test is either enabled in a compat module or listed here as excluded."""

from . import test_count, test_dedup, test_extract, test_group, test_whitelist
from .conftest import _load_tests_yaml

MODULES = [test_count, test_dedup, test_extract, test_group, test_whitelist]

# These compare optparse's help text verbatim; clap renders help differently by design.
EXCLUDED_TESTS = {
    "umi_tools_help",
    "whitelist_help",
    "extract_help",
    "dedup_help",
    "group_help",
    "count_help",
    "count_tab_help",
}


def test_every_upstream_test_is_enabled_or_excluded(umi_tools_tests_dir):
    upstream = set(_load_tests_yaml())
    enabled = set().union(*(set(m.ENABLED_TESTS) for m in MODULES))

    assert not enabled & EXCLUDED_TESTS, "a test is both enabled and excluded"

    missing = upstream - enabled - EXCLUDED_TESTS
    assert not missing, f"upstream tests neither enabled nor excluded: {sorted(missing)}"

    stale = (enabled | EXCLUDED_TESTS) - upstream
    assert not stale, f"tests no longer present upstream: {sorted(stale)}"
