"""Compatibility tests: run umi-tools-rs dedup against the real umi-tools test suite."""

import gzip
import os
import re
import subprocess
import tempfile

import pytest

from .conftest import _get_tests_dir, _load_tests_yaml

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

# Upstream points the CRAM tests at a reference on GitHub; use the local copy.
REFERENCE_URL_RE = re.compile(r"(--reference-file\S*=)https://\S+/(\S+\.fa)")


def _read_cram(path, tests_dir):
    """Decode CRAM to header text plus SAM records, like the upstream harness.

    The reference is looked up next to the test data by the basename of the
    UR field, so the GitHub URL upstream records in the header is never fetched.
    """
    import pysam

    with pysam.AlignmentFile(path) as f:
        urls = {sq["UR"] for sq in f.header.to_dict()["SQ"]}
    assert len(urls) == 1, "harness supports a single reference per CRAM"
    reference = os.path.join(tests_dir, os.path.basename(urls.pop()))
    with pysam.AlignmentFile(path, reference_filename=reference) as f:
        return [f.text] + [read.to_string() for read in f.fetch(until_eof=True)]


def _read(path, cram=False, tests_dir=None):
    """Read file, decode, strip comment lines (matching umi-tools test logic)."""
    if cram:
        return _read_cram(path, tests_dir)
    if path.endswith(".gz"):
        with gzip.open(path) as f:
            data = f.read()
    else:
        with open(path, "rb") as f:
            data = f.read()

    try:
        text = data.decode("ascii")
    except UnicodeDecodeError:
        return data

    return [line for line in text.splitlines() if not line.startswith("#")]


def _build_params():
    """Build pytest parameters from tests.yaml for enabled dedup tests."""
    tests_yaml = _load_tests_yaml()
    if not tests_yaml:
        return []

    params = []
    for name in ENABLED_TESTS:
        if name not in tests_yaml:
            continue
        values = tests_yaml[name]
        params.append(
            pytest.param(
                name,
                values.get("stdin"),
                values["options"],
                values["outputs"],
                values["references"],
                values.get("sort", False),
                id=name,
            )
        )
    return params


@pytest.mark.parametrize("test_name,stdin,options,outputs,references,sort", _build_params())
def test_dedup(
    rust_binary, umi_tools_tests_dir, test_name, stdin, options, outputs, references, sort
):
    tmpdir = tempfile.mkdtemp()
    stdout_path = os.path.join(tmpdir, "stdout")

    # Build --stdin flag
    stdin_flag = ""
    if stdin:
        stdin_flag = f"--stdin={os.path.join(umi_tools_tests_dir, stdin)}"

    # Substitute directory placeholders
    opts = options
    opts = opts.replace("<DIR>", str(umi_tools_tests_dir))
    opts = opts.replace("%DIR%", str(umi_tools_tests_dir))
    opts = opts.replace("<TMP>", tmpdir)
    opts = opts.replace("%TMP%", tmpdir)
    opts = re.sub(r"\n", "", opts)

    opts = REFERENCE_URL_RE.sub(lambda m: f"{m.group(1)}{umi_tools_tests_dir}/{m.group(2)}", opts)

    statement = f"/bin/bash -c '{rust_binary} {opts} {stdin_flag} > {stdout_path}'"

    result = subprocess.run(
        statement, shell=True, capture_output=True, cwd=tmpdir
    )

    assert result.returncode == 0, (
        f"Command failed (rc={result.returncode}):\n"
        f"  cmd: {statement}\n"
        f"  stderr: {result.stderr.decode(errors='replace')}"
    )

    # Compare outputs against references
    for output_name, ref_name in zip(outputs, references):
        if output_name == "stdout":
            output_path = stdout_path
        elif output_name.startswith("<DIR>/") or output_name.startswith("%DIR%/"):
            output_path = os.path.join(str(umi_tools_tests_dir), output_name[6:])
        else:
            output_path = os.path.join(tmpdir, output_name)

        ref_path = os.path.join(str(umi_tools_tests_dir), ref_name)

        assert os.path.exists(output_path), f"Output file missing: {output_path}"
        assert os.path.exists(ref_path), f"Reference file missing: {ref_path}"

        cram = ref_name.endswith(".cram")
        actual = _read(output_path, cram, umi_tools_tests_dir)
        expected = _read(ref_path, cram, umi_tools_tests_dir)

        if sort:
            actual = sorted(actual)
            expected = sorted(expected)

        if actual != expected:
            diffs = []
            for a, b in zip(actual, expected):
                if a != b:
                    diffs.append(f"  got:    {a}\n  expect: {b}")
                    if len(diffs) >= 10:
                        break

            diff_str = "\n---\n".join(diffs)
            pytest.fail(
                f"Output mismatch for {test_name} ({output_name} vs {ref_name}):\n"
                f"  output lines: {len(actual)}, reference lines: {len(expected)}\n"
                f"First differences:\n{diff_str}"
            )
