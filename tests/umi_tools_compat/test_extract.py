"""Compatibility tests: run umi-tools-rs against the real umi-tools test suite."""

import gzip
import os
import re
import subprocess
import tempfile

import pytest

from .conftest import _get_tests_dir, _load_tests_yaml

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

# Flags the Rust binary doesn't support yet — stripped before invocation.
UNSUPPORTED_FLAG_RE = re.compile(r"--log=\S+|-L\s+\S+")


def _read(path):
    """Read file, decode, strip comment lines (matching umi-tools test logic)."""
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
    """Build pytest parameters from tests.yaml for enabled tests."""
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
def test_extract(
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

    # Strip unsupported flags
    opts = UNSUPPORTED_FLAG_RE.sub("", opts)

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

        actual = _read(output_path)
        expected = _read(ref_path)

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
