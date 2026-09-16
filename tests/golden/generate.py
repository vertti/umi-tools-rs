"""Regenerate the golden references in this directory with Python umi_tools.

    uv run --group bench python tests/golden/generate.py [CASE ...]

Inputs are read from the upstream umi-tools tests directory (or
UMI_TOOLS_TESTS_DIR). Each case in cases.yaml is run exactly as umi-tools-rs
runs it, except that umi_tools' log is sent to a file so stdout holds only the
command's output.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
sys.path.insert(0, str(REPO_ROOT / "tests"))

from umi_tools_compat.conftest import _get_tests_dir, substitute_placeholders  # noqa: E402


def generate(name, values, input_dir):
    tmpdir = tempfile.mkdtemp()
    stdin_flag = f"--stdin={input_dir / values['stdin']}" if values.get("stdin") else ""
    opts = substitute_placeholders(values["options"], input_dir, tmpdir)
    stdout_path = os.path.join(tmpdir, "stdout")
    statement = (
        f"/bin/bash -c 'umi_tools {opts} {stdin_flag} -L {tmpdir}/umi_tools.log > {stdout_path}'"
    )
    result = subprocess.run(statement, shell=True, capture_output=True, cwd=tmpdir)
    if result.returncode != 0:
        sys.exit(f"{name}: umi_tools failed\n{statement}\n{result.stderr.decode(errors='replace')}")

    for output_name, ref_name in zip(values["outputs"], values["references"]):
        source = stdout_path if output_name == "stdout" else os.path.join(tmpdir, output_name)
        shutil.copyfile(source, HERE / ref_name)
        print(f"{name}: wrote {ref_name}")
    shutil.rmtree(tmpdir)


def main(selected):
    with open(HERE / "cases.yaml") as f:
        cases = yaml.safe_load(f) or {}
    input_dir = _get_tests_dir()
    if not input_dir.exists():
        sys.exit(f"upstream tests directory not found at {input_dir}")
    unknown = set(selected) - set(cases)
    if unknown:
        sys.exit(f"unknown cases: {sorted(unknown)}")
    for name, values in cases.items():
        if not selected or name in selected:
            generate(name, values, input_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
