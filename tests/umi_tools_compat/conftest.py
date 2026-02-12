import os
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
BINARY_PATH = REPO_ROOT / "target" / "debug" / "umi-tools-rs"

DEFAULT_UMI_TOOLS_TESTS_DIR = REPO_ROOT.parent / "umi-tools" / "tests"


def _get_tests_dir():
    return Path(
        os.environ.get("UMI_TOOLS_TESTS_DIR", str(DEFAULT_UMI_TOOLS_TESTS_DIR))
    )


def _load_tests_yaml():
    tests_dir = _get_tests_dir()
    yaml_path = tests_dir / "tests.yaml"
    if not yaml_path.exists():
        return {}
    with open(yaml_path) as f:
        return yaml.safe_load(f)


@pytest.fixture(scope="session")
def rust_binary():
    if not BINARY_PATH.exists():
        pytest.fail(
            f"Rust binary not found at {BINARY_PATH}. Run `cargo build` first."
        )
    return BINARY_PATH


@pytest.fixture(scope="session")
def umi_tools_tests_dir():
    path = _get_tests_dir()
    if not path.exists():
        pytest.skip(f"umi-tools tests directory not found at {path}")
    return path
