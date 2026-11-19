"""Whitelists that umi_tools cannot build, checked against umi_tools master."""

import subprocess

import pytest


def reads(cells):
    return "".join(f"@r{i}\n{cell}ACGTAC\n+\nIIIIIIIIII\n" for i, cell in enumerate(cells))


INPUTS = {
    "empty": "",
    "one_read": reads(["AAAA"]),
    "one_barcode": reads(["AAAA"] * 5),
    "two_barcodes": reads(["AAAA", "AAAA", "CCCC"]),
    "four_tied": reads([cell for cell in ["AAAA", "CCCC", "GGGG", "TTTT"] for _ in range(3)]),
    "three_barcodes": reads(["AAAA"] * 3 + ["CCCC"] * 2 + ["GGGG"]),
}
DENSITY = ["--knee-method=density"]
LENIENT = [*DENSITY, "--allow-threshold-error"]

FAILS = [
    ("empty", []), ("empty", DENSITY), ("empty", LENIENT), ("empty", ["--set-cell-number=2"]),
    ("one_read", []), ("one_read", DENSITY), ("one_read", LENIENT),
    ("one_barcode", []), ("one_barcode", DENSITY), ("one_barcode", LENIENT),
    ("two_barcodes", DENSITY), ("two_barcodes", ["--set-cell-number=2"]),
    ("four_tied", DENSITY), ("four_tied", LENIENT), ("four_tied", ["--set-cell-number=2"]),
]

SUCCEEDS = [
    ("two_barcodes", [], "AAAA\t\t2\t\nCCCC\t\t1\t\n"),
    ("two_barcodes", LENIENT, ""),
    ("four_tied", [], "AAAA\t\t3\t\nCCCC\t\t3\t\nGGGG\t\t3\t\nTTTT\t\t3\t\n"),
    ("three_barcodes", [], "AAAA\t\t3\t\nCCCC\t\t2\t\n"),
    ("three_barcodes", DENSITY, "AAAA\t\t3\t\nCCCC\t\t2\t\n"),
    ("three_barcodes", ["--set-cell-number=2"], "AAAA\t\t3\t\nCCCC\t\t2\t\n"),
]


def whitelist(binary, tmp_path, name, options):
    source = tmp_path / "reads.fq"
    source.write_text(INPUTS[name])
    return subprocess.run(
        [str(binary), "whitelist", "--bc-pattern=CCCCNN", "-I", str(source), *options],
        capture_output=True, text=True,
    )


@pytest.mark.parametrize("name,options", FAILS)
def test_whitelist_fails(rust_binary, tmp_path, name, options):
    result = whitelist(rust_binary, tmp_path, name, options)
    assert result.returncode == 1, result.stdout
    assert "panicked" not in result.stderr


@pytest.mark.parametrize("name,options,expected", SUCCEEDS)
def test_whitelist_succeeds(rust_binary, tmp_path, name, options, expected):
    result = whitelist(rust_binary, tmp_path, name, options)
    assert result.returncode == 0, result.stderr
    assert result.stdout == expected
