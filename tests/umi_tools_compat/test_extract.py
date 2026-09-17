"""Compatibility tests: run umi-tools-rs against the real umi-tools test suite."""

import gzip
import subprocess

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


def run_paired_extract(binary, tmp_path, options, reads):
    inputs = [tmp_path / f"input{mate}.fq.gz" for mate in (1, 2)]
    outputs = [tmp_path / f"output{mate}.fq.gz" for mate in (1, 2)]
    for path, read in zip(inputs, reads):
        with gzip.open(path, "wt") as stream:
            stream.write(read)
    result = subprocess.run(
        [str(binary), "extract", "-v", "0", *options,
         "-I", str(inputs[0]), "--read2-in", str(inputs[1]),
         "--stdout", str(outputs[0]), "--read2-out", str(outputs[1])],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    actual = []
    for path in outputs:
        with gzip.open(path, "rt") as stream:
            actual.append(stream.read())
    return actual


@pytest.mark.parametrize("method", ["string", "regex"])
@pytest.mark.parametrize("read2_only", [False, True])
@pytest.mark.parametrize("strip_suffixes", [False, True])
@pytest.mark.parametrize("three_prime", [False, True])
def test_paired_barcodes(
    rust_binary, tmp_path, method, read2_only, strip_suffixes, three_prime
):
    """Issue #64: concatenate both UMIs, trim both mates, and honor suffixes."""
    if method == "string":
        pattern = "NNNNNN"
    else:
        pattern = r".*(?P<umi_1>.{6})$" if three_prime else r"^(?P<umi_1>.{6})"
    options = [f"--extract-method={method}", f"--bc-pattern2={pattern}"]
    if read2_only:
        options.append("--read2-only")
    else:
        options.append(f"--bc-pattern={pattern}")
    if strip_suffixes:
        options.append("--ignore-read-pair-suffixes")
    if three_prime:
        options.append("--3prime")

    reads, expected = [], []
    umi = "AGCGTG" if read2_only else "ACCACGAGCGTG"
    for mate, barcode, sequence in [(1, "ACCACG", "ATTTCT"), (2, "AGCGTG", "TTTCTG")]:
        full_seq = sequence + barcode if three_prime else barcode + sequence
        quality = "ABCDEF" + "GHIJKL"
        # Without the flag, upstream requires identical names (including suffix).
        input_suffix = f"/{mate}" if strip_suffixes else "/1"
        reads.append(f"@name{input_suffix} comment{mate}\n{full_seq}\n+\n{quality}\n")
        if mate == 1 and read2_only:
            sequence, trimmed_quality = full_seq, quality
        else:
            trimmed_quality = quality[:6] if three_prime else quality[6:]
        suffix = "" if strip_suffixes else input_suffix
        expected.append(
            f"@name{suffix}_{umi} comment{mate}\n{sequence}\n+\n{trimmed_quality}\n"
        )
    assert run_paired_extract(rust_binary, tmp_path, options, reads) == expected


@pytest.mark.parametrize("low_mate", [1, 2])
@pytest.mark.parametrize("quality_option", ["threshold", "mask"])
def test_paired_umi_quality(rust_binary, tmp_path, low_mate, quality_option):
    options = ["--bc-pattern=NNNNNN", "--bc-pattern2=NNNNNN",
               "--ignore-read-pair-suffixes", "--quality-encoding=phred33",
               f"--quality-filter-{quality_option}=20"]
    reads = []
    for mate, sequence in [(1, "ACCACGATTTCT"), (2, "AGCGTGTTTCTG")]:
        quality = "!" + "I" * 11 if mate == low_mate else "I" * 12
        reads.append(f"@name/{mate}\n{sequence}\n+\n{quality}\n")
    umi = "NCCACGAGCGTG" if low_mate == 1 else "ACCACGNGCGTG"
    expected = ["", ""] if quality_option == "threshold" else [
        f"@name_{umi}\n{sequence}\n+\nIIIIII\n" for sequence in ("ATTTCT", "TTTCTG")
    ]
    assert run_paired_extract(rust_binary, tmp_path, options, reads) == expected
