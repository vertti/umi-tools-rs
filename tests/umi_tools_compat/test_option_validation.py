"""Option combinations umi_tools rejects before reading any input.

Each case exits 1 in umi_tools master with the quoted message.
"""

import subprocess

import pytest

READ1 = "@a\nACGTACGTAC\n+\nIIIIIIIIII\n@b\nTTGGCCAATT\n+\nIIIIIIIIII\n"
READ2 = "@a\nGGGGACGTAC\n+\nIIIIIIIIII\n@b\nCCGGCCAATT\n+\nIIIIIIIIII\n"

REJECTED = [
    ("extract", ["--bc-pattern=CCCC"], "do not include any umi bases"),
    ("extract", ["--bc-pattern=XXXX"], "do not include any umi bases"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<cell_1>.{4})"], "do not include any umi groups"),
    ("extract", ["--bc-pattern=NNNN", "--whitelist=wl.tsv"], "do not include any cell bases"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})", "--whitelist=wl.tsv"], "do not include any cell groups"),
    ("extract", ["--bc-pattern=NNNN", "--filtered-out=f.fq"], "Reads will not be filtered"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})", "--read2-in=r2.fq", "--filtered-out=f.fq"], "Must supply both --filtered-out and --filtered-out2"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})", "--read2-in=r2.fq", "--filtered-out2=f.fq"], "Must supply both --filtered-out and --filtered-out2"),
    ("extract", ["--bc-pattern=NNNN", "--quality-filter-threshold=10"], "must provide a quality encoding"),
    ("extract", ["--bc-pattern=NNNN", "--quality-filter-mask=10"], "must provide a quality encoding"),
    ("extract", ["--bc-pattern=NNNN", "--bc-pattern2=NNNN", "--read2-in=r2.fq", "--either-read"], "requires --extract-method=regex"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})(?P<cell_1>.{4})", "--bc-pattern2=(?P<umi_1>.{4})", "--read2-in=r2.fq", "--either-read"], "not currently compatible with cell barcode extraction"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})", "--read2-in=r2.fq", "--either-read"], "requires --bc-pattern=[PATTERN1]"),
    ("extract", ["--bc-pattern2=NNNN"], "must specify a paired fastq"),
    ("extract", ["--read2-in=r2.fq"], "Must supply --bc-pattern and/or --bc-pattern2"),
    ("extract", [], "Must supply --bc-pattern for single-end"),
    ("extract", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})", "--filtered-out=f.fq", "--filtered-out2=g.fq"], "Cannot use --filtered-out2 without read2 input"),
    ("whitelist", ["--bc-pattern=NNNN"], "do not include any cell bases"),
    ("whitelist", ["--bc-pattern=CCCC"], "do not include any umi bases"),
    ("whitelist", ["--extract-method=regex", "--bc-pattern=(?P<umi_1>.{4})"], "do not include any cell groups"),
    ("whitelist", ["--bc-pattern=CCCCNN", "--expect-cells=5"], "Cannot use --expect-cells with 'distance' knee"),
    ("whitelist", ["--bc-pattern=CCCCNN", "--knee-method=density", "--expect-cells=5", "--set-cell-number=5"], "Cannot supply both --expect-cells and"),
    ("whitelist", ["--bc-pattern=CCCCNN", "--filtered-out=f.fq"], "Reads will not be filtered"),
    ("whitelist", ["--bc-pattern=CCCCNN", "--bc-pattern2=CC"], "must specify a paired fastq"),
]

ACCEPTED = [
    ("extract", ["--bc-pattern=NNNN", "--quality-filter-threshold=0"]),
    ("extract", ["--bc-pattern=NNNN", "--quality-filter-mask=0"]),
    ("extract", ["--bc-pattern=NNNN", "--quality-filter-threshold=10", "--quality-encoding=phred33"]),
    ("extract", ["--bc-pattern=CCNN", "--whitelist=wl.tsv", "--filtered-out=f.fq"]),
]


def run(binary, tmp_path, command, options):
    (tmp_path / "r1.fq").write_text(READ1)
    (tmp_path / "r2.fq").write_text(READ2)
    (tmp_path / "wl.tsv").write_text("ACGT\n")
    return subprocess.run(
        [str(binary), command, *options, "-I", "r1.fq"],
        capture_output=True, text=True, cwd=tmp_path,
    )


@pytest.mark.parametrize("command,options,message", REJECTED)
def test_rejected_option_combinations(rust_binary, tmp_path, command, options, message):
    result = run(rust_binary, tmp_path, command, options)
    assert result.returncode == 1, result.stderr
    assert message in result.stderr


@pytest.mark.parametrize("command,options", ACCEPTED)
def test_accepted_option_combinations(rust_binary, tmp_path, command, options):
    result = run(rust_binary, tmp_path, command, options)
    assert result.returncode == 0, result.stderr
