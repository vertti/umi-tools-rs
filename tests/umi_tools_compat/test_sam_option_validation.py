"""Option combinations umi_tools rejects for BAM input (validateSamOptions).

Each case exits 1 or 0 in umi_tools master as asserted. Several upstream checks
look for the flag in the command line as typed, so an abbreviated flag passes.
"""

import subprocess

import pytest

COMMANDS = {
    "dedup": [],
    "group": ["--group-out=groups.tsv"],
    "count": ["--per-contig", "--per-gene"],
}

TAGS = "--umi-tag and/or --cell-tag options provided. Need to set --extract-umi-method=tag"
PAIRED_ONLY = "is only compatible with paired end reads (--paired)"

REJECTED = [
    (["--umi-tag=RX"], TAGS),
    (["--cell-tag=CB"], TAGS),
    (["--umi-tag-split=-"], TAGS),
    (["--unmapped-reads=use"], f"--unmapped-reads=use {PAIRED_ONLY}"),
    (["--chimeric-pairs=discard"], f"--chimeric-pairs {PAIRED_ONLY}"),
    (["--unpaired-reads=discard"], f"--unpaired-reads {PAIRED_ONLY}"),
    (["--ignore-tlen"], f"--ignore-tlen {PAIRED_ONLY}"),
]

ACCEPTED = [
    ["--chimeric=discard"],
    ["--unpaired=discard"],
    ["--paired", "--chimeric-pairs=discard", "--unpaired-reads=discard", "--ignore-tlen",
     "--unmapped-reads=use"],
]


def run(binary, tests_dir, tmp_path, command, options):
    return subprocess.run(
        [str(binary), command, *COMMANDS[command], *options, "-I", str(tests_dir / "chr19.bam")],
        capture_output=True, cwd=tmp_path,
    )


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("options,message", REJECTED)
def test_rejected(rust_binary, umi_tools_tests_dir, tmp_path, command, options, message):
    result = run(rust_binary, umi_tools_tests_dir, tmp_path, command, options)
    assert result.returncode == 1, result.stderr
    assert message in result.stderr.decode()


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("options", ACCEPTED)
def test_accepted(rust_binary, umi_tools_tests_dir, tmp_path, command, options):
    result = run(rust_binary, umi_tools_tests_dir, tmp_path, command, options)
    assert result.returncode == 0, result.stderr


def test_dedup_statistics_need_umis(rust_binary, umi_tools_tests_dir, tmp_path):
    result = run(rust_binary, umi_tools_tests_dir, tmp_path, "dedup",
                 ["--ignore-umi", "--output-stats=stats"])
    assert result.returncode == 1, result.stderr
    assert ("'--output-stats' and '--ignore-umi' options cannot be used together"
            in result.stderr.decode())
