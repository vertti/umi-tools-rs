"""Gzip streams preserve content and report input and output failures."""

import gzip
import resource
import signal
import subprocess

import pytest

from .test_alignment_regressions import write_bam


@pytest.mark.parametrize("split_members", [False, True])
def test_extract_reads_gzip_members_including_empty_members(rust_binary, tmp_path, split_members):
    data = b"@r1\nAATT\n+\nIIII\n@r2\nCCGG\n+\nJJJJ\n"
    parts = [data[:7], b"", data[7:]] if split_members else [data]
    source = tmp_path / "input.fastq.gz"
    source.write_bytes(b"".join(gzip.compress(part) for part in parts))
    result = subprocess.run(
        [str(rust_binary), "extract", "--bc-pattern=NN", "--stdin", str(source)],
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == b"@r1_AA\nTT\n+\nII\n@r2_CC\nGG\n+\nJJ\n"


@pytest.mark.parametrize("damage", ["crc", "size", "truncated"])
@pytest.mark.parametrize("after_valid_member", [False, True])
def test_extract_rejects_damaged_gzip_members(rust_binary, tmp_path, damage, after_valid_member):
    member = bytearray(gzip.compress(b"@r\nAATT\n+\nIIII\n"))
    if damage == "truncated":
        del member[-4:]
    else:
        member[-8 if damage == "crc" else -4] ^= 1
    prefix = gzip.compress(b"@first\nCCGG\n+\nIIII\n") if after_valid_member else b""
    source = tmp_path / "input.fastq.gz"
    source.write_bytes(prefix + member)
    result = subprocess.run(
        [str(rust_binary), "extract", "--bc-pattern=NN", "--stdin", str(source)],
        capture_output=True,
    )
    assert result.returncode != 0, "damaged gzip input must not be reported as success"
    assert b"Error:" in result.stderr


@pytest.mark.parametrize("level", [1, 3, 9])
def test_gzip_output_decodes_with_python_at_each_level(rust_binary, tmp_path, level):
    output = tmp_path / "output.fastq.gz"
    result = subprocess.run(
        [str(rust_binary), "extract", "--bc-pattern=NN", "--compresslevel", str(level),
         "--stdout", str(output)],
        input=b"@r\nAATT\n+\nIIII\n", capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert gzip.decompress(output.read_bytes()) == b"@r_AA\nTT\n+\nII\n"


@pytest.mark.parametrize("command,options,data", [
    ("count_tab", [], "r_AAAA\tGENE1\n"),
    ("extract", ["--bc-pattern=NN"], "@r\nAATT\n+\nIIII\n"),
    ("whitelist", ["--bc-pattern=CC", "--set-cell-number=1"],
     "@r1\nAATT\n+\nIIII\n@r2\nAATT\n+\nIIII\n@r3\nCCTT\n+\nIIII\n"),
])
def test_gzip_finalization_failure_is_reported(rust_binary, tmp_path, command, options, data):
    path = tmp_path / "out.gz"
    args = [str(rust_binary), command, *options, "--stdout", str(path)]
    complete = subprocess.run(args, input=data, capture_output=True, text=True)
    assert complete.returncode == 0, complete.stderr
    assert gzip.decompress(path.read_bytes())
    limit = path.stat().st_size - 1

    def limit_output():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (limit, limit))

    failed = subprocess.run(args, input=data, capture_output=True, text=True,
                            preexec_fn=limit_output)
    assert failed.returncode != 0, "a truncated gzip file must not be reported as success"
    assert "Error:" in failed.stderr


@pytest.mark.parametrize("flag,expected", [
    ("--whitelist", "@r_AA_\nTT\n+\nII\n"),
    ("--blacklist", ""),
])
def test_extract_reads_compressed_barcode_lists(rust_binary, tmp_path, flag, expected):
    path = tmp_path / "barcodes.tsv.gz"
    path.write_bytes(gzip.compress(b"AA\n"))
    result = subprocess.run(
        [str(rust_binary), "extract", "--bc-pattern=CC", flag, str(path)],
        input="@r\nAATT\n+\nIIII\n", capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == expected


@pytest.mark.parametrize("paired", [False, True])
def test_dedup_reads_compressed_umi_lists(rust_binary, tmp_path, paired):
    source, output, whitelist = [tmp_path / name for name in ("in.bam", "out.bam", "umis.gz")]
    write_bam(source, [])
    whitelist.write_bytes(gzip.compress(b"AA\n" if paired else b"AAAA\n"))
    args = [str(rust_binary), "dedup", "--stdin", str(source), "--stdout", str(output),
            "--filter-umi", "--umi-whitelist", str(whitelist)]
    if paired:
        args += ["--umi-whitelist-paired", str(whitelist)]
    result = subprocess.run(args, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    import pysam
    with pysam.AlignmentFile(output) as stream:
        assert [read.query_name for read in stream] == ["r_AAAA"]
