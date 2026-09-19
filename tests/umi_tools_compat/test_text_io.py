"""Compressed output must report failures, including at finalization."""

import gzip
import resource
import signal
import subprocess

import pytest


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
