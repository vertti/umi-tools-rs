"""Output routing must not change extraction, filtering or pair matching.

Expectations checked against Python umi_tools 1.1.6. These synthetic cases do
not depend on the upstream fixture checkout.
"""

import gzip
import subprocess

import pytest

PATTERN = r"^(?P<cell_1>.{2})(?P<umi_1>.{2})(?P<discard_1>GG)"
MODES = ["single", "r1", "r2", "both"]
ROUTES = [(mode, stdout2) for mode in MODES for stdout2 in (False, True)
          if mode != "single" or not stdout2]


def fastq(name, sequence, mate=1, quality=None):
    quality = quality or "I" * len(sequence)
    return f"@{name}/{mate} comment{mate}\n{sequence}\n+\n{quality}\n"


def run_extract(binary, tmp_path, options, reads, *, stdout2=False, emit_read2=True):
    inputs = [tmp_path / f"input{mate}.fq.gz" for mate in (1, 2)]
    outputs = [tmp_path / name for name in
               ("output1.fq.gz", "output2.fq.gz", "filtered1.fq.gz", "filtered2.fq.gz")]
    args = [str(binary), "extract", "-v", "0", "--ignore-read-pair-suffixes", *options,
            "-I", str(inputs[0]), "--stdout", str(outputs[0]),
            "--filtered-out", str(outputs[2])]
    if len(reads) == 2:
        args += ["--read2-in", str(inputs[1]), "--filtered-out2", str(outputs[3])]
        if emit_read2:
            args += ["--read2-out", str(outputs[1])]
        if stdout2:
            args += ["--read2-stdout"]
    for path, content in zip(inputs, reads):
        with gzip.open(path, "wt") as stream:
            stream.write(content)
    result = subprocess.run(args, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    actual = []
    for path in outputs:
        if path.exists():
            with gzip.open(path, "rt") as stream:
                actual.append(stream.read())
        else:
            actual.append("")
    return actual


def mode_options(mode):
    options = ["--extract-method=regex"]
    # Python forbids cell extraction in either-read mode.
    pattern = PATTERN.replace("cell_1", "discard_2") if mode.startswith("either") else PATTERN
    if mode != "r2":
        options += [f"--bc-pattern={pattern}"]
    if mode in ("r2", "both") or mode.startswith("either"):
        options += [f"--bc-pattern2={pattern}"]
    if mode.startswith("either"):
        options += ["--either-read", "--either-read-resolve=quality", "--quality-encoding=phred33"]
    return options


def pair_sequences(mode, cell):
    active1 = mode not in ("r2", "either-r2")
    active2 = mode in ("r2", "both", "either-r2", "either-both")
    cell2 = cell[2:] if mode == "both" else cell
    return [cell[:2] + "ACGGTATA" if active1 else "TTTTTTTTTT",
            cell2 + "GTGGCGCG" if active2 else "TTTTTTTTTT"]


@pytest.mark.parametrize("mode,stdout2", ROUTES)
def test_filters_independent_of_output(rust_binary, tmp_path, mode, stdout2):
    combined = mode == "both"
    good = "AATT" if combined else "AA"
    corrected = "ACTT" if combined else "AC"
    bad = "GGTT" if combined else "GG"
    blocked = "CCTT" if combined else "CC"
    whitelist = tmp_path / "whitelist.tsv"
    whitelist.write_text(f"{good}\t{corrected}\n{blocked}\t{blocked}\n")
    blacklist = tmp_path / "blacklist.tsv"
    blacklist.write_text(blocked + "\n")
    options = mode_options(mode) + [f"--whitelist={whitelist}",
                                  "--error-correct-cell", f"--blacklist={blacklist}"]
    reads = ["", ""]
    accepted, rejected = ["", ""], ["", ""]
    for name, cell in [("good", good), ("corrected", corrected), ("bad", bad), ("blocked", blocked)]:
        sequences = pair_sequences(mode, cell)
        for index, sequence in enumerate(sequences):
            reads[index] += fastq(name, sequence, index + 1)
            if name in ("bad", "blocked"):
                rejected[index] += fastq(name, sequence, index + 1).replace(f"/{index + 1} ", " ")
                continue
            trimmed = sequence[6:] if "GG" == sequence[4:6] else sequence
            umi = "ACGT" if combined else "GT" if mode == "r2" else "AC"
            header_mate = index + 1
            accepted[index] += f"@{name}_{good}_{umi} comment{header_mate}\n{trimmed}\n+\n{'I' * len(trimmed)}\n"
    if mode == "single":
        reads = reads[:1]
        accepted[1] = rejected[1] = ""
    expected = ([accepted[1], ""] if stdout2 else accepted) + rejected
    assert run_extract(rust_binary, tmp_path, options, reads, stdout2=stdout2) == expected


@pytest.mark.parametrize("mode", ["r1", "r2", "both", "either-r1", "either-r2"])
@pytest.mark.parametrize("stdout2", [False, True])
def test_reconciliation_independent_of_output(rust_binary, tmp_path, mode, stdout2):
    sequences = pair_sequences(mode, "AATT" if mode == "both" else "AA")
    first = fastq("kept", sequences[0])
    second = fastq("extra", sequences[1], 2) + fastq("kept", sequences[1], 2)
    options = mode_options(mode) + ["--reconcile-pairs"]
    actual = run_extract(rust_binary, tmp_path, options, [first, second], stdout2=stdout2)
    baseline = run_extract(rust_binary, tmp_path, mode_options(mode),
                           [first, fastq("kept", sequences[1], 2)], stdout2=stdout2)
    assert actual == baseline
    assert "@kept_" in actual[0]
    assert "extra" not in "".join(actual)


@pytest.mark.parametrize("mode", MODES)
def test_corrected_cell_is_checked_against_blacklist(rust_binary, tmp_path, mode):
    good, bad = ("AATT", "ACTT") if mode == "both" else ("AA", "AC")
    whitelist = tmp_path / "whitelist.tsv"
    whitelist.write_text(f"{good}\t{bad}\n")
    blacklist = tmp_path / "blacklist.tsv"
    blacklist.write_text(good + "\n")
    sequences = pair_sequences(mode, bad)
    reads = [fastq("bad", seq, mate) for mate, seq in enumerate(sequences, 1)]
    if mode == "single":
        reads = reads[:1]
    options = mode_options(mode) + [f"--whitelist={whitelist}", "--error-correct-cell",
                                  f"--blacklist={blacklist}"]
    actual = run_extract(rust_binary, tmp_path, options, reads)
    assert actual[:2] == ["", ""]
    assert actual[2] == reads[0].replace("/1 ", " ")
    assert actual[3] == (reads[1].replace("/2 ", " ") if len(reads) == 2 else "")


@pytest.mark.parametrize("mode", ["either-r1", "either-r2", "either-both"])
@pytest.mark.parametrize("stdout2", [False, True])
@pytest.mark.parametrize("quality_option", ["threshold", "mask"])
def test_either_read_filters_and_headers(rust_binary, tmp_path, mode, stdout2, quality_option):
    options = mode_options(mode) + [f"--quality-filter-{quality_option}=20"]
    reads, accepted, rejected = ["", ""], ["", ""], ["", ""]
    for name in ("good", "low"):
        sequences = pair_sequences(mode, "AA")
        umi = "GT" if mode == "either-r2" else "AC"
        if name == "low" and quality_option == "mask":
            umi = "N" + umi[1:]
        for index, sequence in enumerate(sequences):
            matched = sequence[4:6] == "GG"
            quality = "II!IIIIIII" if matched and name == "low" else "I" * 10
            original = fastq(name, sequence, index + 1, quality)
            reads[index] += original
            if name == "low" and quality_option == "threshold":
                rejected[index] += original.replace(f"/{index + 1} ", " ")
                continue
            if matched and mode != "either-both":
                sequence, quality = sequence[6:], quality[6:]
            accepted[index] += f"@{name}_{umi} comment1\n{sequence}\n+\n{quality}\n"
    assert run_extract(rust_binary, tmp_path, options, reads, stdout2=stdout2) == (
        ([accepted[1], ""] if stdout2 else accepted) + rejected
    )


@pytest.mark.parametrize("mode", ["r1", "both", "either-r2"])
def test_read2_output_is_optional(rust_binary, tmp_path, mode):
    sequences = pair_sequences(mode, "AATT" if mode == "both" else "AA")
    reads = [fastq("kept", seq, mate) for mate, seq in enumerate(sequences, 1)]
    actual = run_extract(rust_binary, tmp_path, mode_options(mode), reads, emit_read2=False)
    expected = run_extract(rust_binary, tmp_path, mode_options(mode), reads)
    assert actual == [expected[0], "", expected[2], expected[3]]


@pytest.mark.parametrize("mode", ["r1", "r2", "both", "either-r1"])
def test_regex_rejection_writes_both_original_reads(rust_binary, tmp_path, mode):
    reads = [fastq("unmatched", "TTTTTTTTTT", mate) for mate in (1, 2)]
    actual = run_extract(rust_binary, tmp_path, mode_options(mode), reads)
    assert actual == ["", "", reads[0].replace("/1 ", " "), reads[1].replace("/2 ", " ")]


@pytest.mark.parametrize("mode", ["r1", "both", "either-r1"])
def test_mismatched_names_are_rejected(rust_binary, tmp_path, mode):
    sequences = pair_sequences(mode, "AATT" if mode == "both" else "AA")
    reads = [fastq("first", sequences[0]), fastq("different", sequences[1], 2)]
    # The helper asserts success; the command must instead reject this pair.
    with pytest.raises(AssertionError, match="(?i)read pairs do not match"):
        run_extract(rust_binary, tmp_path, mode_options(mode), reads)
