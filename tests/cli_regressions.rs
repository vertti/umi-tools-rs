use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run_unchecked(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_umi-tools-rs"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn run(args: &[&str], input: &[u8]) -> Output {
    let output = run_unchecked(args, input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

#[test]
fn extraction_preserves_the_complete_separator() {
    let output = run(
        &["extract", "--bc-pattern=CCNN", "--umi-separator=::_"],
        b"@r comment\nCCAATT\n+\nIIIIII\n",
    );
    assert_eq!(output.stdout, b"@r::_CC::_AA comment\nTT\n+\nII\n");
}

#[test]
fn count_tab_preserves_cell_names_with_multibyte_separators() {
    let output = run(
        &["count_tab", "--per-cell", "--barcode-separator=::_"],
        b"r::_CELL1::_AAAA\tGENE1\nr::_CELL2::_AAAA\tGENE1\n",
    );
    assert_eq!(
        output.stdout,
        b"cell\tgene\tcount\nCELL1\tGENE1\t1\nCELL2\tGENE1\t1\n"
    );
}

const MGI_READ1: &[u8] = b"@DS100001332L1C002R00300000244/1\nACGTGCCATGCTAGC\n+\nIIIIIII;;;;;;;;\n";
const MGI_READ2: &[u8] = b"@DS100001332L1C002R00300000244/2\nCGAATCCTTAACCGG\n+\nIIIIIII;;;;;;;;\n";

fn mgi_extract_args(test: &str, extra: &[&str]) -> (Vec<String>, String) {
    let dir = env!("CARGO_TARGET_TMPDIR");
    let read2_in = format!("{dir}/{test}_R2.fastq");
    let read2_out = format!("{dir}/{test}_R2_out.fastq");
    std::fs::write(&read2_in, MGI_READ2).unwrap();
    let _ = std::fs::remove_file(&read2_out);
    let mut args = vec![
        "extract".to_string(),
        "--bc-pattern=NNNNNXX".to_string(),
        "--bc-pattern2=NNNNNXX".to_string(),
        format!("--read2-in={read2_in}"),
        format!("--read2-out={read2_out}"),
    ];
    args.extend(extra.iter().map(ToString::to_string));
    (args, read2_out)
}

#[test]
fn extraction_rejects_mgi_read_pair_suffixes_by_default() {
    let (args, _) = mgi_extract_args("mgi_default", &[]);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let output = run_unchecked(&args, MGI_READ1);
    assert!(!output.status.success());
    assert_eq!(output.stdout, b"");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "read pairs do not match: DS100001332L1C002R00300000244/1 != DS100001332L1C002R00300000244/2"
        ),
        "{stderr}"
    );
}

#[test]
fn extraction_gives_mgi_mates_one_name_when_ignoring_read_pair_suffixes() {
    let (args, read2_out) = mgi_extract_args("mgi_ignore", &["--ignore-read-pair-suffixes"]);
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let output = run(&args, MGI_READ1);
    assert_eq!(
        output.stdout,
        b"@DS100001332L1C002R00300000244_ACGTGCGAAT\nCCATGCTAGC\n+\nII;;;;;;;;\n"
    );
    assert_eq!(
        std::fs::read(read2_out).unwrap(),
        b"@DS100001332L1C002R00300000244_ACGTGCGAAT\nCCTTAACCGG\n+\nII;;;;;;;;\n"
    );
}
