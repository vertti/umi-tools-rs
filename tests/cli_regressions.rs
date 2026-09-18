use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_umi-tools-rs"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    let output = child.wait_with_output().unwrap();
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
