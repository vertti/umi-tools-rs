use std::fs::File;

use flate2::read::MultiGzDecoder;
use umi_core::extract::{ExtractConfig, extract_reads};
use umi_core::pattern::{BarcodePattern, PrimeEnd};

fn read_gzipped_fastq(path: &str) -> String {
    let file = File::open(path).unwrap();
    let mut decoder = MultiGzDecoder::new(file);
    let mut content = String::new();
    std::io::Read::read_to_string(&mut decoder, &mut content).unwrap();
    content
}

#[test]
fn extract_matches_umi_tools_reference_output() {
    let pattern = BarcodePattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
    let config = ExtractConfig {
        pattern,
        umi_separator: b'_',
    };

    let input = File::open("tests/data/slim.fastq.gz").unwrap();
    let reader = MultiGzDecoder::new(input);
    let mut output = Vec::new();

    let stats = extract_reads(&config, reader, &mut output).unwrap();

    assert!(stats.input_reads > 0, "should have read some records");
    assert_eq!(
        stats.input_reads, stats.output_reads,
        "no reads should be filtered"
    );
    assert_eq!(stats.too_short, 0);

    let actual = String::from_utf8(output).unwrap();
    let expected = read_gzipped_fastq("tests/data/slim_extracted_NNNXXXXNN.fastq.gz");

    assert_eq!(actual, expected, "output should match umi-tools reference");
}
