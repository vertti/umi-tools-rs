use std::fs::File;

use flate2::read::MultiGzDecoder;
use umi_core::extract::{ExtractConfig, QualityEncoding, extract_reads};
use umi_core::pattern::{BarcodePattern, PrimeEnd, RegexPattern, StringPattern};

fn read_gzipped_fastq(path: &str) -> String {
    let file = File::open(path).unwrap();
    let mut decoder = MultiGzDecoder::new(file);
    let mut content = String::new();
    std::io::Read::read_to_string(&mut decoder, &mut content).unwrap();
    content
}

#[test]
fn extract_string_matches_umi_tools_reference() {
    let pattern =
        BarcodePattern::String(StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap());
    let config = ExtractConfig {
        pattern: Some(pattern),
        pattern2: None,
        umi_separator: b'_',
        quality_filter_threshold: None,
        quality_encoding: QualityEncoding::default(),
        whitelist: None,
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

#[test]
fn extract_regex_matches_string_method() {
    // Per umi-tools tests.yaml: regex ^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2}) == string NNNXXXXNN
    let string_pattern =
        BarcodePattern::String(StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap());
    let regex_pattern =
        BarcodePattern::Regex(RegexPattern::parse(r"^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})").unwrap());

    let string_config = ExtractConfig {
        pattern: Some(string_pattern),
        pattern2: None,
        umi_separator: b'_',
        quality_filter_threshold: None,
        quality_encoding: QualityEncoding::default(),
        whitelist: None,
    };
    let regex_config = ExtractConfig {
        pattern: Some(regex_pattern),
        pattern2: None,
        umi_separator: b'_',
        quality_filter_threshold: None,
        quality_encoding: QualityEncoding::default(),
        whitelist: None,
    };

    let string_input = File::open("tests/data/slim.fastq.gz").unwrap();
    let regex_input = File::open("tests/data/slim.fastq.gz").unwrap();

    let mut string_output = Vec::new();
    let mut regex_output = Vec::new();

    let string_stats = extract_reads(
        &string_config,
        MultiGzDecoder::new(string_input),
        &mut string_output,
    )
    .unwrap();
    let regex_stats = extract_reads(
        &regex_config,
        MultiGzDecoder::new(regex_input),
        &mut regex_output,
    )
    .unwrap();

    assert_eq!(string_stats.input_reads, regex_stats.input_reads);
    assert_eq!(string_stats.output_reads, regex_stats.output_reads);
    assert_eq!(regex_stats.no_match, 0);

    let string_out = String::from_utf8(string_output).unwrap();
    let regex_out = String::from_utf8(regex_output).unwrap();

    assert_eq!(
        string_out, regex_out,
        "regex and string methods should produce identical output"
    );
}
