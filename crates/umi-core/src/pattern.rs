use std::borrow::Cow;
use std::ops::Range;

use crate::error::ExtractError;

/// Result of extracting barcodes from a single read's sequence and quality.
#[derive(Debug, Clone)]
pub struct ExtractionResult {
    pub umi: Vec<u8>,
    pub umi_quality: Vec<u8>,
    pub cell_barcode: Vec<u8>,
    pub trimmed_sequence: Vec<u8>,
    pub trimmed_quality: Vec<u8>,
}

/// Borrow unchanged slices while the input record remains in the parser buffer.
pub(crate) struct ExtractionView<'a> {
    pub umi: Cow<'a, [u8]>,
    pub umi_quality: Cow<'a, [u8]>,
    pub cell_barcode: Cow<'a, [u8]>,
    pub trimmed_sequence: Cow<'a, [u8]>,
    pub trimmed_quality: Cow<'a, [u8]>,
}

impl From<ExtractionResult> for ExtractionView<'_> {
    fn from(result: ExtractionResult) -> Self {
        Self {
            umi: Cow::Owned(result.umi),
            umi_quality: Cow::Owned(result.umi_quality),
            cell_barcode: Cow::Owned(result.cell_barcode),
            trimmed_sequence: Cow::Owned(result.trimmed_sequence),
            trimmed_quality: Cow::Owned(result.trimmed_quality),
        }
    }
}

impl ExtractionView<'_> {
    fn into_owned(self) -> ExtractionResult {
        ExtractionResult {
            umi: self.umi.into_owned(),
            umi_quality: self.umi_quality.into_owned(),
            cell_barcode: self.cell_barcode.into_owned(),
            trimmed_sequence: self.trimmed_sequence.into_owned(),
            trimmed_quality: self.trimmed_quality.into_owned(),
        }
    }
}

/// Which end of the read to extract the barcode from (string method only).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimeEnd {
    Five,
    Three,
}

/// A parsed barcode pattern that knows how to extract UMI/cell/sample bases from a read.
#[derive(Debug, Clone)]
pub enum BarcodePattern {
    String(StringPattern),
    Regex(RegexPattern),
}

impl BarcodePattern {
    pub(crate) fn extract_view<'a>(
        &self,
        sequence: &'a [u8],
        quality: &'a [u8],
    ) -> Result<ExtractionView<'a>, ExtractError> {
        match self {
            Self::String(pattern) => pattern.extract_view(sequence, quality),
            Self::Regex(pattern) => pattern.extract(sequence, quality).map(ExtractionView::from),
        }
    }

    /// # Errors
    /// Returns error if the read is too short (string method) or doesn't match (regex method).
    pub fn extract(
        &self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<ExtractionResult, ExtractError> {
        match self {
            Self::String(p) => p.extract(sequence, quality),
            Self::Regex(p) => p.extract(sequence, quality),
        }
    }
}

/// String-method pattern using fixed-position characters.
///
/// Pattern characters:
/// - `N` — UMI base (extracted to read name)
/// - `C` — Cell barcode base (extracted to read name)
/// - `X` — Sample/discard base (stays in output sequence, removed from barcode region)
#[derive(Debug, Clone)]
pub struct StringPattern {
    umi_positions: Vec<usize>,
    cell_positions: Vec<usize>,
    sample_positions: Vec<usize>,
    umi_range: Option<Range<usize>>,
    cell_range: Option<Range<usize>>,
    sample_range: Option<Range<usize>>,
    pattern_length: usize,
    prime_end: PrimeEnd,
}

impl StringPattern {
    /// Parse a string-method pattern like `NNNXXXXNN`.
    ///
    /// # Errors
    /// Returns error if pattern is empty or contains characters other than N, X, C.
    pub fn parse(pattern_str: &str, prime_end: PrimeEnd) -> Result<Self, ExtractError> {
        if pattern_str.is_empty() {
            return Err(ExtractError::InvalidPattern(
                "pattern must not be empty".into(),
            ));
        }

        let mut umi_positions = Vec::new();
        let mut cell_positions = Vec::new();
        let mut sample_positions = Vec::new();

        for (i, ch) in pattern_str.chars().enumerate() {
            match ch {
                'N' => umi_positions.push(i),
                'C' => cell_positions.push(i),
                'X' => sample_positions.push(i),
                other => {
                    return Err(ExtractError::InvalidPattern(format!(
                        "pattern contains invalid character '{other}' at position {i}; \
                         only N, X, C are allowed"
                    )));
                }
            }
        }

        let umi_range = as_contiguous_range(&umi_positions);
        let cell_range = as_contiguous_range(&cell_positions);
        let sample_range = as_contiguous_range(&sample_positions);

        Ok(Self {
            umi_positions,
            cell_positions,
            sample_positions,
            umi_range,
            cell_range,
            sample_range,
            pattern_length: pattern_str.len(),
            prime_end,
        })
    }

    /// Extract barcodes from a read's sequence and quality strings.
    ///
    /// # Errors
    /// Returns error if the read is shorter than the pattern.
    pub fn extract(
        &self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<ExtractionResult, ExtractError> {
        self.extract_view(sequence, quality)
            .map(ExtractionView::into_owned)
    }

    fn extract_view<'a>(
        &self,
        sequence: &'a [u8],
        quality: &'a [u8],
    ) -> Result<ExtractionView<'a>, ExtractError> {
        if sequence.len() < self.pattern_length {
            return Err(ExtractError::ReadTooShort {
                read_len: sequence.len(),
                pattern_len: self.pattern_length,
            });
        }

        let (barcode_region, remaining_seq, barcode_qual, remaining_qual) = match self.prime_end {
            PrimeEnd::Five => (
                &sequence[..self.pattern_length],
                &sequence[self.pattern_length..],
                &quality[..self.pattern_length],
                &quality[self.pattern_length..],
            ),
            PrimeEnd::Three => (
                &sequence[sequence.len() - self.pattern_length..],
                &sequence[..sequence.len() - self.pattern_length],
                &quality[quality.len() - self.pattern_length..],
                &quality[..quality.len() - self.pattern_length],
            ),
        };

        let umi = extract_slice(barcode_region, self.umi_range.as_ref(), &self.umi_positions);
        let umi_quality = extract_slice(barcode_qual, self.umi_range.as_ref(), &self.umi_positions);
        let cell_barcode = extract_slice(
            barcode_region,
            self.cell_range.as_ref(),
            &self.cell_positions,
        );

        let (trimmed_sequence, trimmed_quality) = if self.sample_positions.is_empty() {
            (Cow::Borrowed(remaining_seq), Cow::Borrowed(remaining_qual))
        } else {
            let sample_seq = extract_slice(
                barcode_region,
                self.sample_range.as_ref(),
                &self.sample_positions,
            );
            let sample_qual = extract_slice(
                barcode_qual,
                self.sample_range.as_ref(),
                &self.sample_positions,
            );
            match self.prime_end {
                PrimeEnd::Five => (
                    Cow::Owned(join_slices(&sample_seq, remaining_seq)),
                    Cow::Owned(join_slices(&sample_qual, remaining_qual)),
                ),
                PrimeEnd::Three => (
                    Cow::Owned(join_slices(remaining_seq, &sample_seq)),
                    Cow::Owned(join_slices(remaining_qual, &sample_qual)),
                ),
            }
        };

        Ok(ExtractionView {
            umi,
            umi_quality,
            cell_barcode,
            trimmed_sequence,
            trimmed_quality,
        })
    }
}

/// Regex-method pattern using named capture groups.
///
/// Groups starting with `umi_` are extracted as UMI, `cell_` as cell barcode,
/// `discard_` as bases to remove. Everything else is kept in the output sequence.
#[derive(Debug, Clone)]
pub struct RegexPattern {
    pattern: regex::Regex,
    relaxed: Option<RelaxedPattern>,
}

/// The pattern with each fuzzy group widened to any text of its length.
///
/// Resolving captures with the exact pattern, an alternation over every choice of
/// substituted positions, is slow. The widened match is the fuzzy match whenever its
/// groups are within their substitution limits. Otherwise the groups are pinned to the
/// only places they fit, and the exact pattern is left for reads where they fit in several.
#[derive(Debug, Clone)]
struct RelaxedPattern {
    pattern: regex::Regex,
    /// As `pattern`, with group N starting on the byte N + 1 instead of any character.
    pinned: regex::Regex,
    groups: Vec<WidenedGroup>,
}

#[derive(Debug, Clone)]
struct WidenedGroup {
    capture: String,
    accepted_bytes: Vec<[bool; 128]>,
    limit: usize,
}

enum Placements {
    Impossible,
    Unique(Vec<usize>),
    Ambiguous,
}

impl RelaxedPattern {
    const MAX_GROUPS: usize = 8;

    /// `None` without fuzzy groups, or when one of them does not take part in every match.
    fn parse(pattern_str: &str) -> Result<Option<Self>, ExtractError> {
        let mut groups = Vec::new();
        let mut pins = Vec::new();
        let mut all_required = true;
        let rewritten = rewrite_fuzzy(pattern_str, &mut |group, limit, required| {
            let (header, items) = fuzzy_group_items(group)?;
            if limit >= items.len() {
                return Ok(format!("({header}.{{{}}})", items.len()));
            }
            all_required &= required;
            let capture = format!("__fuzzy_{}", groups.len());
            let widened = format!("(?P<{capture}>.{{{}}})", items.len());
            pins.push((
                widened.clone(),
                format!(
                    "(?P<{capture}>\\x{:02x}.{{{}}})",
                    groups.len() + 1,
                    items.len() - 1
                ),
            ));
            groups.push(WidenedGroup {
                capture,
                accepted_bytes: items
                    .iter()
                    .map(|item| accepted_bytes(item))
                    .collect::<Result<_, _>>()?,
                limit,
            });
            Ok(format!("({header}{widened})"))
        })?;
        if groups.is_empty()
            || groups.len() > Self::MAX_GROUPS
            || !all_required
            || rewritten.alternates_at_top_level
        {
            return Ok(None);
        }
        let compile = |pattern: &str| {
            regex::Regex::new(pattern)
                .map_err(|e| ExtractError::InvalidPattern(format!("invalid regex: {e}")))
        };
        let pinned = pins
            .iter()
            .fold(rewritten.pattern.clone(), |pattern, (widened, pinned)| {
                pattern.replace(widened, pinned)
            });
        Ok(Some(Self {
            pattern: compile(&rewritten.pattern)?,
            pinned: compile(&pinned)?,
            groups,
        }))
    }

    fn within_limits(&self, captures: &regex::Captures) -> bool {
        self.groups.iter().all(|group| {
            captures
                .name(&group.capture)
                .is_some_and(|text| group.accepts(text.as_str().as_bytes()))
        })
    }

    fn placements(&self, sequence: &[u8]) -> Placements {
        let mut starts = Vec::with_capacity(self.groups.len());
        for group in &self.groups {
            let mut fits = sequence
                .windows(group.accepted_bytes.len())
                .enumerate()
                .filter(|(_, window)| group.accepts(window))
                .map(|(start, _)| start);
            match (fits.next(), fits.next()) {
                (None, _) => return Placements::Impossible,
                (Some(start), None) => starts.push(start),
                (Some(_), Some(_)) => return Placements::Ambiguous,
            }
        }
        Placements::Unique(starts)
    }

    /// `sequence` with the byte N + 1 where group N starts.
    fn pin(sequence: &str, starts: &[usize]) -> Option<String> {
        let mut pinned = sequence.as_bytes().to_vec();
        for (marker, &start) in (1_u8..).zip(starts) {
            pinned[start] = marker;
        }
        String::from_utf8(pinned).ok()
    }
}

impl WidenedGroup {
    fn accepts(&self, text: &[u8]) -> bool {
        text.len() == self.accepted_bytes.len()
            && text
                .iter()
                .zip(&self.accepted_bytes)
                .filter(|&(&byte, accepted)| !accepted.get(usize::from(byte)).is_some_and(|&ok| ok))
                .count()
                <= self.limit
    }
}

type GroupSpan<'p> = (&'p str, usize, usize);

/// Spans of the named groups in the match of `pattern` that starts `haystack`.
fn anchored_spans<'p>(pattern: &'p regex::Regex, haystack: &str) -> Option<Vec<GroupSpan<'p>>> {
    let captures = pattern
        .captures(haystack)
        .filter(|captures| captures.get(0).is_some_and(|matched| matched.start() == 0))?;
    Some(named_spans(pattern, &captures))
}

fn named_spans<'p>(pattern: &'p regex::Regex, captures: &regex::Captures) -> Vec<GroupSpan<'p>> {
    pattern
        .capture_names()
        .flatten()
        .filter_map(|name| {
            captures
                .name(name)
                .map(|matched| (name, matched.start(), matched.end()))
        })
        .collect()
}

fn accepted_bytes(item: &str) -> Result<[bool; 128], ExtractError> {
    let item = regex::Regex::new(&format!("^(?:{item})$"))
        .map_err(|e| ExtractError::InvalidPattern(format!("invalid regex: {e}")))?;
    let mut accepted = [false; 128];
    for (byte, slot) in (0_u8..).zip(&mut accepted) {
        *slot = item.is_match(char::from(byte).encode_utf8(&mut [0; 4]));
    }
    Ok(accepted)
}

impl RegexPattern {
    /// Parse a regex pattern string.
    ///
    /// # Errors
    /// Returns error if the regex is invalid or has no `umi_` or `cell_` groups.
    pub fn parse(pattern_str: &str) -> Result<Self, ExtractError> {
        let processed = preprocess_fuzzy(pattern_str)?;

        let pattern = regex::RegexBuilder::new(&processed)
            .size_limit(EXACT_PATTERN_SIZE_LIMIT)
            .build()
            .map_err(|e| ExtractError::InvalidPattern(format!("invalid regex: {e}")))?;

        let has_barcode_group = pattern
            .capture_names()
            .flatten()
            .any(|name| name.starts_with("umi_") || name.starts_with("cell_"));

        if !has_barcode_group {
            return Err(ExtractError::InvalidPattern(
                "regex must contain at least one named group starting with 'umi_' or 'cell_'"
                    .into(),
            ));
        }

        Ok(Self {
            pattern,
            relaxed: RelaxedPattern::parse(pattern_str)?,
        })
    }

    /// Spans of the named groups in the match that starts `sequence`.
    fn group_spans(&self, sequence: &str) -> Option<Vec<GroupSpan<'_>>> {
        let Some(relaxed) = &self.relaxed else {
            return anchored_spans(&self.pattern, sequence);
        };
        let captures = relaxed
            .pattern
            .captures(sequence)
            .filter(|captures| captures.get(0).is_some_and(|matched| matched.start() == 0))?;
        if relaxed.within_limits(&captures) {
            return Some(named_spans(&relaxed.pattern, &captures));
        }
        match relaxed.placements(sequence.as_bytes()) {
            Placements::Impossible => None,
            Placements::Unique(starts) => RelaxedPattern::pin(sequence, &starts).map_or_else(
                || anchored_spans(&self.pattern, sequence),
                |pinned| anchored_spans(&relaxed.pinned, &pinned),
            ),
            Placements::Ambiguous => anchored_spans(&self.pattern, sequence),
        }
    }

    /// Extract barcodes from a read's sequence and quality strings.
    ///
    /// # Errors
    /// Returns `RegexNoMatch` if the regex doesn't match the sequence.
    pub fn extract(
        &self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<ExtractionResult, ExtractError> {
        let seq_str = std::str::from_utf8(sequence)
            .map_err(|e| ExtractError::FastqParse(format!("non-UTF8 sequence: {e}")))?;

        let spans = self
            .group_spans(seq_str)
            .ok_or(ExtractError::RegexNoMatch)?;

        // Collect named group spans into (name, start, end) sorted by name
        let mut umi_spans: Vec<(&str, usize, usize)> = Vec::new();
        let mut cell_spans: Vec<(&str, usize, usize)> = Vec::new();
        let mut discard_spans: Vec<(usize, usize)> = Vec::new();

        for (name, start, end) in spans {
            if name.starts_with("umi_") {
                umi_spans.push((name, start, end));
            } else if name.starts_with("cell_") {
                cell_spans.push((name, start, end));
            } else if name.starts_with("discard_") {
                discard_spans.push((start, end));
            }
        }

        // Sort by group name for deterministic concatenation
        umi_spans.sort_by_key(|&(name, _, _)| name);
        cell_spans.sort_by_key(|&(name, _, _)| name);

        // Build extracted-position bitmask (O(n) lookup instead of O(n*m) Vec::contains)
        let mut extracted = vec![false; sequence.len()];
        for &(_, start, end) in &umi_spans {
            extracted[start..end].fill(true);
        }
        for &(_, start, end) in &cell_spans {
            extracted[start..end].fill(true);
        }
        for &(start, end) in &discard_spans {
            extracted[start..end].fill(true);
        }

        // Build UMI and cell by concatenating group values in sorted name order
        let mut umi = Vec::new();
        let mut umi_quality = Vec::new();
        for &(_, start, end) in &umi_spans {
            umi.extend_from_slice(&sequence[start..end]);
            umi_quality.extend_from_slice(&quality[start..end]);
        }

        let mut cell_barcode = Vec::new();
        for &(_, start, end) in &cell_spans {
            cell_barcode.extend_from_slice(&sequence[start..end]);
        }

        // Build trimmed sequence/quality: keep positions not in any extraction set
        let mut trimmed_sequence = Vec::new();
        let mut trimmed_quality = Vec::new();

        for (i, &is_extracted) in extracted.iter().enumerate() {
            if !is_extracted {
                trimmed_sequence.push(sequence[i]);
                trimmed_quality.push(quality[i]);
            }
        }

        Ok(ExtractionResult {
            umi,
            umi_quality,
            cell_barcode,
            trimmed_sequence,
            trimmed_quality,
        })
    }
}

const MAX_FUZZY_ALTERNATIVES: usize = 10_000;
/// Room for a group at the alternative limit; the default 10 MiB holds about a third of that.
const EXACT_PATTERN_SIZE_LIMIT: usize = 64 << 20;

enum Atom {
    Single,
    Group,
}

/// Rewrite the `{s<=N}` substitution quantifiers of Python's `regex` module as plain regex.
///
/// A quantified single character, class or escape matches any character. A quantified
/// group must have a fixed length and becomes an alternation over the substituted positions.
fn preprocess_fuzzy(pattern_str: &str) -> Result<String, ExtractError> {
    rewrite_fuzzy(pattern_str, &mut |group, limit, _| {
        fuzzy_group(group, limit)
    })
    .map(|rewritten| rewritten.pattern)
}

struct Rewritten {
    pattern: String,
    alternates_at_top_level: bool,
}

/// `rewrite_group` receives each fuzzy group, its substitution limit, and whether the
/// group sits outside every other group with no quantifier of its own.
fn rewrite_fuzzy(
    pattern_str: &str,
    rewrite_group: &mut dyn FnMut(&str, usize, bool) -> Result<String, ExtractError>,
) -> Result<Rewritten, ExtractError> {
    let mut alternates_at_top_level = false;
    let chars: Vec<char> = pattern_str.chars().collect();
    let mut out = String::with_capacity(pattern_str.len());
    let mut open_groups: Vec<usize> = Vec::new();
    let mut last: Option<(usize, Atom)> = None;
    let mut i = 0;

    while i < chars.len() {
        let start = out.len();
        match chars[i] {
            '{' => {
                let Some(close) = chars[i..].iter().position(|&c| c == '}') else {
                    out.extend(&chars[i..]);
                    break;
                };
                let quantifier: String = chars[i..=i + close].iter().collect();
                i += close + 1;
                match substitution_limit(&quantifier)? {
                    None => out.push_str(&quantifier),
                    Some(0) => {}
                    Some(limit) => match last {
                        None => {
                            return Err(ExtractError::InvalidPattern(format!(
                                "fuzzy quantifier {quantifier} has nothing to apply to"
                            )));
                        }
                        Some((atom_start, Atom::Single)) => {
                            out.truncate(atom_start);
                            out.push('.');
                        }
                        Some((atom_start, Atom::Group)) => {
                            let required = open_groups.is_empty()
                                && !matches!(chars.get(i), Some('?' | '*' | '+' | '{'));
                            let group = rewrite_group(&out[atom_start..], limit, required)?;
                            out.truncate(atom_start);
                            out.push_str(&group);
                        }
                    },
                }
                last = None;
            }
            '\\' => {
                out.extend(&chars[i..chars.len().min(i + 2)]);
                i += 2;
                last = Some((start, Atom::Single));
            }
            '[' => {
                let end = class_end(&chars, i);
                out.extend(&chars[i..end]);
                i = end;
                last = Some((start, Atom::Single));
            }
            '(' => {
                open_groups.push(start);
                out.push('(');
                i += 1;
                last = None;
            }
            ')' => {
                out.push(')');
                i += 1;
                last = open_groups
                    .pop()
                    .map(|group_start| (group_start, Atom::Group));
            }
            c @ ('*' | '+' | '?' | '|' | '^' | '$') => {
                alternates_at_top_level |= c == '|' && open_groups.is_empty();
                out.push(c);
                i += 1;
                last = None;
            }
            c => {
                out.push(c);
                i += 1;
                last = Some((start, Atom::Single));
            }
        }
    }

    Ok(Rewritten {
        pattern: out,
        alternates_at_top_level,
    })
}

/// The N of `{s<=N}`, or `None` for a counted repetition such as `{3}` or `{8,12}`.
fn substitution_limit(quantifier: &str) -> Result<Option<usize>, ExtractError> {
    let body = &quantifier[1..quantifier.len() - 1];
    if body.bytes().all(|b| b.is_ascii_digit() || b == b',') {
        return Ok(None);
    }
    body.strip_prefix("s<=")
        .filter(|digits| !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()))
        .ok_or_else(|| {
            ExtractError::InvalidPattern(format!(
                "fuzzy quantifier {quantifier} is not supported: only substitutions, {{s<=N}}, are"
            ))
        })?
        .parse()
        .map(Some)
        .map_err(|_| {
            ExtractError::InvalidPattern(format!(
                "fuzzy quantifier {quantifier} exceeds the supported integer range"
            ))
        })
}

/// Index just past the `]` closing the character class that opens at `start`.
fn class_end(chars: &[char], start: usize) -> usize {
    let mut i = start + 1;
    if chars.get(i) == Some(&'^') {
        i += 1;
    }
    if chars.get(i) == Some(&']') {
        i += 1;
    }
    while i < chars.len() {
        match chars[i] {
            '\\' => i += 2,
            ']' => return i + 1,
            _ => i += 1,
        }
    }
    chars.len()
}

/// The header, such as `?P<name>`, and the single-character items of a parenthesised group.
fn fuzzy_group_items(group: &str) -> Result<(&str, Vec<String>), ExtractError> {
    let unsupported = || {
        ExtractError::InvalidPattern(format!(
            "fuzzy matching of {group} is not supported: the group must be a fixed-length run of \
             characters, classes and escapes"
        ))
    };
    let inner = &group[1..group.len() - 1];
    let header_len = if inner.starts_with("?:") {
        2
    } else if inner.starts_with("?P<") || inner.starts_with("?<") {
        inner.find('>').ok_or_else(unsupported)? + 1
    } else if inner.starts_with('?') {
        return Err(unsupported());
    } else {
        0
    };
    let (header, body) = inner.split_at(header_len);

    let chars: Vec<char> = body.chars().collect();
    let mut items: Vec<String> = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let end = match chars[i] {
            '\\' => chars.len().min(i + 2),
            '[' => class_end(&chars, i),
            '(' | ')' | '*' | '+' | '?' | '|' | '^' | '$' | '{' => return Err(unsupported()),
            _ => i + 1,
        };
        let item: String = chars[i..end].iter().collect();
        i = end;
        let mut count = 1;
        if chars.get(i) == Some(&'{') {
            let close = chars[i..]
                .iter()
                .position(|&c| c == '}')
                .ok_or_else(unsupported)?;
            let digits: String = chars[i + 1..i + close].iter().collect();
            count = digits.parse().map_err(|_| unsupported())?;
            i += close + 1;
        }
        if items.len() + count > MAX_FUZZY_ALTERNATIVES {
            return Err(unsupported());
        }
        items.extend(std::iter::repeat_n(item, count));
    }
    Ok((header, items))
}

/// Rewrite `group` as an alternation that lets up to `limit` of its items match any character.
fn fuzzy_group(group: &str, limit: usize) -> Result<String, ExtractError> {
    let (header, items) = fuzzy_group_items(group)?;
    let wildcards = limit.min(items.len());
    let mut alternatives: Vec<String> = Vec::new();
    let mut chosen: Vec<usize> = (0..wildcards).collect();
    loop {
        if alternatives.len() == MAX_FUZZY_ALTERNATIVES {
            return Err(ExtractError::InvalidPattern(format!(
                "fuzzy matching of {group} with {limit} substitutions needs more than \
                 {MAX_FUZZY_ALTERNATIVES} alternatives"
            )));
        }
        alternatives.push(
            items
                .iter()
                .enumerate()
                .map(|(position, item)| {
                    if chosen.contains(&position) {
                        "."
                    } else {
                        item
                    }
                })
                .collect(),
        );
        // Advance to the next combination of wildcard positions.
        let Some(movable) =
            (0..wildcards).rfind(|&slot| chosen[slot] < items.len() - wildcards + slot)
        else {
            break;
        };
        chosen[movable] += 1;
        for slot in movable + 1..wildcards {
            chosen[slot] = chosen[slot - 1] + 1;
        }
    }

    Ok(format!("({header}(?:{}))", alternatives.join("|")))
}

/// If `positions` is a contiguous ascending sequence [a, a+1, ..., b-1], return Some(a..b).
fn as_contiguous_range(positions: &[usize]) -> Option<Range<usize>> {
    let start = *positions.first()?;
    let is_contiguous = positions
        .iter()
        .enumerate()
        .skip(1)
        .all(|(i, &pos)| pos == start + i);
    is_contiguous.then(|| start..start + positions.len())
}

fn extract_slice<'a>(
    source: &'a [u8],
    range: Option<&Range<usize>>,
    positions: &[usize],
) -> Cow<'a, [u8]> {
    match range {
        Some(range) => Cow::Borrowed(&source[range.clone()]),
        None if positions.is_empty() => Cow::Borrowed(&[]),
        None => Cow::Owned(positions.iter().map(|&i| source[i]).collect()),
    }
}

fn join_slices(a: &[u8], b: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(a.len() + b.len());
    result.extend_from_slice(a);
    result.extend_from_slice(b);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- StringPattern tests ---

    #[test]
    fn contiguous_extraction_borrows_all_five_fields() {
        let sequence = b"ACGTACGT";
        let quality = b"12345678";
        for (end, cell, umi, remaining) in [
            (PrimeEnd::Five, 0..2, 2..4, 4..8),
            (PrimeEnd::Three, 4..6, 6..8, 0..4),
        ] {
            let pattern = StringPattern::parse("CCNN", end).unwrap();
            let view = pattern.extract_view(sequence, quality).unwrap();
            for (actual, expected) in [
                (&view.cell_barcode, &sequence[cell]),
                (&view.umi, &sequence[umi.clone()]),
                (&view.umi_quality, &quality[umi]),
                (&view.trimmed_sequence, &sequence[remaining.clone()]),
                (&view.trimmed_quality, &quality[remaining]),
            ] {
                assert_eq!(actual.as_ref(), expected);
                assert!(
                    matches!(actual, Cow::Borrowed(_)),
                    "copied an unchanged slice"
                );
                assert_eq!(actual.as_ptr(), expected.as_ptr());
            }
        }
    }

    #[test]
    fn extraction_matches_position_reference_for_all_short_patterns() {
        for length in 1..=4u32 {
            for mut code in 0..3usize.pow(length) {
                let pattern: String = (0..length)
                    .map(|_| {
                        let ch = char::from(b"NCX"[code % 3]);
                        code /= 3;
                        ch
                    })
                    .collect();
                for end in [PrimeEnd::Five, PrimeEnd::Three] {
                    for extra in [0, 3] {
                        let size = pattern.len() + extra;
                        let sequence = &b"ACGTNAC"[..size];
                        let quality = &b"1234567"[..size];
                        let start = if end == PrimeEnd::Five { 0 } else { extra };
                        let marker = |index: usize| {
                            index
                                .checked_sub(start)
                                .and_then(|offset| pattern.as_bytes().get(offset))
                                .copied()
                                .unwrap_or(b'X')
                        };
                        let select = |source: &[u8], kind| -> Vec<u8> {
                            source
                                .iter()
                                .enumerate()
                                .filter(|&(index, _)| marker(index) == kind)
                                .map(|(_, &base)| base)
                                .collect()
                        };
                        let result = StringPattern::parse(&pattern, end)
                            .unwrap()
                            .extract_view(sequence, quality)
                            .unwrap();
                        assert_eq!(result.umi.as_ref(), select(sequence, b'N'));
                        assert_eq!(result.umi_quality.as_ref(), select(quality, b'N'));
                        assert_eq!(result.cell_barcode.as_ref(), select(sequence, b'C'));
                        assert_eq!(result.trimmed_sequence.as_ref(), select(sequence, b'X'));
                        assert_eq!(result.trimmed_quality.as_ref(), select(quality, b'X'));
                    }
                }
            }
        }
    }

    #[test]
    fn parse_valid_pattern() {
        let pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        assert_eq!(pat.umi_positions, vec![0, 1, 2, 7, 8]);
        assert_eq!(pat.sample_positions, vec![3, 4, 5, 6]);
        assert_eq!(pat.cell_positions, []);
        assert_eq!(pat.pattern_length, 9);
    }

    #[test]
    fn parse_pattern_with_cell() {
        let pat = StringPattern::parse("CCCNNNNXXXX", PrimeEnd::Five).unwrap();
        assert_eq!(pat.cell_positions, vec![0, 1, 2]);
        assert_eq!(pat.umi_positions, vec![3, 4, 5, 6]);
        assert_eq!(pat.sample_positions, vec![7, 8, 9, 10]);
    }

    #[test]
    fn parse_invalid_pattern() {
        assert!(StringPattern::parse("NNNZXXNN", PrimeEnd::Five).is_err());
        assert!(StringPattern::parse("", PrimeEnd::Five).is_err());
    }

    #[test]
    fn extract_5prime_nnnxxxxnn() {
        let pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        let seq = b"CAGGTTCAATCTCGGTGGGACCTC";
        let qual = b"1=DFFFFHHHHHJJJFGIJIJJIJ";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.umi, b"CAGAA");
        assert_eq!(result.umi_quality, b"1=DHH");
        assert_eq!(result.cell_barcode, b"");
        assert_eq!(result.trimmed_sequence, b"GTTCTCTCGGTGGGACCTC");
        assert_eq!(result.trimmed_quality, b"FFFFHHHJJJFGIJIJJIJ");
    }

    #[test]
    fn extract_read_too_short() {
        let pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        assert!(pat.extract(b"ACGT", b"IIII").is_err());
    }

    #[test]
    fn extract_3prime() {
        let pat = StringPattern::parse("NNXX", PrimeEnd::Three).unwrap();
        let seq = b"ACGTAATTGG";
        let qual = b"IIIIIIIIII";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.umi, b"TT");
        assert_eq!(result.trimmed_sequence, b"ACGTAAGG");
    }

    // --- RegexPattern tests ---

    #[test]
    fn regex_matching_starts_at_the_beginning_of_the_read() {
        let pattern = RegexPattern::parse("(?P<umi_1>AA)").unwrap();
        assert!(matches!(
            pattern.extract(b"TTAA", b"IIII"),
            Err(ExtractError::RegexNoMatch)
        ));
        assert_eq!(pattern.extract(b"AATT", b"IIII").unwrap().umi, b"AA");

        let pattern = RegexPattern::parse(".*(?P<umi_1>AA)").unwrap();
        assert_eq!(pattern.extract(b"TTAA", b"IIII").unwrap().umi, b"AA");
    }

    #[test]
    fn regex_parse_valid() {
        let pat = RegexPattern::parse(r"^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})").unwrap();
        assert!(pat.pattern.is_match("CAGGTTCAATCTCGGTGGGACCTC"));
    }

    #[test]
    fn regex_parse_no_barcode_groups() {
        assert!(RegexPattern::parse(r"^(.{3}).{4}(.{2})").is_err());
    }

    #[test]
    fn regex_parse_invalid_regex() {
        assert!(RegexPattern::parse(r"^(?P<umi_1>.{3").is_err());
    }

    #[test]
    fn oversized_fuzzy_quantifier_returns_an_error() {
        assert!(matches!(
            RegexPattern::parse("(?P<umi_1>A{s<=9999999999999999999999999999999})"),
            Err(ExtractError::InvalidPattern(_))
        ));
    }

    fn fuzzy(pattern: &str) -> regex::Regex {
        regex::Regex::new(&format!("^(?:{})$", preprocess_fuzzy(pattern).unwrap())).unwrap()
    }

    #[test]
    fn fuzzy_group_allows_substitutions_anywhere_in_the_group() {
        let pattern = fuzzy("(?P<discard_1>GATTACA){s<=2}");
        for accepted in ["GATTACA", "CATTACA", "GATTACT", "TATTACT", "GAAAACA"] {
            assert!(pattern.is_match(accepted), "{accepted}");
        }
        for rejected in ["CTTTACT", "GATTAC", "GATTACAA", "ATTACAG"] {
            assert!(!pattern.is_match(rejected), "{rejected}");
        }
        let captures = pattern.captures("CATTACT").unwrap();
        assert_eq!(&captures["discard_1"], "CATTACT");
    }

    #[test]
    fn fuzzy_group_keeps_its_kind() {
        assert!(fuzzy("(?:ACG){s<=1}").is_match("ATG"));
        assert_eq!(&fuzzy("(ACG){s<=1}").captures("ATG").unwrap()[1], "ATG");
        assert_eq!(
            &fuzzy("(?<x>ACG){s<=1}").captures("ATG").unwrap()["x"],
            "ATG"
        );
    }

    #[test]
    fn fuzzy_group_handles_classes_escapes_and_counts() {
        let pattern = fuzzy(r"([AC]G{2}\dT){s<=1}");
        for accepted in ["AGG1T", "CGG9T", "TGG1T", "AGGxT", "AGT1T"] {
            assert!(pattern.is_match(accepted), "{accepted}");
        }
        for rejected in ["TGGxT", "ATT1T"] {
            assert!(!pattern.is_match(rejected), "{rejected}");
        }
    }

    #[test]
    fn fuzzy_group_with_enough_substitutions_matches_any_text_of_its_length() {
        assert!(fuzzy("(ACG){s<=3}").is_match("TTT"));
        assert!(fuzzy("(ACG){s<=9}").is_match("TTT"));
        assert!(!fuzzy("(ACG){s<=9}").is_match("TTTT"));
    }

    #[test]
    fn fuzzy_single_atom_matches_any_character() {
        for pattern in [
            "AC{s<=1}",
            r"A\d{s<=1}",
            "A[CG]{s<=1}",
            "A.{s<=1}",
            r"A[\]C]{s<=2}",
        ] {
            let pattern = fuzzy(pattern);
            assert!(pattern.is_match("AT"));
            assert!(!pattern.is_match("TT"));
        }
    }

    #[test]
    fn zero_substitutions_leave_the_pattern_exact() {
        assert!(!fuzzy("(ACG){s<=0}").is_match("ATG"));
        assert!(!fuzzy("AC{s<=0}").is_match("AT"));
        assert!(fuzzy(r"A\d{s<=0}").is_match("A1"));
    }

    #[test]
    fn widened_fast_path_agrees_with_the_exact_pattern() {
        for pattern in [
            "(?P<cell_1>.{2,5})(?P<discard_1>GA[TC]TA){s<=1}(?P<umi_1>.{3})(AC){s<=1}T*",
            "(?P<cell_1>.{2,5}?)(?P<discard_1>GATTACA){s<=2}(?P<umi_1>.{3})",
            "(?P<cell_1>[AC]{0,3})(GATT){s<=1}(?P<umi_1>.{2})(?:CAT){s<=1}",
        ] {
            assert_fast_path_agrees(pattern);
        }
    }

    fn assert_fast_path_agrees(pattern: &str) {
        let fast = RegexPattern::parse(pattern).unwrap();
        assert!(fast.relaxed.is_some());
        let exact = RegexPattern {
            relaxed: None,
            ..fast.clone()
        };
        let mut state = 7_u32;
        let mut matched = 0;
        for _ in 0..20_000 {
            let sequence: Vec<u8> = (0..18)
                .map(|_| {
                    state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    b"GATTACAC"[(state >> 29) as usize]
                })
                .collect();
            let quality = vec![b'I'; sequence.len()];
            let expected = exact.extract(&sequence, &quality).ok();
            matched += usize::from(expected.is_some());
            assert_eq!(
                fast.extract(&sequence, &quality).ok().map(|r| (
                    r.cell_barcode,
                    r.umi,
                    r.trimmed_sequence
                )),
                expected.map(|r| (r.cell_barcode, r.umi, r.trimmed_sequence)),
                "{pattern} on {}",
                String::from_utf8_lossy(&sequence)
            );
        }
        assert!(matched > 100, "{pattern}: {matched}");
    }

    #[test]
    fn fast_path_needs_every_fuzzy_group_in_every_match() {
        for pattern in [
            "(?P<umi_1>.{3})(ACG){s<=1}?",
            "(?P<umi_1>.{3})(ACG){s<=1}{2}",
            "(?P<umi_1>.{3})((ACG){s<=1})",
            "(?P<umi_1>.{3})T|(ACG){s<=1}",
            "(?P<umi_1>.{3})(ACG){s<=0}",
            "(?P<umi_1>.{3})(ACG){s<=3}",
            "(?P<umi_1>.{3})ACG{s<=1}",
        ] {
            assert!(
                RelaxedPattern::parse(pattern).unwrap().is_none(),
                "{pattern}"
            );
            assert!(RegexPattern::parse(pattern).is_ok(), "{pattern}");
        }
        assert!(
            RelaxedPattern::parse("(?P<umi_1>.{3})(ACG){s<=1}")
                .unwrap()
                .is_some()
        );
    }

    #[test]
    fn counted_repetition_is_not_a_fuzzy_quantifier() {
        assert_eq!(
            preprocess_fuzzy("(?P<umi_1>.{6})T{3,}.{2,4}").unwrap(),
            "(?P<umi_1>.{6})T{3,}.{2,4}"
        );
    }

    #[test]
    fn unsupported_fuzzy_forms_are_named_in_the_error() {
        for (pattern, named) in [
            ("(ACG){e<=1}", "{e<=1}"),
            ("(ACG){i<=1}", "{i<=1}"),
            ("(ACG){d<=1,s<=1}", "{d<=1,s<=1}"),
            ("(AC+G){s<=1}", "(AC+G)"),
            ("(A|C){s<=1}", "(A|C)"),
            ("(A(C)G){s<=1}", "(A(C)G)"),
            ("{s<=1}", "{s<=1}"),
            ("A*{s<=1}", "{s<=1}"),
        ] {
            let Err(ExtractError::InvalidPattern(message)) = preprocess_fuzzy(pattern) else {
                panic!("{pattern} was accepted");
            };
            assert!(message.contains(named), "{pattern}: {message}");
        }
    }

    #[test]
    fn regex_extract_equivalent_to_string() {
        // Regex ^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2}) should produce same result as NNNXXXXNN
        let string_pat = StringPattern::parse("NNNXXXXNN", PrimeEnd::Five).unwrap();
        let regex_pat = RegexPattern::parse(r"^(?P<umi_1>.{3}).{4}(?P<umi_2>.{2})").unwrap();

        let seq = b"CAGGTTCAATCTCGGTGGGACCTC";
        let qual = b"1=DFFFFHHHHHJJJFGIJIJJIJ";

        let string_result = string_pat.extract(seq, qual).unwrap();
        let regex_result = regex_pat.extract(seq, qual).unwrap();

        assert_eq!(string_result.umi, regex_result.umi);
        assert_eq!(string_result.cell_barcode, regex_result.cell_barcode);
        assert_eq!(
            string_result.trimmed_sequence,
            regex_result.trimmed_sequence
        );
        assert_eq!(string_result.trimmed_quality, regex_result.trimmed_quality);
    }

    #[test]
    fn regex_extract_with_cell() {
        let pat =
            RegexPattern::parse(r"^(?P<cell_1>.{3})(?P<umi_1>.{4})(?P<discard_1>.{2})").unwrap();

        let seq = b"ABCDEFGHIJKLM";
        let qual = b"1234567890ABC";

        let result = pat.extract(seq, qual).unwrap();

        assert_eq!(result.cell_barcode, b"ABC");
        assert_eq!(result.umi, b"DEFG");
        // Positions 0-8 extracted/discarded, remaining: JKLM (positions 9-12)
        assert_eq!(result.trimmed_sequence, b"JKLM");
        assert_eq!(result.trimmed_quality, b"0ABC");
    }

    #[test]
    fn regex_no_match() {
        let pat = RegexPattern::parse(r"^(?P<umi_1>ZZZZZ)").unwrap();
        let result = pat.extract(b"ACGTACGT", b"IIIIIIII");
        assert!(matches!(result, Err(ExtractError::RegexNoMatch)));
    }
}
