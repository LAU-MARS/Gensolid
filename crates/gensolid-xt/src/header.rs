//! Text XT header.
//!
//! A text transmit file starts with two character-set check lines (so transcoding between
//! code pages is detected), then three header parts of `KEY=VALUE;` items, then the
//! end-of-header marker, then the transmit data.

use gensolid_core::{Error, Result};

const CHARSET_LINE_1: &str = "**ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
const CHARSET_LINE_2: &str = "**PARASOLID !\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~0123456789";
const PART_MARKERS: [&str; 3] = ["**PART1;", "**PART2;", "**PART3;"];
const END_OF_HEADER: &str = "**END_OF_HEADER";
const LINE_WIDTH: usize = 80;

/// Header items, in file order, of the three header parts.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TextHeader {
    /// Part 1 items (machine and application information).
    pub part1: Vec<(String, String)>,
    /// Part 2 items (schema and user field information).
    pub part2: Vec<(String, String)>,
    /// Part 3 items (application-defined).
    pub part3: Vec<(String, String)>,
}

impl TextHeader {
    /// Value of the first item named `key` in any part.
    pub fn get(&self, key: &str) -> Option<&str> {
        [&self.part1, &self.part2, &self.part3]
            .into_iter()
            .flatten()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

/// A parsed text transmit file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TextFile<'a> {
    /// Header.
    pub header: TextHeader,
    /// Modeller version from the transmit preamble, if present.
    pub modeller_version: Option<String>,
    /// Schema key (`SCH_...`) from the transmit preamble, if present.
    pub schema: Option<String>,
    /// Everything after the header, starting with the transmit preamble.
    pub data: &'a str,
}

fn parse_items(text: &str) -> Result<Vec<(String, String)>> {
    text.split(';')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(|item| {
            let (key, value) = item
                .split_once('=')
                .ok_or_else(|| Error::parse(format!("header item without '=': {item:?}")))?;
            Ok((key.trim().to_owned(), value.trim().to_owned()))
        })
        .collect()
}

/// Parses a text transmit file's header and preamble.
pub fn read_text(input: &str) -> Result<TextFile<'_>> {
    let mut parts: [String; 3] = Default::default();
    let mut current: Option<usize> = None;
    let mut data_start = None;
    let mut offset = 0;
    for (index, line) in input.split_inclusive('\n').enumerate() {
        offset += line.len();
        let trimmed = line.trim_end();
        let expected_charset = match index {
            0 => Some(CHARSET_LINE_1),
            1 => Some(CHARSET_LINE_2),
            _ => None,
        };
        if let Some(expected) = expected_charset {
            if !trimmed
                .strip_prefix(expected)
                .is_some_and(|padding| padding.bytes().all(|byte| byte == b'*'))
            {
                return Err(Error::parse(
                    "character-set check lines missing or altered; not text XT, or transcoded",
                ));
            }
            continue;
        }
        if let Some(padding) = trimmed.strip_prefix(END_OF_HEADER) {
            if !padding.bytes().all(|byte| byte == b'*') {
                return Err(Error::parse("altered end-of-header marker"));
            }
            data_start = Some(offset);
            break;
        }
        if let Some(index) = PART_MARKERS.iter().position(|m| trimmed.starts_with(m)) {
            if current.is_some_and(|c| index <= c) {
                return Err(Error::parse("header parts out of order"));
            }
            current = Some(index);
            parts[index].push_str(&trimmed[PART_MARKERS[index].len()..]);
            continue;
        }
        let Some(index) = current else {
            return Err(Error::parse(format!("unexpected header line: {trimmed:?}")));
        };
        parts[index].push_str(trimmed);
    }
    let start = data_start.ok_or_else(|| Error::parse("missing end-of-header marker"))?;
    let data = &input[start..];

    let [p1, p2, p3] = &parts;
    let header = TextHeader {
        part1: parse_items(p1)?,
        part2: parse_items(p2)?,
        part3: parse_items(p3)?,
    };

    let mut modeller_version = None;
    let mut schema = None;
    let mut previous: [Option<&str>; 2] = [None, None];
    for word in data.split_whitespace().take(64) {
        if previous == [Some("modeller"), Some("version")] && modeller_version.is_none() {
            modeller_version = Some(word.to_owned());
        }
        if word.starts_with("SCH_") {
            schema = Some(word.to_owned());
            break;
        }
        previous = [previous[1], Some(word)];
    }
    // Real text XT packs length-prefixed strings together and wraps physical
    // lines inside those strings. Ordinary whitespace words are only a legacy
    // header diagnostic; the strict stream reader supplies authoritative values.
    if let Ok(internal) = parasolid_core::inspect_xt(
        input.as_bytes(),
        parasolid_core::InspectionLimits::default(),
    ) {
        modeller_version = internal
            .modeller_version
            .split_whitespace()
            .last()
            .map(str::to_owned);
        schema = Some(internal.schema_key);
    }
    Ok(TextFile {
        header,
        modeller_version,
        schema,
        data,
    })
}

fn write_part(out: &mut String, marker: &str, items: &[(String, String)]) {
    let mut line = String::from(marker);
    for (key, value) in items {
        let item = format!("{key}={value};");
        if line.len() + item.len() > LINE_WIDTH {
            out.push_str(&line);
            out.push('\n');
            line.clear();
        }
        line.push_str(&item);
    }
    out.push_str(&line);
    out.push('\n');
}

fn padded(prefix: &str) -> String {
    let mut line = prefix.to_owned();
    while line.len() < LINE_WIDTH {
        line.push('*');
    }
    line
}

/// Writes the character-set check lines, the three header parts and the end-of-header
/// marker, wrapping items at 80 columns. Items must not contain `;` or newlines.
pub fn write_text_header(header: &TextHeader) -> Result<String> {
    for (key, value) in [&header.part1, &header.part2, &header.part3]
        .into_iter()
        .flatten()
    {
        let key_ok = !key.is_empty() && !key.contains(['=', ';', '\n', '\r']);
        if !key_ok || value.contains([';', '\n', '\r']) {
            return Err(Error::invalid_argument(format!(
                "header item {key:?} cannot be encoded"
            )));
        }
    }
    let mut out = String::new();
    out.push_str(&padded(CHARSET_LINE_1));
    out.push('\n');
    out.push_str(&padded(CHARSET_LINE_2));
    out.push('\n');
    write_part(&mut out, PART_MARKERS[0], &header.part1);
    write_part(&mut out, PART_MARKERS[1], &header.part2);
    write_part(&mut out, PART_MARKERS[2], &header.part3);
    out.push_str(&padded(END_OF_HEADER));
    out.push('\n');
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_header() -> TextHeader {
        let item = |k: &str, v: &str| (k.to_owned(), v.to_owned());
        TextHeader {
            part1: vec![
                item("MC", "x86_64"),
                item("MC_MODEL", "generic"),
                item("OS", "windows"),
                item("APPL", "gensolid"),
                item("SITE", "test"),
                item("USER", "verifier"),
                item("FORMAT", "text"),
                item("GUISE", "transmit"),
                item("KEY", "fixture"),
                item("FILE", "block.x_t"),
            ],
            part2: vec![item("SCH", "SCH_0000000_00000"), item("USFLD_SIZE", "0")],
            part3: Vec::new(),
        }
    }

    #[test]
    fn header_round_trips() {
        let header = sample_header();
        let mut text = write_text_header(&header).unwrap();
        text.push_str(
            "T51 : TRANSMIT FILE created by modeller version 0000000\nSCH_0000000_00000\n1 2 3\n",
        );
        let file = read_text(&text).unwrap();
        assert_eq!(file.header, header);
        assert_eq!(file.modeller_version.as_deref(), Some("0000000"));
        assert_eq!(file.schema.as_deref(), Some("SCH_0000000_00000"));
        assert!(file.data.starts_with("T51 : TRANSMIT FILE"));
        assert_eq!(file.header.get("USFLD_SIZE"), Some("0"));
        assert!(text.lines().all(|l| l.len() <= LINE_WIDTH));
    }

    #[test]
    fn rejects_transcoded_files() {
        let text = write_text_header(&sample_header())
            .unwrap()
            .replacen("[\\]", "[/]", 1);
        assert!(read_text(&text).is_err());
    }

    #[test]
    fn rejects_missing_end_marker() {
        let text = write_text_header(&sample_header()).unwrap();
        let truncated = &text[..text.find(END_OF_HEADER).unwrap()];
        assert!(read_text(truncated).is_err());
    }
}
