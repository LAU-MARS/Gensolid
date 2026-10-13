//! A validity report for any XT file, including files this build cannot decode yet.
//!
//! A parsed header never implies a valid body. Only complete mapping followed by
//! `check_body` without faults or skipped checks produces `Checked`.

use crate::header::read_text;
use crate::import::{inspect, map_model, source_builtin};
use crate::tokens::tokenize;
use gensolid_core::{ErrorCode, Precision};

/// Why a transmit file did or did not yield bodies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XtStatus {
    /// The header parsed, but its exact schema or an encountered node type is unsupported.
    UnsupportedSchema,
    /// A known schema decoded, but the kernel cannot represent or check its geometry.
    UnsupportedGeometry,
    /// Every decoded body passed the kernel checker without skipped checks.
    Checked,
    /// Decoded topology or geometry failed validation.
    Invalid,
    /// The text is not a transmit file.
    Parse,
}

/// Validity report for one transmit file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XtValidityReport {
    /// Outcome.
    pub status: XtStatus,
    /// Schema key from the preamble, when the header parsed.
    pub schema: Option<String>,
    /// Modeller version from the preamble, when present.
    pub modeller_version: Option<String>,
    /// Tokens in the node stream after the header.
    pub data_tokens: usize,
    /// Bodies that were built and checked.
    pub decoded_bodies: usize,
    /// What was checked, or why it was not.
    pub message: String,
}

/// Reports on `input` without failing to produce a report.
pub fn validity_report(input: &str) -> XtValidityReport {
    let strict = validity_report_bytes(input.as_bytes());
    if strict.schema.is_some() || strict.status != XtStatus::Parse {
        return strict;
    }
    // Preserve the historical header-only diagnostic for an unknown schema, even for
    // legacy synthetic preambles that do not contain XT length prefixes. Known schemas
    // always go through the strict reader below. The byte API is strict for all inputs.
    if let Ok(file) = read_text(input)
        && let Some(schema) = &file.schema
        && schema != crate::schema::ABC_KEY
        && let Ok(key) = parasolid_core::SchemaKey::parse(schema)
        && let Ok(registry) = parasolid_core::BuiltinProfileRegistry::compiled()
        && registry.provider_for_key(&key).is_none()
    {
        return unsupported_legacy(input);
    }
    strict
}

/// Reports on arbitrary bytes, including neutral-binary XT and malformed files.
pub fn validity_report_bytes(input: &[u8]) -> XtValidityReport {
    validity_report_bytes_with_precision(input, Precision::DEFAULT)
}

/// Reports with a caller-selected session precision in source-file units.
pub fn validity_report_bytes_with_precision(
    input: &[u8],
    precision: Precision,
) -> XtValidityReport {
    let header = match inspect(input) {
        Ok(header) => header,
        Err(error) => {
            return XtValidityReport {
                status: XtStatus::Parse,
                schema: None,
                modeller_version: None,
                data_tokens: 0,
                decoded_bodies: 0,
                message: error.to_string(),
            };
        }
    };
    let mut report = XtValidityReport {
        status: XtStatus::Parse,
        schema: Some(header.schema.clone()),
        modeller_version: header.modeller.split_whitespace().last().map(str::to_owned),
        data_tokens: if header.binary {
            0
        } else {
            std::str::from_utf8(&input[header.data_start..])
                .map_or(0, |text| tokenize(text).count())
        },
        decoded_bodies: 0,
        message: String::new(),
    };
    let source = match source_builtin(input, &header) {
        Ok(source) => source,
        Err(error) => {
            report.status = match error.code() {
                ErrorCode::Unimplemented | ErrorCode::Unsupported => XtStatus::UnsupportedSchema,
                ErrorCode::InvalidArgument | ErrorCode::Degenerate => XtStatus::Invalid,
                _ => XtStatus::Parse,
            };
            report.message = error.to_string();
            return report;
        }
    };
    let decoded = match map_model(&source, precision) {
        Ok(decoded) => decoded,
        Err(error) => {
            report.status = if matches!(
                error.code(),
                ErrorCode::Unimplemented | ErrorCode::Unsupported
            ) {
                XtStatus::UnsupportedGeometry
            } else {
                XtStatus::Invalid
            };
            report.message = error.to_string();
            return report;
        }
    };
    report.decoded_bodies = decoded.bodies.len();
    let faults: usize = decoded.reports.iter().map(|r| r.faults.len()).sum();
    let unchecked: usize = decoded.reports.iter().map(|r| r.unchecked).sum();
    report.status = if faults != 0 {
        XtStatus::Invalid
    } else if unchecked != 0 {
        XtStatus::UnsupportedGeometry
    } else {
        XtStatus::Checked
    };
    report.message = format!(
        "checked {} bodies: {faults} faults, {unchecked} unchecked geometric checks",
        report.decoded_bodies
    );
    for (index, body) in decoded.reports.iter().enumerate() {
        for fault in &body.faults {
            report.message.push_str(&format!(
                "; body {index}: {:?}: {}",
                fault.kind, fault.message
            ));
        }
    }
    report
}

fn unsupported_legacy(input: &str) -> XtValidityReport {
    let file = match read_text(input) {
        Ok(file) => file,
        Err(error) => {
            return XtValidityReport {
                status: XtStatus::Parse,
                schema: None,
                modeller_version: None,
                data_tokens: 0,
                decoded_bodies: 0,
                message: error.to_string(),
            };
        }
    };
    let data_tokens = tokenize(file.data).count();
    let schema = file.schema.clone();
    let message = match &schema {
        Some(schema) => {
            format!("schema {schema} has no registered decoder; no body was built or checked")
        }
        None => "transmit file has no schema key; no body was built or checked".to_owned(),
    };
    XtValidityReport {
        status: XtStatus::UnsupportedSchema,
        schema,
        modeller_version: file.modeller_version,
        data_tokens,
        decoded_bodies: 0,
        message,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::{TextHeader, write_text_header};

    #[test]
    fn any_text_file_gets_a_report() {
        let garbage = validity_report("this is not an XT file");
        assert_eq!(garbage.status, XtStatus::Parse);
        assert_eq!(garbage.decoded_bodies, 0);

        let mut header = TextHeader::default();
        header.part1.push(("FORMAT".to_owned(), "text".to_owned()));
        let mut text = write_text_header(&header).unwrap();
        text.push_str(
            "T51 : TRANSMIT FILE created by modeller version 3500000\nSCH_3500000_15000\n1 2 3\n",
        );
        let report = validity_report(&text);
        assert_eq!(report.status, XtStatus::UnsupportedSchema);
        assert_eq!(report.schema.as_deref(), Some("SCH_3500000_15000"));
        assert_eq!(report.modeller_version.as_deref(), Some("3500000"));
        assert!(report.data_tokens >= 3);
        assert_eq!(report.decoded_bodies, 0);
        assert!(report.message.contains("no body was built"));
    }
}
