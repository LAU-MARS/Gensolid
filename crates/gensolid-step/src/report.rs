//! Validity reporting that distinguishes unsupported inputs from invalid bodies.

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};

/// What was established about a STEP input by this build.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepStatus {
    /// Every decoded body passed every enabled check.
    Valid,
    /// At least one decoded body failed a check.
    Invalid,
    /// No faults found, but a required geometric check was unavailable.
    Unchecked,
    /// The file needs an unsupported or unimplemented feature.
    Unsupported,
    /// The file is not a well-formed supported STEP document.
    Parse,
}

/// A report for any STEP text, including unreadable and unsupported documents.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StepValidityReport {
    /// Overall outcome; never treats an unsupported input as invalid geometry.
    pub status: StepStatus,
    /// Bodies decoded into the model.
    pub decoded_bodies: usize,
    /// Bodies on which the immutable checker ran.
    pub checked_bodies: usize,
    /// Number of explicit checker faults.
    pub faults: usize,
    /// Number of unavailable geometric checks.
    pub unchecked: usize,
    /// Diagnostics including source entity identifiers where available.
    pub messages: Vec<String>,
}

/// Produces a validity report without conflating unsupported input with invalid geometry.
pub fn validity_report(text: &str, precision: Precision) -> StepValidityReport {
    let mut report = StepValidityReport {
        status: StepStatus::Valid,
        decoded_bodies: 0,
        checked_bodies: 0,
        faults: 0,
        unchecked: 0,
        messages: Vec::new(),
    };
    let model = match crate::read_unchecked(text, precision) {
        Ok(model) => model,
        Err(error) => {
            report.status = match error.code() {
                ErrorCode::Unimplemented | ErrorCode::Unsupported => StepStatus::Unsupported,
                ErrorCode::InvalidArgument | ErrorCode::Degenerate | ErrorCode::OutOfRange => {
                    StepStatus::Invalid
                }
                _ => StepStatus::Parse,
            };
            report.messages.push(error.to_string());
            return report;
        }
    };
    let mut invalid_relationship = false;
    let mut unsupported_relationship = false;
    for body in model.bodies() {
        report.decoded_bodies += 1;
        match check_body(&model, body, &CheckOptions::default()) {
            Ok(checked) => {
                if checked.is_valid()
                    && let Err(error) = crate::voids::validate(&model, body)
                {
                    match error.code() {
                        ErrorCode::InvalidArgument
                        | ErrorCode::Degenerate
                        | ErrorCode::OutOfRange => invalid_relationship = true,
                        ErrorCode::Unimplemented | ErrorCode::Unsupported => {
                            unsupported_relationship = true
                        }
                        _ => report.unchecked += 1,
                    }
                    report.messages.push(error.to_string());
                }
                report.checked_bodies += 1;
                report.faults += checked.faults.len();
                report.unchecked += checked.unchecked;
                report
                    .messages
                    .extend(checked.faults.into_iter().map(|fault| format!("{fault:?}")));
            }
            Err(error) => {
                report.unchecked += 1;
                report.messages.push(error.to_string());
            }
        }
    }
    if report.faults != 0 || invalid_relationship {
        report.status = StepStatus::Invalid;
    } else if unsupported_relationship {
        report.status = StepStatus::Unsupported;
    } else if report.unchecked != 0 {
        report.status = StepStatus::Unchecked;
    }
    report
}
