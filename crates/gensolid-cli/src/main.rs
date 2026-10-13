//! File diagnostics and development scores. The human-owned L0 verifier remains the
//! authority for milestone acceptance; this command never declares L0 acceptance.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use gensolid_core::{Error, Precision, Result};
use gensolid_step::StepStatus;
use gensolid_xt::XtStatus;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

mod catalog;
mod m2;
mod manifest;
mod oracle;
mod paired;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
enum Verdict {
    Valid,
    Invalid,
    Unsupported,
    Unchecked,
    Parse,
    Io,
    Integrity,
}

#[derive(Debug, Serialize)]
struct FileReport {
    path: String,
    format: String,
    status: Verdict,
    decoded_bodies: usize,
    detail: Value,
}

#[derive(Debug, Default, Serialize)]
struct Coverage {
    files: usize,
    valid: usize,
    invalid: usize,
    unsupported: usize,
    unchecked: usize,
    parse: usize,
    io: usize,
    integrity: usize,
    definitive_verdicts: usize,
    verdict_coverage: Option<f64>,
}

/// Aggregate-only receipt from an independently hosted private scoring service.
#[derive(Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct PrivateScore {
    suite_id: String,
    total: usize,
    passed: usize,
}

impl PrivateScore {
    fn read(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .map_err(|error| Error::parse(format!("reading private score: {error}")))?;
        let score: Self = serde_json::from_str(&text)
            .map_err(|error| Error::parse(format!("parsing private score: {error}")))?;
        score.validate()?;
        Ok(score)
    }

    fn validate(&self) -> Result<()> {
        if self.suite_id.trim().is_empty() || self.total == 0 || self.passed > self.total {
            return Err(Error::invalid_argument(
                "private score requires a nonempty suite_id, total > 0 and passed <= total",
            ));
        }
        Ok(())
    }
}

impl Coverage {
    fn from_reports(reports: &[FileReport]) -> Self {
        let mut result = Self::default();
        for report in reports {
            result.files += 1;
            match report.status {
                Verdict::Valid => result.valid += 1,
                Verdict::Invalid => result.invalid += 1,
                Verdict::Unsupported => result.unsupported += 1,
                Verdict::Unchecked => result.unchecked += 1,
                Verdict::Parse => result.parse += 1,
                Verdict::Io => result.io += 1,
                Verdict::Integrity => result.integrity += 1,
            }
        }
        result.definitive_verdicts = result.valid + result.invalid;
        result.verdict_coverage = ratio(result.definitive_verdicts, result.files);
        result
    }
}

fn ratio(numerator: usize, denominator: usize) -> Option<f64> {
    (denominator != 0).then(|| numerator as f64 / denominator as f64)
}

fn format(path: &Path) -> Option<&'static str> {
    let extension = path.extension()?.to_str()?.to_ascii_lowercase();
    match extension.as_str() {
        "step" | "stp" => Some("step"),
        "x_t" | "x_b" | "xmt_txt" | "xmt_bin" => Some("xt"),
        _ => None,
    }
}

fn inspect_file(path: &Path) -> FileReport {
    let mut report = FileReport {
        path: path.display().to_string(),
        format: format(path).unwrap_or("unknown").to_owned(),
        status: Verdict::Io,
        decoded_bodies: 0,
        detail: Value::Null,
    };
    let bytes = match fs::read(path) {
        Ok(bytes) => bytes,
        Err(error) => {
            report.detail = json!({"message": error.to_string()});
            return report;
        }
    };
    match report.format.as_str() {
        "xt" => {
            let checked = gensolid_xt::validity_report_bytes(&bytes);
            report.status = match checked.status {
                XtStatus::Checked => Verdict::Valid,
                XtStatus::Invalid => Verdict::Invalid,
                XtStatus::UnsupportedSchema | XtStatus::UnsupportedGeometry => Verdict::Unsupported,
                XtStatus::Parse => Verdict::Parse,
            };
            report.decoded_bodies = checked.decoded_bodies;
            report.detail = json!({
                "xt_status": format!("{:?}", checked.status),
                "schema": checked.schema,
                "modeller_version": checked.modeller_version,
                "data_tokens": checked.data_tokens,
                "message": checked.message,
            });
        }
        "step" => {
            let text = match std::str::from_utf8(&bytes) {
                Ok(text) => text,
                Err(error) => {
                    report.status = Verdict::Parse;
                    report.detail = json!({"message": error.to_string()});
                    return report;
                }
            };
            let checked = gensolid_step::validity_report(text, Precision::DEFAULT);
            report.status = match checked.status {
                StepStatus::Valid => Verdict::Valid,
                StepStatus::Invalid => Verdict::Invalid,
                StepStatus::Unchecked => Verdict::Unchecked,
                StepStatus::Unsupported => Verdict::Unsupported,
                StepStatus::Parse => Verdict::Parse,
            };
            report.decoded_bodies = checked.decoded_bodies;
            report.detail = json!({
                "checked_bodies": checked.checked_bodies,
                "faults": checked.faults,
                "unchecked": checked.unchecked,
                "messages": checked.messages,
            });
        }
        _ => {
            report.status = Verdict::Unsupported;
            report.detail =
                json!({"message": "expected .step, .stp, .x_t, .x_b, .xmt_txt or .xmt_bin"});
        }
    }
    report
}

fn file_paths(path: &Path, paths: &mut Vec<PathBuf>) -> Result<()> {
    let metadata = fs::symlink_metadata(path)
        .map_err(|error| Error::parse(format!("reading {}: {error}", path.display())))?;
    // Directory links can introduce cycles. Explicit file links are still readable.
    if metadata.file_type().is_symlink() {
        if path.is_file() && format(path).is_some() {
            paths.push(path.to_path_buf());
        }
        return Ok(());
    }
    if metadata.is_file() {
        if format(path).is_some() {
            paths.push(path.to_path_buf());
        }
        return Ok(());
    }
    for entry in fs::read_dir(path)
        .map_err(|error| Error::parse(format!("listing {}: {error}", path.display())))?
    {
        let entry = entry.map_err(|error| Error::parse(error.to_string()))?;
        file_paths(&entry.path(), paths)?;
    }
    Ok(())
}

fn output(value: &Value, destination: Option<&Path>) -> Result<()> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|error| Error::internal(format!("serializing report: {error}")))?;
    if let Some(path) = destination {
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)
                .map_err(|error| Error::parse(format!("creating {}: {error}", parent.display())))?;
        }
        fs::write(path, format!("{text}\n"))
            .map_err(|error| Error::parse(format!("writing {}: {error}", path.display())))?;
    } else {
        println!("{text}");
    }
    Ok(())
}

fn check(paths: &[String]) -> Result<u8> {
    if paths.is_empty() {
        return Err(Error::invalid_argument(
            "check requires at least one file or directory",
        ));
    }
    let mut files = Vec::new();
    for path in paths {
        let path = Path::new(path);
        if path.is_dir() {
            file_paths(path, &mut files)?;
        } else {
            files.push(path.to_path_buf());
        }
    }
    files.sort();
    files.dedup();
    let reports: Vec<_> = files.iter().map(|path| inspect_file(path)).collect();
    let coverage = Coverage::from_reports(&reports);
    let status = if coverage.invalid != 0 {
        1
    } else if coverage.files == 0 || coverage.valid != coverage.files {
        2
    } else {
        0
    };
    output(&json!({"coverage": coverage, "reports": reports}), None)?;
    Ok(status)
}

fn benchmark(arguments: &[String]) -> Result<u8> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut corpus = root.join("corpus");
    let mut destination = None;
    let mut private_receipt = None;
    let mut iter = arguments.iter();
    while let Some(flag) = iter.next() {
        let value = iter
            .next()
            .ok_or_else(|| Error::invalid_argument(format!("missing value for {flag}")))?;
        match flag.as_str() {
            "--corpus" => corpus = PathBuf::from(value),
            "--output" => destination = Some(PathBuf::from(value)),
            "--private-score" => private_receipt = Some(PrivateScore::read(Path::new(value))?),
            _ => {
                return Err(Error::invalid_argument(format!(
                    "unknown benchmark option {flag}"
                )));
            }
        }
    }
    let mut files = Vec::new();
    file_paths(&corpus, &mut files)?;
    let frozen = manifest::FrozenInputs::read(&corpus, &mut files)?;
    let violations = frozen.violations();
    files.sort();
    files.dedup();
    let reports: Vec<_> = files.iter().map(|path| {
        let mut report=inspect_file(path);
        if let Some(message)=violations.get(path)
            && report.status != Verdict::Io
        {
            let import_status=report.status;
            report.status=Verdict::Integrity;
            report.detail=json!({"integrity_error":message,"import_status":import_status,"import_report":report.detail});
        }
        report
    }).collect();
    let coverage = Coverage::from_reports(&reports);
    let occt_available = gensolid_oracle::occt_step_available();
    let mut step_inputs = 0;
    let mut compared = 0;
    let mut agreed = 0;
    let mut oracle_errors = Vec::new();
    let mut oracle_rows = Vec::new();
    let mut oracle_readable = 0;
    let mut oracle_unreadable = 0;
    let mut valid_valid = 0;
    let mut invalid_invalid = 0;
    let mut valid_invalid = 0;
    let mut invalid_valid = 0;
    for (path, report) in files.iter().zip(&reports) {
        if report.format != "step" {
            continue;
        }
        step_inputs += 1;
        if !occt_available {
            continue;
        }
        match oracle::step(path) {
            Ok(observation) => {
                let theirs = observation.read.then_some(observation.valid);
                if theirs.is_some() {
                    oracle_readable += 1;
                } else {
                    oracle_unreadable += 1;
                }
                let ours = match report.status {
                    Verdict::Valid => Some(true),
                    Verdict::Invalid => Some(false),
                    _ => None,
                };
                let verdict_match = match (ours, theirs) {
                    (Some(ours), Some(theirs)) => {
                        compared += 1;
                        if ours == theirs {
                            agreed += 1;
                        }
                        match (ours, theirs) {
                            (true, true) => valid_valid += 1,
                            (false, false) => invalid_invalid += 1,
                            (true, false) => valid_invalid += 1,
                            (false, true) => invalid_valid += 1,
                        }
                        Some(ours == theirs)
                    }
                    _ => None,
                };
                oracle_rows.push(json!({"path":path,"gensolid_status":report.status,"occt_read":theirs.is_some(),"occt_valid":theirs,"occt_transfer":observation,"verdict_match":verdict_match}));
            }
            Err(error) => {
                oracle_rows.push(
                    json!({"path":path,"gensolid_status":report.status,"error":error.to_string()}),
                );
                oracle_errors.push(format!("{}: {error}", path.display()));
            }
        }
    }
    let blocks = gensolid_oracle::block_cases()?;
    let block_agreement = if occt_available {
        match gensolid_oracle::agree_blocks(&blocks) {
            Ok(score) => {
                json!({"compared":score.compared,"agreed":score.agreed,"unreadable":score.unreadable,"agreement":score.rate()})
            }
            Err(error) => json!({"error":error.to_string()}),
        }
    } else {
        json!({"status":"oracle_unavailable"})
    };
    let m1 = gensolid_oracle::m1_intersection_score(&Precision::DEFAULT)?;
    let seed = 20261013;
    let randomized = catalog::score(seed, 64, &Precision::DEFAULT)?;
    let extended = catalog::extended_score(seed, 64, &Precision::DEFAULT)?;
    let m2 = m2::run(&Precision::DEFAULT)?;
    let paired_directory = if corpus.join("abc-smoke").is_dir() {
        corpus.join("abc-smoke")
    } else {
        corpus.clone()
    };
    let paired = match paired::compare(&paired_directory, &Precision::DEFAULT) {
        Ok(score) => score,
        Err(error) => json!({"status":"unavailable","error":error.to_string()}),
    };
    output(
        &json!({
            "report_kind":"development_milestone_report",
            "acceptance_scope":"public_development_benchmark",
            "l0_private_acceptance":"not_established",
            "precision":{"linear":Precision::DEFAULT.linear(),"angular":Precision::DEFAULT.angular(),"size":Precision::DEFAULT.size()},
            "m2":m2,
            "m0":{
                "corpus":corpus,
                "coverage":coverage,
                "frozen_inputs":{"status":if !violations.is_empty() || !frozen.errors.is_empty() {"failed"} else if frozen.expected.is_empty() {"unfrozen"} else {"verified"},
                    "expected_files":frozen.expected.len(),"expected_cad_files":frozen.expected.keys().filter(|path|format(path).is_some()).count(),"integrity_failures":violations.len(),"errors":frozen.errors,
                    "violations":violations.iter().map(|(path,message)|json!({"path":path,"message":message})).collect::<Vec<_>>()},
                "files":reports,
                "occt_step":{"available":occt_available,"inputs":step_inputs,"compared":compared,"agreed":agreed,
                    "oracle_readable":oracle_readable,"oracle_unreadable":oracle_unreadable,"oracle_error_inputs":oracle_errors.len(),
                    "confusion_matrix":{"valid_valid":valid_valid,"invalid_invalid":invalid_invalid,"gensolid_valid_occt_invalid":valid_invalid,"gensolid_invalid_occt_valid":invalid_valid},
                "comparison_coverage":ratio(compared,step_inputs),"agreement":ratio(agreed,compared),"errors":oracle_errors,"files":oracle_rows},
                "occt_xt":{"status":"not_supported_by_occt","message":"OCCT is a STEP oracle; native XT needs independent producer/paired-file evidence."},
                "paired_xt_step":paired,
                "blocks":{"cases":blocks.len(),"occt":block_agreement},
                "private_blocks":{"status":"reserved","available":false,"score":null},
                "private_score_interface":{
                    "status":if private_receipt.is_some() { "aggregate_receipt_supplied" } else { "reserved" },
                    "receipt":private_receipt,
                    "message":"Receipts carry only suite_id/total/passed; private cases stay with the independent scoring service.",
                },
            },
            "m1":{
                "scope":"fixed analytic known-answer catalog in two rigid frames",
                "pairs":m1.pairs,"matched":m1.matched,"completeness":m1.completeness(),"max_point_error":m1.max_point_error,
                "general_surface_pairs":"unimplemented_pairs_remain",
                "randomized_known_answers":{
                    "scope":"fourteen known relative configurations at 64 random scales and common rigid frames",
                    "seed":seed,"frames":64,"pairs":randomized.pairs,"matched":randomized.matched,
                    "unimplemented":randomized.unimplemented,"failed":randomized.failed,
                    "completeness":randomized.completeness(),"max_point_error":randomized.max_point_error,
                },
                "extended_known_answers":{
                    "scope":"cubic curve-plane, bilinear surface-plane and oblique cone-plane; 64 random scales and common rigid frames",
                    "seed":seed,"pairs":extended.pairs,"matched":extended.matched,"unimplemented":extended.unimplemented,"failed":extended.failed,
                    "completeness":extended.completeness(),"max_point_error":(extended.pairs>extended.unimplemented+extended.failed).then_some(extended.max_point_error),
                },
                "combined_randomized_coverage":{
                    "pairs":randomized.pairs+extended.pairs,"matched":randomized.matched+extended.matched,
                    "unimplemented":randomized.unimplemented+extended.unimplemented,"failed":randomized.failed+extended.failed,
                    "completeness":ratio(randomized.matched+extended.matched,randomized.pairs+extended.pairs),
                },
            },
        }),
        destination.as_deref(),
    )?;
    // A successful report generation is separate from the recorded coverage/verdicts.
    Ok(if violations.is_empty() && frozen.errors.is_empty() {
        0
    } else {
        2
    })
}

fn run() -> Result<u8> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let Some((command, rest)) = arguments.split_first() else {
        return Err(Error::invalid_argument(
            "usage: gensolid-cli check <files/directories...> | benchmark [--corpus directory] [--output report.json] [--private-score receipt.json]",
        ));
    };
    match command.as_str() {
        "check" => check(rest),
        "benchmark" => benchmark(rest),
        _ => Err(Error::invalid_argument(format!(
            "unknown command {command}"
        ))),
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(code) => ExitCode::from(code),
        Err(error) => {
            eprintln!("{error}");
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_corpus_and_unsupported_files_have_no_success_score() {
        assert_eq!(Coverage::from_reports(&[]).verdict_coverage, None);
        let reports = [FileReport {
            path: "future.x_t".into(),
            format: "xt".into(),
            status: Verdict::Unsupported,
            decoded_bodies: 0,
            detail: Value::Null,
        }];
        let coverage = Coverage::from_reports(&reports);
        assert_eq!(coverage.definitive_verdicts, 0);
        assert_eq!(coverage.unsupported, 1);
        assert_eq!(coverage.verdict_coverage, Some(0.0));
    }

    #[test]
    fn valid_files_and_missing_files_remain_distinct() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
        let missing = inspect_file(&root.join("target/nonexistent-fixture.x_b"));
        assert_eq!(missing.status, Verdict::Io);
        let fixture = root.join("crates/gensolid-step/tests/fixtures/block.step");
        let valid = inspect_file(&fixture);
        assert_eq!(valid.status, Verdict::Valid);
        assert_eq!(valid.decoded_bodies, 1);
    }

    #[test]
    fn private_receipts_reject_empty_scores_and_embedded_cases() {
        let empty = PrivateScore {
            suite_id: "held-out-v1".into(),
            total: 0,
            passed: 0,
        };
        assert!(empty.validate().is_err());
        let impossible = PrivateScore {
            suite_id: "held-out-v1".into(),
            total: 2,
            passed: 3,
        };
        assert!(impossible.validate().is_err());
        assert!(
            serde_json::from_str::<PrivateScore>(
                r#"{"suite_id":"held-out-v1","total":2,"passed":1,"cases":[]}"#
            )
            .is_err()
        );
        assert!(
            serde_json::from_str::<PrivateScore>(
                r#"{"suite_id":"held-out-v1","total":2,"passed":1}"#
            )
            .is_ok_and(|receipt| receipt.validate().is_ok())
        );
    }
}
