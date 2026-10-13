//! Invoke OCCT in a separate process so invalid supported imports can be compared
//! without a checked-import shortcut discarding them. The helper reuses the
//! immutable verifier's binding selector and BRepCheck_Analyzer.

use std::path::Path;
use std::process::Command;

use gensolid_core::{Error, Result};
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
pub(crate) struct Observation {
    pub(crate) read: bool,
    pub(crate) valid: bool,
    expected_solids: usize,
    transferred_solids: usize,
    partial: bool,
}

pub(crate) fn step(path: &Path) -> Result<Observation> {
    let python = ["python", "python3"]
        .into_iter()
        .find(|name| {
            Command::new(name)
                .arg("--version")
                .output()
                .is_ok_and(|output| output.status.success())
        })
        .ok_or_else(|| Error::unsupported("Python executable not found"))?;
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occt-diagnostic.py");
    let output = Command::new(python)
        .arg(script)
        .arg("--step")
        .arg(path)
        .output()
        .map_err(|error| Error::internal(format!("starting OCCT helper: {error}")))?;
    if !output.status.success() {
        return Err(Error::internal(format!(
            "OCCT helper failed ({}): {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    let report: Observation = serde_json::from_slice(&output.stdout)
        .map_err(|error| Error::parse(format!("OCCT helper returned malformed JSON: {error}")))?;
    if report.read
        && (report.partial
            || report.expected_solids == 0
            || report.expected_solids != report.transferred_solids)
    {
        return Err(Error::parse(
            "OCCT reported a complete read without transferring every source solid",
        ));
    }
    Ok(report)
}
