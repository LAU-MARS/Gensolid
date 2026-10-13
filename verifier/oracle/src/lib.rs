//! Differential tests against Open CASCADE, run as a separate process.
//!
//! The public cases live in `verifier/corpus`. Each case is built with Gensolid and, when
//! pythonocc is installed, with OCCT. OCCT is never linked: the helper
//! `occt_block.py` is a subprocess, so the LGPL kernel stays outside this crate.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use gensolid_api::{Axis2, Session, SessionOptions};
use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::analytic_intersection_score;
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_step::{read_planar_solids, write_block_step};
use gensolid_xt::{XtValidityReport, validity_report};
use serde::Deserialize;

/// One solid block in the public corpus.
#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct BlockCase {
    /// Case name.
    pub name: String,
    /// Size along the frame X axis.
    pub x: f64,
    /// Size along the frame Y axis.
    pub y: f64,
    /// Size along the frame Z axis.
    pub z: f64,
    /// Frame origin. The block is centred here.
    pub location: [f64; 3],
    /// Frame Z axis.
    pub axis: [f64; 3],
    /// Frame X reference direction.
    pub ref_direction: [f64; 3],
    /// Expected face count.
    pub faces: usize,
    /// Expected edge count.
    pub edges: usize,
    /// Expected vertex count.
    pub vertices: usize,
}

#[derive(Debug, Deserialize)]
struct CorpusFile {
    cases: Vec<BlockCase>,
}

/// Counts, volume and axis-aligned box of one block.
#[derive(Clone, Debug, PartialEq)]
pub struct BlockReport {
    /// Face count.
    pub faces: usize,
    /// Edge count.
    pub edges: usize,
    /// Vertex count.
    pub vertices: usize,
    /// Solid volume.
    pub volume: f64,
    /// Minimum corner of the axis-aligned box.
    pub box_min: [f64; 3],
    /// Maximum corner of the axis-aligned box.
    pub box_max: [f64; 3],
}

#[derive(Debug, Deserialize)]
struct OcctReport {
    faces: usize,
    edges: usize,
    vertices: usize,
    volume: f64,
    box_min: [f64; 3],
    box_max: [f64; 3],
}

/// Public block corpus.
pub fn block_cases() -> Result<Vec<BlockCase>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/blocks.json");
    let text = std::fs::read_to_string(&path).map_err(|err| {
        Error::new(
            gensolid_core::ErrorCode::Parse,
            format!("reading {}: {err}", path.display()),
        )
    })?;
    let file: CorpusFile = serde_json::from_str(&text)
        .map_err(|err| Error::parse(format!("parsing {}: {err}", path.display())))?;
    if file.cases.is_empty() {
        return Err(Error::parse("block corpus is empty"));
    }
    Ok(file.cases)
}

/// Build `case` in a fresh Gensolid session and measure it.
pub fn measure_gensolid(case: &BlockCase) -> Result<BlockReport> {
    let mut session = Session::start(SessionOptions::default());
    let basis = Axis2 {
        location: case.location,
        axis: case.axis,
        ref_direction: case.ref_direction,
    };
    let body = session.body_create_solid_block(case.x, case.y, case.z, Some(&basis))?;
    if !session.body_check(body)?.is_empty() {
        return Err(Error::new(
            gensolid_core::ErrorCode::Internal,
            format!("{} failed body_check", case.name),
        ));
    }
    let (box_min, box_max) = session.body_ask_bounding_box(body)?;
    Ok(BlockReport {
        faces: session.body_ask_faces(body)?.len(),
        edges: session.body_ask_edges(body)?.len(),
        vertices: session.body_ask_vertices(body)?.len(),
        volume: session.body_ask_volume(body)?,
        box_min,
        box_max,
    })
}

fn python() -> Option<&'static str> {
    for name in ["python", "python3"] {
        let found = Command::new(name)
            .arg("--version")
            .output()
            .is_ok_and(|output| output.status.success());
        if found {
            return Some(name);
        }
    }
    None
}

fn script() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("occt_block.py")
}

/// Whether the OCCT helper can import pythonocc.
pub fn occt_available() -> bool {
    let Some(python) = python() else {
        return false;
    };
    Command::new(python)
        .arg(script())
        .arg("--probe")
        .output()
        .is_ok_and(|output| output.status.success())
}

/// Measure the same block with OCCT.
pub fn measure_occt(case: &BlockCase) -> Result<BlockReport> {
    let python = python().ok_or_else(|| Error::unsupported("python was not found"))?;
    let output = Command::new(python)
        .arg(script())
        .arg("--x")
        .arg(case.x.to_string())
        .arg("--y")
        .arg(case.y.to_string())
        .arg("--z")
        .arg(case.z.to_string())
        .arg("--location")
        .arg(format_xyz(case.location))
        .arg("--axis")
        .arg(format_xyz(case.axis))
        .arg("--ref")
        .arg(format_xyz(case.ref_direction))
        .output()
        .map_err(|err| {
            Error::new(
                gensolid_core::ErrorCode::Internal,
                format!("starting OCCT helper: {err}"),
            )
        })?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(Error::new(
            gensolid_core::ErrorCode::Internal,
            format!("OCCT helper failed for {}: {stderr}", case.name),
        ));
    }
    let parsed: OcctReport = serde_json::from_slice(&output.stdout).map_err(|err| {
        Error::parse(format!(
            "OCCT helper returned {}: {err}",
            String::from_utf8_lossy(&output.stdout)
        ))
    })?;
    Ok(BlockReport {
        faces: parsed.faces,
        edges: parsed.edges,
        vertices: parsed.vertices,
        volume: parsed.volume,
        box_min: parsed.box_min,
        box_max: parsed.box_max,
    })
}

fn format_xyz(value: [f64; 3]) -> String {
    format!("{},{},{}", value[0], value[1], value[2])
}

/// How often Gensolid and OCCT agree that a STEP solid is valid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValidityAgreement {
    /// Files both kernels read.
    pub compared: usize,
    /// Files whose valid/invalid verdict matched.
    pub agreed: usize,
    /// Files at least one kernel could not read.
    pub unreadable: usize,
}

impl ValidityAgreement {
    /// `agreed / compared`, or `None` when nothing was compared.
    pub fn rate(self) -> Option<f64> {
        if self.compared == 0 {
            None
        } else {
            Some(self.agreed as f64 / self.compared as f64)
        }
    }
}

/// STEP files in `dir` (`.step` and `.stp`).
pub fn step_files(dir: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    let Ok(entries) = fs::read_dir(dir) else {
        return files;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let step = path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("step") || ext.eq_ignore_ascii_case("stp"));
        if step {
            files.push(path);
        }
    }
    files.sort();
    files
}

/// Directory for an ABC STEP subset. Empty until files are added.
pub fn abc_directory() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../corpus/abc")
}

/// Validity report for one XT document. Every input produces a report.
pub fn report_xt(input: &str) -> XtValidityReport {
    validity_report(input)
}

/// Pass rate of `verifier/private/blocks.json`, if that file exists.
///
/// The cases themselves are not returned.
pub fn private_block_score() -> Result<Option<f64>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../private/blocks.json");
    if !path.exists() {
        return Ok(None);
    }
    let text = fs::read_to_string(&path)
        .map_err(|err| Error::parse(format!("reading private corpus: {err}")))?;
    let file: CorpusFile = serde_json::from_str(&text)
        .map_err(|err| Error::parse(format!("parsing private corpus: {err}")))?;
    if file.cases.is_empty() {
        return Ok(Some(1.0));
    }
    let precision = Precision::DEFAULT;
    let mut passed = 0usize;
    for case in &file.cases {
        let report = measure_gensolid(case)?;
        let expected = case.x * case.y * case.z;
        let counts = report.faces == case.faces
            && report.edges == case.edges
            && report.vertices == case.vertices;
        if counts && (report.volume - expected).abs() <= volume_tolerance(case, &precision) {
            passed += 1;
        }
    }
    Ok(Some(passed as f64 / file.cases.len() as f64))
}

fn check_script() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("occt_check.py")
}

/// Whether OCCT's STEP reader can be imported.
pub fn occt_step_available() -> bool {
    let Some(python) = python() else {
        return false;
    };
    Command::new(python)
        .arg(check_script())
        .arg("--probe")
        .output()
        .is_ok_and(|output| output.status.success())
}

#[derive(Deserialize)]
struct OcctValidity {
    read: bool,
    valid: bool,
}

fn occt_step_valid(path: &Path) -> Result<Option<bool>> {
    let python = python().ok_or_else(|| Error::unsupported("python was not found"))?;
    let output = Command::new(python)
        .arg(check_script())
        .arg("--step")
        .arg(path)
        .output()
        .map_err(|err| {
            Error::new(
                gensolid_core::ErrorCode::Internal,
                format!("starting OCCT STEP check: {err}"),
            )
        })?;
    if !output.status.success() {
        return Err(Error::new(
            gensolid_core::ErrorCode::Internal,
            format!(
                "OCCT STEP check failed: {}",
                String::from_utf8_lossy(&output.stderr)
            ),
        ));
    }
    let parsed: OcctValidity = serde_json::from_slice(&output.stdout).map_err(|err| {
        Error::parse(format!(
            "OCCT STEP check returned {}: {err}",
            String::from_utf8_lossy(&output.stdout)
        ))
    })?;
    Ok(parsed.read.then_some(parsed.valid))
}

fn gensolid_step_valid(text: &str) -> Result<Option<bool>> {
    let model = match read_planar_solids(text, Precision::DEFAULT) {
        Ok(model) => model,
        Err(_) => return Ok(None),
    };
    let mut valid = true;
    let mut bodies = 0usize;
    for body in model.bodies() {
        bodies += 1;
        if !check_body(&model, body, &CheckOptions::default())?.is_valid() {
            valid = false;
        }
    }
    Ok((bodies > 0).then_some(valid))
}

/// Compare Gensolid and OCCT on STEP files already on disk.
pub fn agree_step_files(paths: &[PathBuf]) -> Result<ValidityAgreement> {
    let mut agreement = ValidityAgreement {
        compared: 0,
        agreed: 0,
        unreadable: 0,
    };
    for path in paths {
        let text = fs::read_to_string(path)
            .map_err(|err| Error::parse(format!("reading {}: {err}", path.display())))?;
        let ours = gensolid_step_valid(&text)?;
        let theirs = occt_step_valid(path)?;
        match (ours, theirs) {
            (Some(ours), Some(theirs)) => {
                agreement.compared += 1;
                if ours == theirs {
                    agreement.agreed += 1;
                }
            }
            _ => agreement.unreadable += 1,
        }
    }
    Ok(agreement)
}

/// Write each block to a temporary STEP file and compare validity verdicts.
pub fn agree_blocks(cases: &[BlockCase]) -> Result<ValidityAgreement> {
    let precision = Precision::DEFAULT;
    let mut paths = Vec::new();
    for case in cases {
        let frame = Frame::new(
            Point3::new(case.location[0], case.location[1], case.location[2]),
            Vec3::new(case.axis[0], case.axis[1], case.axis[2]),
            Vec3::new(
                case.ref_direction[0],
                case.ref_direction[1],
                case.ref_direction[2],
            ),
            &precision,
        )?;
        let step = write_block_step(case.x, case.y, case.z, &frame, precision)?;
        let path = std::env::temp_dir().join(format!("gensolid-{}.step", case.name));
        fs::write(&path, step).map_err(|err| {
            Error::new(
                gensolid_core::ErrorCode::Internal,
                format!("writing {}: {err}", path.display()),
            )
        })?;
        paths.push(path);
    }
    let agreement = agree_step_files(&paths);
    for path in &paths {
        let _ = fs::remove_file(path);
    }
    agreement
}

/// M1 score for analytic intersections in two frames. Exposed so the verifier records the
/// same completeness number as the geometry tests.
pub fn m1_intersection_score(precision: &Precision) -> Result<gensolid_geom::IntersectionScore> {
    let frames = [
        Frame::from_axis(Point3::ORIGIN, Vec3::Z, precision)?,
        Frame::from_axis(
            Point3::new(1.0, 2.0, 3.0),
            Vec3::new(0.2, -0.4, 0.9),
            precision,
        )?,
    ];
    analytic_intersection_score(&frames, precision)
}

/// Absolute volume tolerance: one session-linear displacement of the whole surface.
pub fn volume_tolerance(case: &BlockCase, precision: &Precision) -> f64 {
    let area = 2.0 * (case.x * case.y + case.y * case.z + case.z * case.x);
    precision.linear() * area
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn corpus_blocks_match_analytic_volume_and_counts() {
        let precision = Precision::DEFAULT;
        let cases = block_cases().unwrap();
        for case in cases {
            let report = measure_gensolid(&case).unwrap();
            assert_eq!(report.faces, case.faces, "{}", case.name);
            assert_eq!(report.edges, case.edges, "{}", case.name);
            assert_eq!(report.vertices, case.vertices, "{}", case.name);
            let expected = case.x * case.y * case.z;
            assert!(
                (report.volume - expected).abs() <= volume_tolerance(&case, &precision),
                "{} volume {} expected {expected}",
                case.name,
                report.volume
            );
        }
    }

    #[test]
    fn corpus_blocks_match_occt_when_installed() {
        if !occt_available() {
            eprintln!(
                "occt oracle skipped: pythonocc is not importable (pip install cadquery, or conda install pythonocc-core)"
            );
            return;
        }
        let precision = Precision::DEFAULT;
        for case in block_cases().unwrap() {
            let ours = measure_gensolid(&case).unwrap();
            let theirs = measure_occt(&case).unwrap();
            assert_eq!(ours.faces, theirs.faces, "{}", case.name);
            assert_eq!(ours.edges, theirs.edges, "{}", case.name);
            assert_eq!(ours.vertices, theirs.vertices, "{}", case.name);
            let volume_tol = volume_tolerance(&case, &precision);
            assert!(
                (ours.volume - theirs.volume).abs() <= volume_tol,
                "{} volume {} vs OCCT {}",
                case.name,
                ours.volume,
                theirs.volume
            );
            for (label, a, b) in [
                ("box min", ours.box_min, theirs.box_min),
                ("box max", ours.box_max, theirs.box_max),
            ] {
                for axis in 0..3 {
                    assert!(
                        (a[axis] - b[axis]).abs() <= precision.linear(),
                        "{} {label}[{axis}] {} vs OCCT {}",
                        case.name,
                        a[axis],
                        b[axis]
                    );
                }
            }
        }
    }

    #[test]
    fn xt_report_covers_unreadable_input() {
        let report = report_xt("not a transmit file");
        assert_eq!(report.decoded_bodies, 0);
        assert!(!report.message.is_empty());
    }

    #[test]
    fn abc_subset_is_scanned_even_when_empty() {
        let files = step_files(&abc_directory());
        if files.is_empty() {
            assert_eq!(agree_step_files(&files).unwrap().rate(), None);
        }
    }

    #[test]
    fn private_corpus_is_absent_or_scored_without_echoing_cases() {
        let score = private_block_score().unwrap();
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../private/blocks.json");
        if path.exists() {
            let rate = score.unwrap();
            assert!((0.0..=1.0).contains(&rate));
        } else {
            assert!(score.is_none());
        }
    }

    #[test]
    #[allow(clippy::float_cmp)] // Completeness is a ratio of equal integer counts.
    fn analytic_intersection_completeness_is_one() {
        let precision = Precision::DEFAULT;
        let score = m1_intersection_score(&precision).unwrap();
        assert_eq!(score.completeness(), 1.0);
        assert!(score.max_point_error <= precision.linear(), "{score:?}");
    }

    #[test]
    #[allow(clippy::float_cmp)] // The rate is agreed/compared for a full match.
    fn step_validity_matches_occt_when_installed() {
        if !occt_step_available() {
            eprintln!("occt step check skipped: pythonocc STEP reader is not importable");
            return;
        }
        let agreement = agree_blocks(&block_cases().unwrap()).unwrap();
        assert!(agreement.compared >= 4, "{agreement:?}");
        assert_eq!(agreement.unreadable, 0, "{agreement:?}");
        assert_eq!(agreement.rate(), Some(1.0), "{agreement:?}");
    }
}
