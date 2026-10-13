//! Pair the official ABC smoke model's native XT parts with its STEP solids.
//! Counts and vertex residuals are evidence about import fidelity, not a proof that
//! every trimmed face or curve is equivalent.

use std::fs;
use std::path::Path;

use gensolid_core::{Error, Precision, Result};
use gensolid_math::Point3;
use gensolid_topo::{BodyId, Model};
use serde_json::{Value, json};

struct Signature {
    faces: usize,
    edges: usize,
    vertices: Vec<Point3>,
}

fn signature(model: &Model, body: BodyId) -> Result<Signature> {
    Ok(Signature {
        faces: model.body_faces(body)?.len(),
        edges: model.body_edges(body)?.len(),
        vertices: model
            .body_vertices(body)?
            .into_iter()
            .map(|vertex| model.vertex(vertex).map(|vertex| vertex.point))
            .collect::<Result<Vec<_>>>()?,
    })
}

fn vertex_residual(a: &[Point3], b: &[Point3]) -> Option<f64> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let directed = |from: &[Point3], to: &[Point3]| {
        from.iter()
            .map(|point| {
                to.iter()
                    .map(|other| (*point - *other).length())
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max)
    };
    // STEP adds vertices to closed edges which are vertexless in native XT. Source
    // XT vertices must occur in STEP, while additional STEP ring vertices are retained.
    Some(directed(a, b))
}

pub(crate) fn compare(corpus: &Path, precision: &Precision) -> Result<Value> {
    let step = corpus.join("00000050_80d90bfdd2e74e709956122a_step_000.step");
    let parts: Vec<_> = (1..=3)
        .map(|part| corpus.join(format!("Part Studio 1 - Part {part}.x_t")))
        .collect();
    if !step.is_file() || parts.iter().any(|path| !path.is_file()) {
        return Ok(json!({"status":"paired_example_not_present"}));
    }
    let text = fs::read_to_string(&step).map_err(|error| Error::parse(error.to_string()))?;
    let model = gensolid_step::read_solids(&text, *precision)?;
    let step_signatures = model
        .bodies()
        .map(|body| signature(&model, body))
        .collect::<Result<Vec<_>>>()?;
    let mut used = vec![false; step_signatures.len()];
    let mut matched = 0;
    let mut vertex_comparisons = 0;
    let mut rows = Vec::new();
    for path in parts {
        let bytes = fs::read(&path).map_err(|error| Error::parse(error.to_string()))?;
        let imported = gensolid_xt::read_bytes(&bytes, *precision)?;
        for body in imported.bodies {
            let xt = signature(&imported.model, body)?;
            let candidate = step_signatures
                .iter()
                .enumerate()
                .filter(|(index, other)| {
                    !used[*index] && xt.faces == other.faces && xt.edges == other.edges
                })
                .map(|(index, other)| (index, vertex_residual(&xt.vertices, &other.vertices)))
                .min_by(|a, b| {
                    a.1.unwrap_or(f64::INFINITY)
                        .total_cmp(&b.1.unwrap_or(f64::INFINITY))
                });
            let pair = candidate
                .filter(|(_, error)| error.is_none_or(|error| error <= precision.linear()));
            if let Some((index, _)) = pair {
                used[index] = true;
                matched += 1;
            }
            if pair.is_some_and(|(_, error)| error.is_some()) {
                vertex_comparisons += 1;
            }
            rows.push(json!({"xt":path,"faces":xt.faces,"edges":xt.edges,"vertices":xt.vertices.len(),
                "matched_step_body":pair.map(|(index,_)|index),"max_vertex_residual":candidate.and_then(|(_,error)|error),
                "vertex_comparison":if xt.vertices.is_empty() { "no_vertices" } else { "measured" }}));
        }
    }
    Ok(
        json!({"scope":"one-to-one face/edge counts and XT-to-STEP vertex residual; STEP ring-edge vertices are retained and curved interiors are not compared",
        "step_bodies":step_signatures.len(),"xt_bodies":rows.len(),"face_edge_count_matches":matched,
        "step_counts":step_signatures.iter().map(|signature|json!({"faces":signature.faces,"edges":signature.edges,"vertices":signature.vertices.len()})).collect::<Vec<_>>(),
        "vertex_measured_matches":vertex_comparisons,"counts_only_matches":matched-vertex_comparisons,
        "all_face_edge_counts_matched":matched == step_signatures.len() && matched == rows.len(),"parts":rows}),
    )
}
