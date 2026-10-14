//! Minimal Euler-style edge subdivision with all shared co-edge uses preserved.

use std::collections::HashSet;

use gensolid_core::{Error, Result, effective_tolerance};
use gensolid_geom::{Curve, Line, ParamRange, Surface};
use gensolid_math::{Interval, Sign};

use crate::entity::*;
use crate::model::{Capacity, Model};
use crate::polyhedron::{dot, vector};

/// Handles produced by subdividing a bounded straight edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EdgeSplit {
    /// Original edge, retained as the first carrier segment.
    pub original_edge: EdgeId,
    /// Second segment, sharing the original geometric carrier and tolerance.
    pub new_edge: EdgeId,
    /// New vertex common to both segments.
    pub vertex: VertexId,
}

fn carrier(model: &Model, edge: EdgeId) -> Result<(Line, VertexId, VertexId)> {
    let edge = model.edge(edge)?;
    let (Some(start), Some(end)) = (edge.start, edge.end) else {
        return Err(Error::unimplemented(
            "edge subdivision requires two endpoint vertices",
        ));
    };
    let curve = edge
        .curve
        .ok_or_else(|| Error::invalid_argument("edge has no geometric carrier"))?;
    let Curve::Line(line) = model.geom().curve(curve)? else {
        return Err(Error::unimplemented(
            "edge subdivision with curved or periodic parameter bounds is not implemented",
        ));
    };
    Ok((*line, start, end))
}

/// Arc-length parameter bounds of a straight edge on its stored line carrier.
/// Curved edges require explicit trimming/periodic semantics and return Unimplemented.
pub fn edge_parameter_bounds(model: &Model, edge: EdgeId) -> Result<ParamRange> {
    let (line, start, end) = carrier(model, edge)?;
    let start = model.vertex(start)?.point;
    let end = model.vertex(end)?.point;
    let direction = line.direction();
    let delta = dot(
        vector(end, start),
        [
            Interval::point(direction.x),
            Interval::point(direction.y),
            Interval::point(direction.z),
        ],
    );
    if delta.sign() != Some(Sign::Positive) {
        return Err(Error::degenerate(
            "edge endpoints do not have certified increasing carrier parameters",
        ));
    }
    let lo = line.closest_param(start);
    let hi = line.closest_param(end);
    if !(lo.is_finite() && hi.is_finite())
        || (Interval::point(hi) - Interval::point(lo)).sign() != Some(Sign::Positive)
    {
        return Err(Error::unimplemented(
            "edge parameter span is not representable on its line carrier",
        ));
    }
    Ok(ParamRange::bounded(lo, hi))
}

/// Splits a bounded line edge at a strict interior fraction of its carrier range.
///
/// The original edge retains its ID and first segment; the second edge shares its
/// line and tolerance. Every positive fin is followed by its new fin; every negative
/// fin is preceded by the corresponding reversed new fin. Wire-shell edge lists are
/// updated too. Incident faces must be planar and contain the bounded carrier at the
/// edge tolerance. Ring and curved edges are refused until explicit trim semantics
/// exist. All preconditions are checked before mutation; a higher layer must run
/// `check_body` on affected bodies before returning this topology-layer edit.
pub fn split_line_edge(model: &mut Model, edge: EdgeId, fraction: f64) -> Result<EdgeSplit> {
    if !(fraction.is_finite() && fraction > 0.0 && fraction < 1.0) {
        return Err(Error::invalid_argument(
            "split fraction must lie strictly between zero and one",
        ));
    }
    let original = model.edge(edge)?.clone();
    for &fin in &original.fins {
        if model.fin_geometry(fin)?.is_some() {
            return Err(Error::unimplemented(
                "line subdivision with explicit pcurves requires the trimmed-face editor",
            ));
        }
    }
    let (line, start, end) = carrier(model, edge)?;
    let range = edge_parameter_bounds(model, edge)?;
    let precision = *model.precision();
    let parameter = range.lo + (range.hi - range.lo) * fraction;
    let point = line.eval(parameter).point;
    if !point.is_finite() {
        return Err(Error::degenerate("split point is not finite"));
    }
    for vertex in [start, end] {
        let delta = vector(point, model.vertex(vertex)?.point);
        if (dot(delta, delta) - Interval::point(precision.linear()).sqr()).sign()
            != Some(Sign::Positive)
        {
            return Err(Error::degenerate(
                "split point is not separated from both endpoints at session precision",
            ));
        }
    }
    let edge_tolerance = effective_tolerance(original.tolerance, &precision);
    let mut fin_uses = Vec::new();
    let mut seen = HashSet::new();
    for fin in &original.fins {
        if !seen.insert(*fin) {
            return Err(Error::invalid_argument("edge repeats a fin reference"));
        }
        let fin_value = model.fin(*fin)?;
        if fin_value.edge != edge {
            return Err(Error::invalid_argument(
                "edge fin does not point back to the edge",
            ));
        }
        let lp = model.loop_(fin_value.owner)?;
        if lp.fins.iter().filter(|id| **id == *fin).count() != 1 {
            return Err(Error::invalid_argument(
                "fin must occur exactly once in its owner loop",
            ));
        }
        let face = model.face(lp.face)?;
        let surface = face
            .surface
            .ok_or_else(|| Error::invalid_argument("incident face has no surface"))?;
        let Surface::Plane(plane) = model.geom().surface(surface)? else {
            return Err(Error::unimplemented(
                "line subdivision on nonplanar incident faces is not implemented",
            ));
        };
        let normal = plane.frame().z();
        let normal = [
            Interval::point(normal.x),
            Interval::point(normal.y),
            Interval::point(normal.z),
        ];
        for parameter in [range.lo, range.hi, parameter] {
            let distance = dot(
                normal,
                vector(line.eval(parameter).point, plane.frame().origin()),
            );
            match (distance.sqr() - Interval::point(edge_tolerance).sqr()).sign() {
                Some(Sign::Negative | Sign::Zero) => (),
                Some(Sign::Positive) => {
                    return Err(Error::invalid_argument(
                        "bounded edge carrier is outside an incident plane's tolerance",
                    ));
                }
                None => {
                    return Err(Error::unimplemented(
                        "edge carrier membership is unresolved at its tolerance",
                    ));
                }
            }
        }
        fin_uses.push((*fin, fin_value.clone()));
    }
    let wire_shells = model.wire_shells(edge);
    for shell in &wire_shells {
        if model
            .shell(*shell)?
            .wire_edges
            .iter()
            .filter(|id| **id == edge)
            .count()
            != 1
        {
            return Err(Error::invalid_argument("wire shell repeats the edge"));
        }
    }
    model.ensure_capacity(&Capacity {
        fins: fin_uses.len(),
        edges: 1,
        vertices: 1,
        ..Capacity::default()
    })?;
    let mut raw = model.raw();
    let vertex = raw.insert_vertex(Vertex {
        point,
        tolerance: original.tolerance,
    })?;
    let new_edge = raw.insert_edge(Edge {
        fins: Vec::new(),
        start: Some(vertex),
        end: Some(end),
        curve: original.curve,
        tolerance: original.tolerance,
    })?;
    raw.edge_mut(edge)?.end = Some(vertex);
    for (fin_id, fin) in fin_uses {
        let new_fin = raw.insert_fin(Fin {
            owner: fin.owner,
            edge: new_edge,
            sense: fin.sense,
        })?;
        let lp = raw.loop_mut(fin.owner)?;
        let position = lp
            .fins
            .iter()
            .position(|id| *id == fin_id)
            .ok_or_else(|| Error::internal("validated fin disappeared from its loop"))?;
        lp.fins.insert(
            position + usize::from(fin.sense == Sense::Positive),
            new_fin,
        );
        raw.edge_mut(new_edge)?.fins.push(new_fin);
    }
    for shell in wire_shells {
        let shell = raw.shell_mut(shell)?;
        let position = shell
            .wire_edges
            .iter()
            .position(|id| *id == edge)
            .ok_or_else(|| Error::internal("validated edge disappeared from its wire shell"))?;
        shell.wire_edges.insert(position + 1, new_edge);
    }
    Ok(EdgeSplit {
        original_edge: edge,
        new_edge,
        vertex,
    })
}
