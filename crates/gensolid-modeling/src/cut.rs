//! Bounded plane-sheet cuts of ruled NURBS faces, with independent algebraic validation.
mod extruded;
mod rational;
use crate::{body_check_complete, body_split_face_isoparametric, boundary};
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{BSurface, Plane, Surface, backend::GeometryBackend};
use gensolid_math::{Frame, Interval, Point2, Point3, Sign, Vec3, predicates::orient2d};
use gensolid_topo::{BodyId, EdgeId, FaceId, Model};

/// Topology created by a bounded-sheet cut. Source bodies and tags are retained.
#[derive(Clone, Debug)]
pub struct FaceSplit {
    /// New body containing the split faces and propagated adjacent edge subdivisions.
    pub body: BodyId,
    /// Two replacement faces, sharing the new intersection edge.
    pub faces: [FaceId; 2],
    /// Shared three-dimensional intersection edge.
    pub intersection: EdgeId,
}

/// Explicit capability selection for bounded planar-sheet cuts.
#[derive(Clone, Copy, Debug, Default)]
pub struct SheetCutOptions {
    /// Allow certified oblique cuts of positive-weight, single-span Bezier
    /// translation extrusions against an affine planar tool patch. The complete
    /// transverse branch must stay inside the tool; trimmed, partial, tangent
    /// and multi-span branches are refused.
    /// Defaults to false, preserving the legacy cut capability contract.
    pub allow_rational_oblique: bool,
}
fn dot_delta(p: Point3, o: Point3, n: Vec3) -> Interval {
    (Interval::point(p.x) - Interval::point(o.x)) * Interval::point(n.x)
        + (Interval::point(p.y) - Interval::point(o.y)) * Interval::point(n.y)
        + (Interval::point(p.z) - Interval::point(o.z)) * Interval::point(n.z)
}
fn planar_tool(model: &Model, face: FaceId) -> Result<(Plane, Vec<Point2>, BSurface)> {
    let data = boundary::extract(model, face)?;
    if data.edges.len() != 4 {
        return Err(Error::unimplemented(
            "cutting sheet must be a convex quadrilateral",
        ));
    }
    let points = data
        .edges
        .iter()
        .map(|e| e.curve.eval(e.curve.range().lo).point)
        .collect::<Vec<_>>();
    let frame = Frame::new(
        points[0],
        (points[1] - points[0]).cross(points[2] - points[0]),
        points[1] - points[0],
        model.precision(),
    )?;
    let p = model.precision();
    for e in &data.edges {
        if e.curve.degree() != 1 || e.curve.points().len() != 2 {
            return Err(Error::unimplemented(
                "cutting sheet needs straight boundary segments",
            ));
        }
    }
    let surface = match data.surface {
        Surface::Plane(_) => BSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            2,
            2,
            vec![points[0], points[3], points[1], points[2]],
            None,
        )?,
        Surface::BSurface(s) => s,
        _ => {
            return Err(Error::unimplemented(
                "cutting sheet is not a planar NURBS/plane patch",
            ));
        }
    };
    for &q in surface.points().iter().chain(points.iter()) {
        let d = dot_delta(q, frame.origin(), frame.z());
        if d.lo() < -p.solver() || d.hi() > p.solver() {
            return Err(Error::unimplemented(
                "tool is not planar within its geometry tolerance",
            ));
        }
    }
    let polygon: Vec<_> = points
        .iter()
        .map(|&q| {
            let l = frame.to_local(q);
            Point2::new(l.x, l.y)
        })
        .collect();
    for i in 0..4 {
        for j in 0..4 {
            if j == i || j == (i + 1) % 4 {
                continue;
            }
            if orient2d(polygon[i], polygon[(i + 1) % 4], polygon[j]) != Some(Sign::Positive) {
                return Err(Error::unimplemented(
                    "cutting sheet must be strictly convex",
                ));
            }
        }
    }
    Ok((Plane::new(frame), polygon, surface))
}
fn inside_tool(points: &[Point3], plane: Plane, polygon: &[Point2], p: &Precision) -> Result<()> {
    // Convex-hull inclusion of the entire target projection proves containment of
    // its true intersection, even when the isoparametric representative is tolerant.
    // Enclose world subtraction and projection, rather than trusting sampled points.
    let frame = plane.frame();
    for &q in points {
        let x = dot_delta(q, frame.origin(), frame.x());
        let y = dot_delta(q, frame.origin(), frame.y());
        for (i, a) in polygon.iter().enumerate() {
            let b = polygon[(i + 1) % polygon.len()];
            let dx = Interval::point(b.x) - Interval::point(a.x);
            let dy = Interval::point(b.y) - Interval::point(a.y);
            let turn = dx * (y - Interval::point(a.y)) - dy * (x - Interval::point(a.x));
            // Require strict containment so trimmed endpoints/corner contacts cannot
            // masquerade as a full boundary-to-boundary cut of the target.
            let length = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2)).sqrt();
            if (turn - Interval::point(p.linear() * length)).sign() != Some(Sign::Positive) {
                return Err(Error::unimplemented(
                    "intersection reaches or leaves the cutting sheet boundary",
                ));
            }
        }
    }
    Ok(())
}
fn certified_fraction(
    s: &BSurface,
    constant_u: bool,
    plane: Plane,
    proposed: Option<f64>,
    p: &Precision,
) -> Result<f64> {
    let (nu, nv) = s.counts();
    let (pu, pv) = s.degrees();
    if (if constant_u { pu } else { pv }) != 1 || (if constant_u { nu } else { nv }) != 2 {
        return Err(Error::unimplemented(
            "cut requires a linear NURBS parameter direction",
        ));
    }
    let (uk, vk) = s.knots();
    let knots = if constant_u { uk } else { vk };
    let range = if constant_u {
        s.ranges().0
    } else {
        s.ranges().1
    };
    if knots[0].to_bits() != range.lo.to_bits() || knots[3].to_bits() != range.hi.to_bits() {
        return Err(Error::unimplemented("cut needs a clamped linear direction"));
    }
    let indices = |j| {
        if constant_u {
            (j, nv + j)
        } else {
            (j * nv, j * nv + 1)
        }
    };
    let n = plane.frame().z();
    let o = plane.frame().origin();
    let (a, b) = indices(0);
    let wa = s.weights().map_or(1.0, |w| w[a]);
    let wb = s.weights().map_or(1.0, |w| w[b]);
    let da = n.dot(s.points()[a] - o) * wa;
    let db = n.dot(s.points()[b] - o) * wb;
    let t = proposed.unwrap_or(da / (da - db));
    if !t.is_finite() || t <= 0.0 || t >= 1.0 {
        return Err(Error::unimplemented(
            "intersection is not strictly inside the ruled parameter interval",
        ));
    }
    let mut orientation = None;
    let mut residual = 0.0f64;
    let mut minimum_slope = f64::INFINITY;
    let mut min_weight = f64::INFINITY;
    let mut max_weight = 0.0f64;
    for j in 0..if constant_u { nv } else { nu } {
        let (a, b) = indices(j);
        let wa = s.weights().map_or(1.0, |w| w[a]);
        let wb = s.weights().map_or(1.0, |w| w[b]);
        let d0 = dot_delta(s.points()[a], o, n);
        let d1 = dot_delta(s.points()[b], o, n);
        let margin = Interval::point(p.linear());
        let sign = if (d0 + margin).sign() == Some(Sign::Negative)
            && (d1 - margin).sign() == Some(Sign::Positive)
        {
            Sign::Positive
        } else if (d0 - margin).sign() == Some(Sign::Positive)
            && (d1 + margin).sign() == Some(Sign::Negative)
        {
            Sign::Negative
        } else {
            return Err(Error::unimplemented(
                "plane crossing is tangent, clustered or not globally transverse",
            ));
        };
        if orientation.is_some_and(|old| old != sign) {
            return Err(Error::unimplemented(
                "intersection is not one complete isoparametric branch",
            ));
        }
        orientation = Some(sign);
        let h = d0 * Interval::point(wa) * (Interval::point(1.0) - Interval::point(t))
            + d1 * Interval::point(wb) * Interval::point(t);
        let w = Interval::point(wa) * (Interval::point(1.0) - Interval::point(t))
            + Interval::point(wb) * Interval::point(t);
        let slope = d1 * Interval::point(wb) - d0 * Interval::point(wa);
        minimum_slope = minimum_slope.min(if sign == Sign::Positive {
            slope.lo()
        } else {
            -slope.hi()
        });
        residual = residual.max(h.lo().abs().max(h.hi().abs()));
        min_weight = min_weight.min(wa.min(wb));
        max_weight = max_weight.max(wa.max(wb));
        let bound = p.default_local() * w.lo();
        // Explicit tolerance fallback: each homogeneous coefficient lies within
        // the plane tube. Positive weights extend this bound to the whole curve.
        if h.lo() < -bound || h.hi() > bound {
            return Err(Error::unimplemented(
                "candidate is not an isoparametric plane intersection within precision",
            ));
        }
    }
    // Plane residual alone is insufficient near tangency. Bound parameter-root
    // displacement by the complete homogeneous slope, then convert it to spatial
    // displacement with a global rational ruling derivative bound.
    let reciprocal = |x: f64| -> Result<Interval> {
        if Interval::point(x).sign() != Some(Sign::Positive) {
            return Err(Error::unimplemented("intersection conditioning unresolved"));
        }
        Interval::new((1. / x).next_down(), (1. / x).next_up())
            .ok_or_else(|| Error::unimplemented("intersection conditioning overflow"))
    };
    let weight_ratio = Interval::point(max_weight) * reciprocal(min_weight)?;
    let root_error = (Interval::point(residual)
        + Interval::point(p.solver()) * Interval::point(max_weight))
        * reciprocal(minimum_slope)?;
    let mut diameter_squared = Interval::point(0.);
    for axis in 0..3 {
        let coordinate = |q: &Point3| [q.x, q.y, q.z][axis];
        let lo = s
            .points()
            .iter()
            .map(coordinate)
            .fold(f64::INFINITY, f64::min);
        let hi = s
            .points()
            .iter()
            .map(coordinate)
            .fold(f64::NEG_INFINITY, f64::max);
        diameter_squared = diameter_squared + (Interval::point(hi) - Interval::point(lo)).sqr();
    }
    let squared_error = diameter_squared * weight_ratio.sqr().sqr() * root_error.sqr();
    if (squared_error - Interval::point(p.default_local()).sqr()).sign() != Some(Sign::Negative) {
        return Err(Error::unimplemented(
            "isoparametric section lacks a spatial error bound",
        ));
    }
    Ok(t)
}

/// Split a full ruled NURBS face by a bounded convex planar sheet. The entire cut
/// must lie strictly inside the tool. Linear-direction transversality proves one
/// branch; a common isoparameter is accepted only with a complete control-hull bound.
/// Polynomial Bezier extrusions also support oblique sections: spatial and UV curves
/// are constructed algebraically, verified over their whole domains, and used to
/// trim two loops on the original carrier. Already-trimmed target faces are refused.
/// An optional backend proposes the parameter; its success/absence is never taken
/// as proof of intersection completeness. General curved tool faces are unsupported.
pub fn body_split_face_by_sheet(
    model: &mut Model,
    body: BodyId,
    face: FaceId,
    tool: FaceId,
    backend: Option<&dyn GeometryBackend>,
) -> Result<FaceSplit> {
    body_split_face_by_sheet_with_options(
        model,
        body,
        face,
        tool,
        backend,
        SheetCutOptions::default(),
    )
}

/// Cut with an explicit capability selection. Rational oblique cuts are rebuilt
/// from the original carrier and independently certified over their full domain;
/// provider coordinates and reported deviations never establish correctness.
/// An empty provider result retains the legacy refusal behavior. Publication is
/// transactional and every resulting body passes the unchanged complete checker.
pub fn body_split_face_by_sheet_with_options(
    model: &mut Model,
    body: BodyId,
    face: FaceId,
    tool: FaceId,
    backend: Option<&dyn GeometryBackend>,
    options: SheetCutOptions,
) -> Result<FaceSplit> {
    body_check_complete(model, body)?;
    let source_faces = model.body_faces(body)?;
    let position = source_faces
        .iter()
        .position(|&f| f == face)
        .ok_or_else(|| Error::invalid_argument("target face is not in body"))?;
    let tool_bodies = model
        .bodies()
        .filter_map(|b| {
            model
                .body_faces(b)
                .ok()
                .filter(|f| f.contains(&tool))
                .map(|_| b)
        })
        .collect::<Vec<_>>();
    if tool_bodies.is_empty() {
        return Err(Error::invalid_argument(
            "tool face must belong to a checked body",
        ));
    }
    for b in tool_bodies {
        body_check_complete(model, b)?;
    }
    let source = boundary::extract(model, face)?;
    let Surface::BSurface(s) = source.surface else {
        return Err(Error::unimplemented(
            "bounded-sheet cutting currently needs a NURBS target",
        ));
    };
    let (plane, polygon, tool_surface) = planar_tool(model, tool)?;
    let p = *model.precision();
    let candidates = if let Some(backend) = backend {
        if backend.precision() != p {
            return Err(Error::invalid_argument(
                "backend and model precisions differ",
            ));
        }
        Some(backend.intersect_surfaces(&s, &tool_surface)?)
    } else {
        None
    };
    let mut solution = None;
    for constant_u in [true, false] {
        let r = if constant_u {
            s.ranges().0
        } else {
            s.ranges().1
        };
        let proposals: Vec<Option<f64>> = if let Some(candidates) = &candidates {
            candidates
                .iter()
                .map(|c| {
                    let d = c.first_uv.range();
                    let uv = c.first_uv.eval((d.lo + d.hi) * 0.5).point;
                    Some(((if constant_u { uv.x } else { uv.y }) - r.lo) / (r.hi - r.lo))
                })
                .collect()
        } else {
            vec![None]
        };
        for proposed in proposals {
            let Ok(t) = certified_fraction(&s, constant_u, plane, proposed, &p) else {
                continue;
            };
            inside_tool(s.points(), plane, &polygon, &p)?;
            solution = Some((constant_u, t));
            break;
        }
        if solution.is_some() {
            break;
        }
    }
    let mut crossing = None;
    if solution.is_none() && candidates.as_ref().is_none_or(|c| !c.is_empty()) {
        for linear_u in [true, false] {
            let section = extruded::section(&s, linear_u, plane, &p).or_else(|error| {
                if options.allow_rational_oblique {
                    rational::section(&s, &tool_surface, linear_u, plane, &p)
                } else {
                    Err(error)
                }
            });
            if let Ok((curve, fractions)) = section {
                inside_tool(s.points(), plane, &polygon, &p)?;
                crossing = Some((linear_u, curve, fractions));
                break;
            }
        }
    }
    model.transaction(|m| {
        let result = if let Some((axis, t)) = solution {
            body_split_face_isoparametric(m, body, face, axis, t)?
        } else if let Some((linear_u, curve, fractions)) = crossing {
            crate::edit::split_rectangle_crossing(m, body, face, linear_u, curve, fractions)?
        } else {
            return Err(Error::unimplemented("no certified bounded ruled-face cut; numerical candidates do not establish topology"));
        };
        let new = m.body_faces(result)?;
        let faces = [new[position], new[position + 1]];
        let face_edges = |f| -> Result<Vec<EdgeId>> {
            let mut out = Vec::new();
            for &lp in &m.face(f)?.loops {
                for &fin in &m.loop_(lp)?.fins {
                    out.push(m.fin(fin)?.edge);
                }
            }
            Ok(out)
        };
        let a = face_edges(faces[0])?;
        let b = face_edges(faces[1])?;
        let shared = a.into_iter().filter(|e| b.contains(e)).collect::<Vec<_>>();
        let [intersection] = shared.as_slice() else {
            return Err(Error::internal(
                "split faces do not share exactly one intersection edge",
            ));
        };
        Ok(FaceSplit {
            body: result,
            faces,
            intersection: *intersection,
        })
    })
}
