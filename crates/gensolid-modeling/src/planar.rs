//! Convex planar CSG, plane offsets and closed hollow solids.
//! Boolean results are convex-cell decompositions with disjoint interiors.
//! Clipping acts on original B-rep planes and edges; there is no tessellation.
mod fuse;
pub use fuse::body_boolean;

use crate::{body_check_complete, boundary, edit::planar_polygon};
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::Plane;
use gensolid_math::{Frame, Interval, Point2, Point3, Sign, Vec3, predicates::orient2d};
use gensolid_topo::{BodyId, BodyKind, Model, Region, Shell};

/// Regularized solid set operation. Lower-dimensional contact does not create a solid.
#[derive(Clone, Copy, Debug)]
pub enum BooleanOperation {
    /// Material common to both arguments.
    Intersection,
    /// Material in the first argument but not the second.
    Difference,
    /// Material in either argument.
    Union,
}
#[derive(Clone)]
struct Polygon {
    points: Vec<Point3>,
    plane: Plane,
}
#[derive(Clone)]
struct Poly {
    faces: Vec<Polygon>,
}
fn dot(a: Vec3, b: Vec3) -> Interval {
    Interval::point(a.x) * Interval::point(b.x)
        + Interval::point(a.y) * Interval::point(b.y)
        + Interval::point(a.z) * Interval::point(b.z)
}
fn distance(plane: Plane, p: Point3) -> Interval {
    let o = plane.frame().origin();
    let n = plane.frame().z();
    (Interval::point(p.x) - Interval::point(o.x)) * Interval::point(n.x)
        + (Interval::point(p.y) - Interval::point(o.y)) * Interval::point(n.y)
        + (Interval::point(p.z) - Interval::point(o.z)) * Interval::point(n.z)
}
fn side(plane: Plane, p: Point3, precision: &Precision) -> Result<Sign> {
    let d = distance(plane, p);
    let t = Interval::point(precision.default_local());
    if (d - t).sign() == Some(Sign::Positive) {
        return Ok(Sign::Positive);
    }
    if (d + t).sign() == Some(Sign::Negative) {
        return Ok(Sign::Negative);
    }
    // Explicit tolerance fallback for a point unresolved relative to the boundary tube.
    if d.lo() >= -precision.default_local() && d.hi() <= precision.default_local() {
        return Ok(Sign::Zero);
    }
    Err(Error::unimplemented(
        "plane classification unresolved at local tolerance",
    ))
}
fn indexed(poly: &Poly, p: &Precision) -> Result<(Vec<Point3>, Vec<Vec<usize>>)> {
    let mut points = Vec::new();
    let mut faces = Vec::new();
    for face in &poly.faces {
        let mut ids = Vec::new();
        for &q in &face.points {
            let mut found = None;
            for (i, &r) in points.iter().enumerate() {
                if boundary::near(q, r, p)? {
                    found = Some(i);
                    break;
                }
            }
            ids.push(match found {
                Some(i) => i,
                None => {
                    points.push(q);
                    points.len() - 1
                }
            });
        }
        faces.push(ids);
    }
    Ok((points, faces))
}
fn build(m: &mut Model, poly: &Poly) -> Result<BodyId> {
    let (p, f) = indexed(poly, m.precision())?;
    let b = gensolid_topo::make_convex_polyhedron(m, &p, &f)?;
    body_check_complete(m, b)?;
    Ok(b)
}
fn read(model: &Model, body: BodyId) -> Result<Poly> {
    body_check_complete(model, body)?;
    if model.body(body)?.kind != BodyKind::Solid {
        return Err(Error::invalid_argument("operation requires a solid"));
    }
    let precision = model.precision();
    for edge in model.body_edges(body)? {
        if model
            .edge(edge)?
            .tolerance
            .is_some_and(|t| t.value() > precision.default_local())
        {
            return Err(Error::unimplemented(
                "planar operations on widened edge tolerances",
            ));
        }
    }
    for vertex in model.body_vertices(body)? {
        if model
            .vertex(vertex)?
            .tolerance
            .is_some_and(|t| t.value() > precision.default_local())
        {
            return Err(Error::unimplemented(
                "planar operations on widened vertex tolerances",
            ));
        }
    }
    if model.body_faces(body)?.len() > 64 {
        return Err(Error::unimplemented("convex planar face budget exceeded"));
    }
    let mut faces = Vec::new();
    for id in model.body_faces(body)? {
        let points = planar_polygon(model, id)?;
        let a = points[0];
        let b = points[1];
        let c = points[2];
        let frame = Frame::new(a, (b - a).cross(c - a), b - a, model.precision())?;
        faces.push(Polygon {
            points,
            plane: Plane::new(frame),
        });
    }
    let poly = Poly { faces };
    // Independent constructor verifies convexity, winding and manifold topology.
    build(&mut Model::new(*model.precision()), &poly)?;
    Ok(poly)
}
fn unique(points: Vec<Point3>, p: &Precision) -> Result<Vec<Point3>> {
    let mut out = Vec::new();
    for q in points {
        if !out
            .iter()
            .map(|&r| boundary::near(q, r, p))
            .collect::<Result<Vec<_>>>()?
            .contains(&true)
        {
            out.push(q);
        }
    }
    Ok(out)
}
fn hull(points: Vec<Point3>, plane: Plane, p: &Precision) -> Result<Vec<Point3>> {
    let points = unique(points, p)?;
    let mut values: Vec<_> = points
        .into_iter()
        .map(|q| {
            let l = plane.frame().to_local(q);
            (Point2::new(l.x, l.y), q)
        })
        .collect();
    values.sort_by(|a, b| a.0.x.total_cmp(&b.0.x).then(a.0.y.total_cmp(&b.0.y)));
    if values.len() < 3 {
        return Err(Error::degenerate("cut polygon collapsed"));
    }
    let mut out = Vec::new();
    // Monotone-chain hull; orientation is an exact predicate of projected coordinates.
    let chain = |indices: Vec<usize>| -> Result<Vec<usize>> {
        let mut result: Vec<usize> = Vec::new();
        for i in indices {
            while result.len() >= 2 {
                let n = result.len();
                match orient2d(
                    values[result[n - 2]].0,
                    values[result[n - 1]].0,
                    values[i].0,
                ) {
                    Some(Sign::Positive) => break,
                    Some(_) => {
                        result.pop();
                    }
                    None => return Err(Error::degenerate("nonfinite cut polygon")),
                }
            }
            result.push(i);
        }
        Ok(result)
    };
    let mut lower = chain((0..values.len()).collect())?;
    let mut upper = chain((0..values.len()).rev().collect())?;
    lower.pop();
    upper.pop();
    lower.extend(upper);
    for i in lower {
        out.push(values[i].1);
    }
    if out.len() < 3 {
        return Err(Error::degenerate("cut polygon has no area"));
    }
    Ok(out)
}
fn clip(poly: &Poly, plane: Plane, p: &Precision) -> Result<Option<Poly>> {
    let signs = poly
        .faces
        .iter()
        .flat_map(|f| f.points.iter())
        .map(|&q| side(plane, q, p))
        .collect::<Result<Vec<_>>>()?;
    if !signs.contains(&Sign::Negative) {
        return Ok(None);
    }
    if !signs.contains(&Sign::Positive) {
        return Ok(Some(poly.clone()));
    }
    let mut faces = Vec::new();
    let mut cap = Vec::new();
    for face in &poly.faces {
        let mut points = Vec::new();
        for (i, &a) in face.points.iter().enumerate() {
            let b = face.points[(i + 1) % face.points.len()];
            let sa = side(plane, a, p)?;
            let sb = side(plane, b, p)?;
            if sa != Sign::Positive {
                points.push(a);
            }
            if sa == Sign::Zero {
                cap.push(a);
            }
            if (sa == Sign::Positive && sb == Sign::Negative)
                || (sa == Sign::Negative && sb == Sign::Positive)
            {
                let denominator = distance(plane, a) - distance(plane, b);
                if denominator.sign().is_none() || denominator.contains(0.0) {
                    return Err(Error::unimplemented("edge-plane denominator unresolved"));
                }
                let n = plane.frame().z();
                let da = n.dot(a - plane.frame().origin());
                let db = n.dot(b - plane.frame().origin());
                let t = da / (da - db);
                let q = a + (b - a) * t;
                if side(plane, q, p)? != Sign::Zero {
                    return Err(Error::unimplemented(
                        "edge-plane root residual exceeds precision",
                    ));
                }
                points.push(q);
                cap.push(q);
            }
        }
        let points = unique(points, p)?;
        if points.len() >= 3 {
            faces.push(Polygon {
                points,
                plane: face.plane,
            });
        }
    }
    faces.push(Polygon {
        points: hull(cap, plane, p)?,
        plane,
    });
    Ok(Some(Poly { faces }))
}
fn opposite(plane: Plane, p: &Precision) -> Result<Plane> {
    Ok(Plane::new(Frame::new(
        plane.frame().origin(),
        -plane.frame().z(),
        plane.frame().x(),
        p,
    )?))
}
fn difference(a: &Poly, b: &Poly, p: &Precision) -> Result<Vec<Poly>> {
    let mut remainder = Some(a.clone());
    let mut out = Vec::new();
    for face in &b.faces {
        let Some(current) = remainder else {
            break;
        };
        if let Some(piece) = clip(&current, opposite(face.plane, p)?, p)? {
            out.push(piece);
        }
        remainder = clip(&current, face.plane, p)?;
    }
    Ok(out)
}

/// Convex planar Boolean as checked convex cells with disjoint interiors.
/// Adjacent output cells retain their interfaces and are separate bodies. An empty
/// vector is a certified empty regularized result for this supported input family.
/// General curved and nonconvex inputs are explicitly rejected.
pub fn body_boolean_cells(
    model: &mut Model,
    first: BodyId,
    second: BodyId,
    operation: BooleanOperation,
) -> Result<Vec<BodyId>> {
    let a = read(model, first)?;
    let b = read(model, second)?;
    let p = *model.precision();
    let cells = match operation {
        BooleanOperation::Difference => difference(&a, &b, &p)?,
        BooleanOperation::Union => {
            let mut out = difference(&a, &b, &p)?;
            out.push(b);
            out
        }
        BooleanOperation::Intersection => {
            let mut result = Some(a);
            for face in b.faces {
                if let Some(poly) = result {
                    result = clip(&poly, face.plane, &p)?;
                } else {
                    break;
                }
            }
            result.into_iter().collect()
        }
    };
    model.transaction(|m| cells.iter().map(|c| build(m, c)).collect())
}

fn offset(poly: &Poly, distance: f64, p: &Precision) -> Result<Poly> {
    offset_distances(poly, &vec![distance; poly.faces.len()], p)
}
fn offset_distances(poly: &Poly, distances: &[f64], p: &Precision) -> Result<Poly> {
    if distances.iter().any(|d| !d.is_finite()) {
        return Err(Error::invalid_argument("offset must be finite"));
    }
    let (points, indices) = indexed(poly, p)?;
    let mut moved = Vec::new();
    for (i, &point) in points.iter().enumerate() {
        let incident: Vec<_> = indices
            .iter()
            .enumerate()
            .filter(|(_, face)| face.contains(&i))
            .map(|(j, _)| (poly.faces[j].plane.frame().z(), distances[j]))
            .collect();
        let [a, b, c] = incident.as_slice() else {
            return Err(Error::unimplemented(
                "offset requires trihedral vertices and unchanged topology",
            ));
        };
        let cross = b.0.cross(c.0);
        let determinant = dot(a.0, cross);
        if determinant.contains(0.0) {
            return Err(Error::unimplemented(
                "offset planes do not have a certified unique intersection",
            ));
        }
        let delta =
            (cross * a.1 + c.0.cross(a.0) * b.1 + a.0.cross(b.0) * c.1) * (1.0 / a.0.dot(cross));
        let q = point + delta;
        for (n, distance) in &incident {
            let residual = dot(q - point, *n) - Interval::point(*distance);
            if residual.lo() < -p.default_local() || residual.hi() > p.default_local() {
                return Err(Error::unimplemented(
                    "offset vertex residual is outside precision",
                ));
            }
        }
        moved.push(q);
    }
    let mut faces = Vec::new();
    for ((old, ids), &distance) in poly.faces.iter().zip(indices).zip(distances) {
        let f = old.plane.frame();
        let plane = Plane::new(Frame::new(f.origin() + f.z() * distance, f.z(), f.x(), p)?);
        let points = ids.into_iter().map(|i| moved[i]).collect();
        faces.push(Polygon { points, plane });
    }
    for face in &faces {
        for &point in &face.points {
            for bound in &faces {
                if side(bound.plane, point, p)? == Sign::Positive {
                    return Err(Error::degenerate(
                        "offset crossed a supporting plane; topology change required",
                    ));
                }
            }
        }
    }
    Ok(Poly { faces })
}
/// Offset all supporting planes of a simple convex planar solid; positive grows it.
/// Topology changes, collapsed faces and non-trihedral vertices are refused.
pub fn body_offset_planar(model: &mut Model, body: BodyId, distance: f64) -> Result<BodyId> {
    let source = read(model, body)?;
    let result = offset(&source, distance, model.precision())?;
    model.transaction(|m| build(m, &result))
}

/// Hollow a simple convex planar solid with a strictly contained inward plane offset.
/// Creates a finite empty cavity region. This closed-shell operation removes no faces.
pub fn body_hollow_planar(model: &mut Model, body: BodyId, thickness: f64) -> Result<BodyId> {
    let p = *model.precision();
    if !thickness.is_finite()
        || (Interval::point(thickness) - Interval::point(p.linear())).sign() != Some(Sign::Positive)
    {
        return Err(Error::degenerate("wall thickness must exceed precision"));
    }
    let outer = read(model, body)?;
    let inner = offset(&outer, -thickness, &p)?;
    let mut scratch = Model::new(p);
    let inner_body = build(&mut scratch, &inner)?;
    for q in inner.faces.iter().flat_map(|f| f.points.iter()) {
        for f in &outer.faces {
            if side(f.plane, *q, &p)? != Sign::Negative {
                return Err(Error::degenerate("inner shell is not strictly contained"));
            }
        }
    }
    let mut faces = model
        .body_faces(body)?
        .into_iter()
        .map(|f| boundary::extract(model, f))
        .collect::<Result<Vec<_>>>()?;
    let outer_count = faces.len();
    for f in scratch.body_faces(inner_body)? {
        faces.push(boundary::reverse(&boundary::extract(&scratch, f)?)?);
    }
    model.transaction(|m| {
        let result = boundary::assemble(m, &faces, true)?;
        let face_ids = m.body_faces(result)?;
        let regions = m.body(result)?.regions.clone();
        let infinite_shell = m.region(regions[0])?.shells[0];
        let outer_solid_shell = m.region(regions[1])?.shells[0];
        let mut raw = m.raw();
        let cavity = raw.insert_region(Region {
            body: result,
            solid: false,
            shells: Vec::new(),
        })?;
        let front = raw.insert_shell(Shell {
            region: cavity,
            faces: face_ids[outer_count..].to_vec(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        })?;
        let back = raw.insert_shell(Shell {
            region: regions[1],
            faces: face_ids[outer_count..].to_vec(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        })?;
        raw.body_mut(result)?.regions.push(cavity);
        raw.region_mut(cavity)?.shells.push(front);
        raw.region_mut(regions[1])?.shells.push(back);
        raw.shell_mut(infinite_shell)?.faces.truncate(outer_count);
        raw.shell_mut(outer_solid_shell)?
            .faces
            .truncate(outer_count);
        for &f in &face_ids[outer_count..] {
            let face = raw.face_mut(f)?;
            face.front = Some(front);
            face.back = Some(back);
        }
        body_check_complete(m, result)?;
        Ok(result)
    })
}

/// Shell a simple convex planar solid by removing one face and offsetting the others
/// inward. The inner opening stays on the removed face's original plane; planar rim
/// faces connect the two boundary loops. Thickness/topology collapse is rejected.
pub fn body_shell_planar(
    model: &mut Model,
    body: BodyId,
    remove: gensolid_topo::FaceId,
    thickness: f64,
) -> Result<BodyId> {
    let p = *model.precision();
    if !thickness.is_finite()
        || (Interval::point(thickness) - Interval::point(p.linear())).sign() != Some(Sign::Positive)
    {
        return Err(Error::degenerate("wall thickness must exceed precision"));
    }
    let ids = model.body_faces(body)?;
    let index = ids
        .iter()
        .position(|&f| f == remove)
        .ok_or_else(|| Error::invalid_argument("removed face must belong to body"))?;
    let outer = read(model, body)?;
    let mut distances = vec![-thickness; outer.faces.len()];
    distances[index] = 0.0;
    let inner = offset_distances(&outer, &distances, &p)?;
    let mut scratch = Model::new(p);
    let inner_id = build(&mut scratch, &inner)?;
    for &q in inner.faces.iter().flat_map(|f| f.points.iter()) {
        for (i, f) in outer.faces.iter().enumerate() {
            let sign = side(f.plane, q, &p)?;
            if sign == Sign::Positive || (i != index && sign != Sign::Negative) {
                return Err(Error::degenerate("offset shell is not contained"));
            }
        }
    }
    let inner_ids = scratch.body_faces(inner_id)?;
    let mut faces = Vec::new();
    for (i, &id) in ids.iter().enumerate() {
        if i != index {
            faces.push(boundary::extract(model, id)?);
        }
    }
    for (i, &id) in inner_ids.iter().enumerate() {
        if i != index {
            faces.push(boundary::reverse(&boundary::extract(&scratch, id)?)?);
        }
    }
    let a = &outer.faces[index].points;
    let b = &inner.faces[index].points;
    if a.len() != b.len() {
        return Err(Error::internal("offset changed face incidence"));
    }
    for i in 0..a.len() {
        let j = (i + 1) % a.len();
        let points = [a[i], a[j], b[j], b[i]];
        let plane = outer.faces[index].plane;
        let mut edges = Vec::new();
        for k in 0..4 {
            let q = points[k];
            let r = points[(k + 1) % 4];
            let u = plane.frame().to_local(q);
            let v = plane.frame().to_local(r);
            edges.push(boundary::BoundaryEdge {
                curve: boundary::line(q, r)?,
                uv: Some(boundary::line(
                    Point3::new(u.x, u.y, 0.0),
                    Point3::new(v.x, v.y, 0.0),
                )?),
            });
        }
        faces.push(boundary::BoundaryFace {
            surface: gensolid_geom::Surface::Plane(plane),
            sense: gensolid_topo::Sense::Positive,
            edges,
        });
    }
    model.transaction(|m| {
        let b = boundary::assemble(m, &faces, true)?;
        body_check_complete(m, b)?;
        Ok(b)
    })
}
