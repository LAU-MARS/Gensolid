//! Exact bounded NURBS faces and shared topology for revolution and multi-section loft.
use crate::{
    entity::*,
    model::{Capacity, Model},
};
use gensolid_core::{Error, Result};
use gensolid_geom::{BCurve, BSurface, Curve, Line, Plane, Surface};
use gensolid_math::{Frame, Interval, Point2, Point3, Sign, Vec3, predicates::orient2d};
use std::collections::BTreeMap;

/// Constructed solid and the bounded parameter faces for inspection and display.
#[derive(Clone, Debug)]
pub struct ParametricBody {
    /// B-rep body. A checked modelling layer must validate it before returning it.
    pub body: BodyId,
    /// Rectangular NURBS face domains. Planar caps are represented by their loops.
    pub patches: Vec<FaceId>,
}
struct FaceData {
    corners: Vec<usize>,
    surface: Surface,
    curves: Vec<Curve>,
    sense: Sense,
}
fn assemble(model: &mut Model, points: &[Point3], faces: Vec<FaceData>) -> Result<ParametricBody> {
    let mut uses: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
    for face in &faces {
        if face.corners.len() != face.curves.len() {
            return Err(Error::internal("boundary curve count mismatch"));
        }
        for (k, &a) in face.corners.iter().enumerate() {
            let b = face.corners[(k + 1) % face.corners.len()];
            if a == b || a >= points.len() || b >= points.len() {
                return Err(Error::degenerate("collapsed parametric boundary"));
            }
            uses.entry((a.min(b), a.max(b))).or_default().push((a, b));
        }
    }
    if uses
        .values()
        .any(|v| v.len() != 2 || v[0] != (v[1].1, v[1].0))
    {
        return Err(Error::invalid_argument(
            "parametric boundary is not a closed oriented two-manifold",
        ));
    }
    if points.iter().any(|p| {
        !p.is_finite()
            || [p.x, p.y, p.z]
                .into_iter()
                .any(|x| x.abs() > model.precision().size() * 0.5)
    }) {
        return Err(Error::invalid_argument(
            "parametric body outside session size box",
        ));
    }
    model.ensure_capacity(&Capacity {
        bodies: 1,
        regions: 2,
        shells: 2,
        faces: faces.len(),
        loops: faces.len(),
        fins: faces.iter().map(|f| f.corners.len()).sum(),
        edges: uses.len(),
        vertices: points.len(),
        curves: uses.len(),
        surfaces: faces.len(),
    })?;
    let mut raw = model.raw();
    let body = raw.insert_body(Body {
        kind: BodyKind::Solid,
        regions: Vec::new(),
    })?;
    let empty = raw.insert_region(Region {
        body,
        solid: false,
        shells: Vec::new(),
    })?;
    let solid = raw.insert_region(Region {
        body,
        solid: true,
        shells: Vec::new(),
    })?;
    raw.body_mut(body)?.regions = vec![empty, solid];
    let shell = |region| Shell {
        region,
        faces: Vec::new(),
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    };
    let front = raw.insert_shell(shell(empty))?;
    let back = raw.insert_shell(shell(solid))?;
    raw.region_mut(empty)?.shells.push(front);
    raw.region_mut(solid)?.shells.push(back);
    let vertices = points
        .iter()
        .map(|&point| {
            raw.insert_vertex(Vertex {
                point,
                tolerance: None,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut edges: BTreeMap<(usize, usize), (EdgeId, usize)> = BTreeMap::new();
    let mut patches = Vec::new();
    for data in faces {
        let is_patch = matches!(data.surface, Surface::BSurface(_));
        let surface = raw.geom_mut().add_surface(data.surface)?;
        let face = raw.insert_face(Face {
            front: Some(front),
            back: Some(back),
            surface: Some(surface),
            loops: Vec::new(),
            sense: data.sense,
        })?;
        if is_patch {
            patches.push(face);
        }
        let lp = raw.insert_loop(Loop {
            face,
            fins: Vec::new(),
        })?;
        for (i, curve) in data.curves.into_iter().enumerate() {
            let a = data.corners[i];
            let b = data.corners[(i + 1) % data.corners.len()];
            let key = (a.min(b), a.max(b));
            let (edge, start) = if let Some(value) = edges.get(&key) {
                *value
            } else {
                let curve = raw.geom_mut().add_curve(curve)?;
                let edge = raw.insert_edge(Edge {
                    fins: Vec::new(),
                    start: Some(vertices[a]),
                    end: Some(vertices[b]),
                    curve: Some(curve),
                    tolerance: None,
                })?;
                edges.insert(key, (edge, a));
                (edge, a)
            };
            let fin = raw.insert_fin(Fin {
                owner: lp,
                edge,
                sense: if start == a {
                    Sense::Positive
                } else {
                    Sense::Negative
                },
            })?;
            raw.loop_mut(lp)?.fins.push(fin);
            raw.edge_mut(edge)?.fins.push(fin);
        }
        raw.face_mut(face)?.loops.push(lp);
        raw.shell_mut(front)?.faces.push(face);
        raw.shell_mut(back)?.faces.push(face);
    }
    Ok(ParametricBody { body, patches })
}
fn ccw(profile: &[Point2]) -> Result<Vec<Point2>> {
    if profile.len() > 256 {
        return Err(Error::unimplemented(
            "convex profile certification budget exceeded",
        ));
    }
    if profile.len() < 3 || profile.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
        return Err(Error::invalid_argument(
            "profile requires three finite corners",
        ));
    }
    let mut sign = None;
    for i in 0..profile.len() {
        let a = profile[i];
        let b = profile[(i + 1) % profile.len()];
        for (j, &c) in profile.iter().enumerate() {
            if j == i || j == (i + 1) % profile.len() {
                continue;
            }
            let s = orient2d(a, b, c)
                .ok_or_else(|| Error::degenerate("profile orientation unavailable"))?;
            if s == Sign::Zero {
                return Err(Error::degenerate(
                    "profile must have strictly convex corners",
                ));
            }
            if sign.is_some_and(|old| old != s) {
                return Err(Error::invalid_argument("profile must be simple and convex"));
            }
            sign = Some(s);
        }
    }
    let mut p = profile.to_vec();
    if sign == Some(Sign::Negative) {
        p.reverse();
    }
    Ok(p)
}
fn lines(points: &[Point3], corners: &[usize], p: &gensolid_core::Precision) -> Result<Vec<Curve>> {
    corners
        .iter()
        .enumerate()
        .map(|(i, &a)| {
            Line::through(points[a], points[corners[(i + 1) % corners.len()]], p).map(Curve::Line)
        })
        .collect()
}
fn cap(
    points: &[Point3],
    corners: Vec<usize>,
    frame: Frame,
    p: &gensolid_core::Precision,
) -> Result<FaceData> {
    Ok(FaceData {
        curves: lines(points, &corners, p)?,
        corners,
        surface: Surface::Plane(Plane::new(frame)),
        sense: Sense::Positive,
    })
}
/// Revolve a convex (radius, axial-height) profile about frame Z by a positive angle.
/// Rational quadratic arcs retain circular geometry exactly; no mesh is used.
/// Profiles touching/crossing the axis are explicitly unsupported pending pole topology.
pub fn make_revolution(
    model: &mut Model,
    profile: &[Point2],
    angle: f64,
    frame: &Frame,
) -> Result<ParametricBody> {
    let p = *model.precision();
    let profile = ccw(profile)?;
    if !angle.is_finite() || angle <= p.angular() || angle > std::f64::consts::TAU + p.angular() {
        return Err(Error::invalid_argument(
            "revolution angle must lie in (angular precision, full turn]",
        ));
    }
    for point in &profile {
        if (Interval::point(point.x) - Interval::point(p.linear())).sign() != Some(Sign::Positive) {
            return Err(Error::unimplemented(
                "axis-touching revolution needs pole topology",
            ));
        }
    }
    let full = (angle - std::f64::consts::TAU).abs() <= p.angular();
    let angle = if full { std::f64::consts::TAU } else { angle };
    let segments = (angle / std::f64::consts::FRAC_PI_2).ceil() as usize;
    let rings = if full { segments } else { segments + 1 };
    let n = profile.len();
    let mut points = Vec::new();
    for j in 0..rings {
        let a = angle * j as f64 / segments as f64;
        for q in &profile {
            points.push(frame.to_world(Vec3::new(q.x * a.cos(), q.x * a.sin(), q.y)));
        }
    }
    let mut faces = Vec::new();
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    for j in 0..segments {
        let next = (j + 1) % rings;
        let a = angle * j as f64 / segments as f64;
        let b = angle * (j + 1) as f64 / segments as f64;
        let mid = (a + b) * 0.5;
        let w = ((b - a) * 0.5).cos();
        let arc_points = |i: usize| {
            vec![
                points[j * n + i],
                frame.to_world(Vec3::new(
                    profile[i].x * mid.cos() / w,
                    profile[i].x * mid.sin() / w,
                    profile[i].y,
                )),
                points[next * n + i],
            ]
        };
        for i in 0..n {
            let k = (i + 1) % n;
            let aa = arc_points(i);
            let bb = arc_points(k);
            let net = (0..3).flat_map(|s| [aa[s], bb[s]]).collect();
            let surface = BSurface::new(
                2,
                1,
                knots.clone(),
                vec![0.0, 0.0, 1.0, 1.0],
                3,
                2,
                net,
                Some(vec![1.0, 1.0, w, w, 1.0, 1.0]),
            )?;
            let mut reverse = bb;
            reverse.reverse();
            let corners = vec![j * n + i, next * n + i, next * n + k, j * n + k];
            let curves = vec![
                Curve::BCurve(BCurve::new(2, knots.clone(), aa, Some(vec![1.0, w, 1.0]))?),
                Curve::Line(Line::through(
                    points[next * n + i],
                    points[next * n + k],
                    &p,
                )?),
                Curve::BCurve(BCurve::new(
                    2,
                    knots.clone(),
                    reverse,
                    Some(vec![1.0, w, 1.0]),
                )?),
                Curve::Line(Line::through(points[j * n + k], points[j * n + i], &p)?),
            ];
            faces.push(FaceData {
                corners,
                surface: Surface::BSurface(surface),
                curves,
                sense: Sense::Positive,
            });
        }
    }
    if !full {
        let start = Frame::new(frame.origin(), -frame.y(), frame.x(), &p)?;
        let radial = frame.x() * angle.cos() + frame.y() * angle.sin();
        let tangent = frame.y() * angle.cos() - frame.x() * angle.sin();
        let end = Frame::new(frame.origin(), tangent, radial, &p)?;
        faces.push(cap(&points, (0..n).collect(), start, &p)?);
        faces.push(cap(
            &points,
            (0..n).rev().map(|i| (rings - 1) * n + i).collect(),
            end,
            &p,
        )?);
    }
    assemble(model, &points, faces)
}

/// A polygon section in a plane parallel to the loft frame XY plane.
#[derive(Clone, Debug)]
pub struct LoftSection {
    /// Strictly increasing local axial coordinate.
    pub height: f64,
    /// Corresponding convex polygon corners, with common count and orientation.
    pub profile: Vec<Point2>,
}
/// Multi-section ruled loft with independently shaped, offset and rotated profiles.
/// Corresponding corners define exact bilinear NURBS faces. Each interpolated profile
/// is interval-certified convex across the whole slab, preventing interior foldovers.
pub fn make_loft_sections(
    model: &mut Model,
    sections: &[LoftSection],
    frame: &Frame,
) -> Result<ParametricBody> {
    let p = *model.precision();
    if sections.len() < 2 || sections.len() > 64 {
        return Err(Error::invalid_argument("loft needs 2..64 sections"));
    }
    let profiles = sections
        .iter()
        .map(|s| ccw(&s.profile))
        .collect::<Result<Vec<_>>>()?;
    let n = profiles[0].len();
    if profiles.iter().any(|v| v.len() != n) {
        return Err(Error::invalid_argument(
            "loft sections require explicit one-to-one corner correspondence",
        ));
    }
    for pair in sections.windows(2) {
        if (Interval::point(pair[1].height)
            - Interval::point(pair[0].height)
            - Interval::point(p.linear()))
        .sign()
            != Some(Sign::Positive)
        {
            return Err(Error::invalid_argument(
                "loft heights must increase by more than linear precision",
            ));
        }
    }
    for pair in profiles.windows(2) {
        certify_slab(&pair[0], &pair[1], &p)?;
    }
    let points = sections
        .iter()
        .zip(&profiles)
        .flat_map(|(s, profile)| {
            profile
                .iter()
                .map(move |q| frame.to_world(Vec3::new(q.x, q.y, s.height)))
        })
        .collect::<Vec<_>>();
    let mut faces = Vec::new();
    let knots = vec![0.0, 0.0, 1.0, 1.0];
    for j in 0..sections.len() - 1 {
        for i in 0..n {
            let k = (i + 1) % n;
            let a = j * n + i;
            let b = j * n + k;
            let c = (j + 1) * n + k;
            let d = (j + 1) * n + i;
            let corners = vec![a, b, c, d];
            let surface = BSurface::new(
                1,
                1,
                knots.clone(),
                knots.clone(),
                2,
                2,
                vec![points[a], points[d], points[b], points[c]],
                None,
            )?;
            faces.push(FaceData {
                curves: lines(&points, &corners, &p)?,
                corners,
                surface: Surface::BSurface(surface),
                sense: Sense::Positive,
            });
        }
    }
    let first = sections[0].height;
    let last = sections[sections.len() - 1].height;
    faces.push(cap(
        &points,
        (0..n).rev().collect(),
        Frame::new(
            frame.to_world(Vec3::new(0.0, 0.0, first)),
            -frame.z(),
            frame.x(),
            &p,
        )?,
        &p,
    )?);
    faces.push(cap(
        &points,
        (0..n).map(|i| (sections.len() - 1) * n + i).collect(),
        Frame::new(
            frame.to_world(Vec3::new(0.0, 0.0, last)),
            frame.z(),
            frame.x(),
            &p,
        )?,
        &p,
    )?);
    assemble(model, &points, faces)
}
fn certify_slab(a: &[Point2], b: &[Point2], p: &gensolid_core::Precision) -> Result<()> {
    let mut pending = vec![(0.0, 1.0, 0usize)];
    let mut visited = 0;
    while let Some((lo, hi, depth)) = pending.pop() {
        visited += 1;
        if visited > 4096 {
            return Err(Error::unimplemented(
                "loft convexity certification budget exceeded",
            ));
        }
        let t = Interval::new(lo, hi).ok_or_else(|| Error::internal("invalid slab interval"))?;
        let points = a
            .iter()
            .zip(b)
            .map(|(a, b)| {
                [
                    Interval::point(a.x) + (Interval::point(b.x) - Interval::point(a.x)) * t,
                    Interval::point(a.y) + (Interval::point(b.y) - Interval::point(a.y)) * t,
                ]
            })
            .collect::<Vec<_>>();
        let mut certified = true;
        for i in 0..a.len() {
            let j = (i + 1) % a.len();
            for k in 0..a.len() {
                if k == i || k == j {
                    continue;
                }
                let orientation = (points[j][0] - points[i][0]) * (points[k][1] - points[i][1])
                    - (points[j][1] - points[i][1]) * (points[k][0] - points[i][0]);
                match orientation.sign() {
                    Some(Sign::Positive) => (),
                    Some(Sign::Negative | Sign::Zero) => {
                        return Err(Error::degenerate(
                            "loft interpolated profile folds or collapses",
                        ));
                    }
                    None => certified = false,
                }
            }
        }
        if !certified {
            if depth == 32 || hi - lo <= p.angular() {
                return Err(Error::unimplemented(
                    "loft slab convexity cannot be certified",
                ));
            }
            let mid = (lo + hi) * 0.5;
            pending.push((lo, mid, depth + 1));
            pending.push((mid, hi, depth + 1));
        }
    }
    Ok(())
}
