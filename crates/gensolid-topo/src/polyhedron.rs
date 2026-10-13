//! Construction of exact planar solids from an oriented, convex polygon boundary.

use std::collections::{BTreeMap, HashSet};

use gensolid_core::{Error, Result};
use gensolid_geom::{Curve, Line, Plane, Surface};
use gensolid_math::{Frame, Interval, Point2, Point3, Sign, predicates::orient2d};

use crate::entity::*;
use crate::model::{Capacity, Model};

fn point(points: &[Point3], index: usize) -> Result<Point3> {
    points
        .get(index)
        .copied()
        .ok_or_else(|| Error::invalid_argument("polyhedron face index is outside the vertex array"))
}

pub(crate) fn vector(point: Point3, origin: Point3) -> [Interval; 3] {
    [
        Interval::point(point.x) - Interval::point(origin.x),
        Interval::point(point.y) - Interval::point(origin.y),
        Interval::point(point.z) - Interval::point(origin.z),
    ]
}

pub(crate) fn dot(a: [Interval; 3], b: [Interval; 3]) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

struct PreparedFace {
    corners: Vec<usize>,
    plane: Plane,
}

/// Builds a connected convex polyhedron with planar convex polygon faces.
///
/// `faces` contains outward, counter-clockwise vertex indices. Every vertex must be
/// used, every edge must occur in two opposite directions, and no face may have holes.
/// Inputs and all carrier geometry are validated before the first mutation. This is
/// a topology-layer constructor: a higher layer must run `gensolid_check::check_body`
/// before handing the body to an application.
pub fn make_convex_polyhedron(
    model: &mut Model,
    points: &[Point3],
    faces: &[Vec<usize>],
) -> Result<BodyId> {
    let precision = *model.precision();
    if points.len() < 4 || faces.len() < 4 {
        return Err(Error::invalid_argument(
            "a solid polyhedron needs at least four vertices and faces",
        ));
    }
    for p in points {
        if !p.is_finite()
            || [p.x, p.y, p.z]
                .into_iter()
                .any(|coordinate| coordinate.abs() > precision.size() * 0.5)
        {
            return Err(Error::invalid_argument(
                "polyhedron vertex is non-finite or outside the session size box",
            ));
        }
    }
    for (index, a) in points.iter().enumerate() {
        for b in &points[index + 1..] {
            let d = vector(*a, *b);
            let separation = dot(d, d) - Interval::point(precision.linear()).sqr();
            match separation.sign() {
                Some(Sign::Positive) => (),
                Some(_) => {
                    return Err(Error::degenerate(
                        "polyhedron has vertices coincident at session precision",
                    ));
                }
                None => {
                    return Err(Error::unimplemented(
                        "polyhedron vertex separation is unresolved at session precision",
                    ));
                }
            }
        }
    }
    let mut used = HashSet::new();
    let mut edges: BTreeMap<(usize, usize), Vec<(usize, bool)>> = BTreeMap::new();
    let mut prepared = Vec::new();
    let mut fin_count = 0usize;
    let uncertainty = Interval::new(-precision.linear(), precision.linear())
        .ok_or_else(|| Error::internal("invalid session precision interval"))?;
    for (face_index, corners) in faces.iter().enumerate() {
        if corners.len() < 3 {
            return Err(Error::invalid_argument(
                "polyhedron face needs at least three vertices",
            ));
        }
        let mut unique = HashSet::new();
        for index in corners {
            point(points, *index)?;
            if !unique.insert(*index) {
                return Err(Error::invalid_argument("polyhedron face repeats a vertex"));
            }
            used.insert(*index);
        }
        fin_count = fin_count
            .checked_add(corners.len())
            .ok_or_else(|| Error::invalid_argument("polyhedron fin count overflows"))?;
        let origin = point(points, corners[0])?;
        let first_edge = point(points, corners[1])? - origin;
        let first_direction = first_edge
            .normalized(precision.linear())
            .ok_or_else(|| Error::degenerate("polyhedron edge is degenerate"))?;
        let mut normal = None;
        for index in &corners[2..] {
            let Some(direction) = (point(points, *index)? - origin).normalized(precision.linear())
            else {
                continue;
            };
            if let Some(candidate) = first_direction
                .cross(direction)
                .normalized(precision.angular())
            {
                normal = Some(candidate);
                break;
            }
        }
        let normal = normal.ok_or_else(|| {
            Error::degenerate("polyhedron face is collinear at angular precision")
        })?;
        let frame = Frame::new(origin, normal, first_direction, &precision)?;
        let polygon = corners
            .iter()
            .map(|index| {
                let local = frame.to_local(point(points, *index)?);
                Ok(Point2::new(local.x, local.y))
            })
            .collect::<Result<Vec<_>>>()?;
        let n = [
            Interval::point(frame.z().x),
            Interval::point(frame.z().y),
            Interval::point(frame.z().z),
        ];
        for index in corners {
            let distance = dot(n, vector(point(points, *index)?, origin));
            match (distance + uncertainty).sign() {
                None => {
                    // Explicit tolerance fallback for uncertain coplanarity.
                    if distance.lo() < -precision.linear() || distance.hi() > precision.linear() {
                        return Err(Error::unimplemented(
                            "polyhedron face planarity is unresolved at session precision",
                        ));
                    }
                }
                _ => {
                    return Err(Error::invalid_argument(
                        "polyhedron face is not planar at session precision",
                    ));
                }
            }
        }
        for (k, a) in polygon.iter().enumerate() {
            let b = polygon[(k + 1) % polygon.len()];
            for p in &polygon {
                match orient2d(*a, b, *p) {
                    Some(Sign::Positive | Sign::Zero) => (),
                    Some(Sign::Negative) => {
                        return Err(Error::invalid_argument(
                            "polyhedron face must be a convex counter-clockwise polygon",
                        ));
                    }
                    None => {
                        return Err(Error::degenerate(
                            "polyhedron face projection is not finite",
                        ));
                    }
                }
            }
            let start = corners[k];
            let end = corners[(k + 1) % corners.len()];
            edges
                .entry((start.min(end), start.max(end)))
                .or_default()
                .push((face_index, start < end));
        }
        for p in points {
            let distance = dot(n, vector(*p, origin));
            match (distance + uncertainty).sign() {
                Some(Sign::Negative | Sign::Zero) => (),
                Some(Sign::Positive) => {
                    return Err(Error::invalid_argument(
                        "polyhedron is nonconvex or a face points inward",
                    ));
                }
                None => {
                    // Explicit tolerance fallback for points on a rounded boundary.
                    if distance.hi() > precision.linear() {
                        return Err(Error::unimplemented(
                            "polyhedron convexity is unresolved at session precision",
                        ));
                    }
                }
            }
        }
        prepared.push(PreparedFace {
            corners: corners.clone(),
            plane: Plane::new(frame),
        });
    }
    if used.len() != points.len() {
        return Err(Error::invalid_argument("polyhedron has unused vertices"));
    }
    let mut adjacency = vec![Vec::new(); faces.len()];
    let mut carriers = BTreeMap::new();
    for (key, uses) in &edges {
        let [(a, sense_a), (b, sense_b)] = uses.as_slice() else {
            return Err(Error::invalid_argument(
                "polyhedron edge must have exactly two face uses",
            ));
        };
        if a == b || sense_a == sense_b {
            return Err(Error::invalid_argument(
                "polyhedron edge uses must be in opposite directions",
            ));
        }
        adjacency[*a].push(*b);
        adjacency[*b].push(*a);
        carriers.insert(
            *key,
            Line::through(point(points, key.0)?, point(points, key.1)?, &precision)?,
        );
    }
    let mut seen = HashSet::new();
    let mut pending = vec![0usize];
    while let Some(face) = pending.pop() {
        if seen.insert(face) {
            pending.extend(&adjacency[face]);
        }
    }
    if seen.len() != faces.len() {
        return Err(Error::invalid_argument(
            "polyhedron boundary is disconnected",
        ));
    }
    let reference = point(points, 0)?;
    let mut volume = Interval::point(0.0);
    for face in &prepared {
        let origin = point(points, face.corners[0])?;
        for pair in face.corners[1..].windows(2) {
            volume = volume
                + dot(
                    vector(origin, reference),
                    cross(
                        vector(point(points, pair[0])?, origin),
                        vector(point(points, pair[1])?, origin),
                    ),
                );
        }
    }
    if volume.sign() != Some(Sign::Positive) {
        return Err(Error::degenerate(
            "polyhedron does not have a certified positive enclosed volume",
        ));
    }
    model.ensure_capacity(&Capacity {
        bodies: 1,
        regions: 2,
        shells: 2,
        faces: faces.len(),
        loops: faces.len(),
        fins: fin_count,
        edges: edges.len(),
        vertices: points.len(),
        curves: edges.len(),
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
    let material = raw.insert_region(Region {
        body,
        solid: true,
        shells: Vec::new(),
    })?;
    raw.body_mut(body)?.regions = vec![empty, material];
    let empty_shell = |region| Shell {
        region,
        faces: Vec::new(),
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    };
    let front = raw.insert_shell(empty_shell(empty))?;
    let back = raw.insert_shell(empty_shell(material))?;
    raw.region_mut(empty)?.shells.push(front);
    raw.region_mut(material)?.shells.push(back);
    let vertices = points
        .iter()
        .map(|p| {
            raw.insert_vertex(Vertex {
                point: *p,
                tolerance: None,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut edge_ids = BTreeMap::new();
    for (key, line) in carriers {
        let curve = raw.geom_mut().add_curve(Curve::Line(line))?;
        let start = vertices
            .get(key.0)
            .copied()
            .ok_or_else(|| Error::internal("prepared polyhedron start vertex is missing"))?;
        let end = vertices
            .get(key.1)
            .copied()
            .ok_or_else(|| Error::internal("prepared polyhedron end vertex is missing"))?;
        edge_ids.insert(
            key,
            raw.insert_edge(Edge {
                fins: Vec::new(),
                start: Some(start),
                end: Some(end),
                curve: Some(curve),
                tolerance: None,
            })?,
        );
    }
    for prepared in prepared {
        let surface = raw.geom_mut().add_surface(Surface::Plane(prepared.plane))?;
        let face = raw.insert_face(Face {
            front: Some(front),
            back: Some(back),
            loops: Vec::new(),
            surface: Some(surface),
            sense: Sense::Positive,
        })?;
        let lp = raw.insert_loop(Loop {
            face,
            fins: Vec::new(),
        })?;
        for (k, start) in prepared.corners.iter().enumerate() {
            let end = prepared.corners[(k + 1) % prepared.corners.len()];
            let edge = *edge_ids
                .get(&((*start).min(end), (*start).max(end)))
                .ok_or_else(|| Error::internal("prepared polyhedron edge is missing"))?;
            let fin = raw.insert_fin(Fin {
                owner: lp,
                edge,
                sense: if *start < end {
                    Sense::Positive
                } else {
                    Sense::Negative
                },
            })?;
            raw.edge_mut(edge)?.fins.push(fin);
            raw.loop_mut(lp)?.fins.push(fin);
        }
        raw.face_mut(face)?.loops.push(lp);
        raw.shell_mut(front)?.faces.push(face);
        raw.shell_mut(back)?.faces.push(face);
    }
    Ok(body)
}
