//! Conservative certification of planar BREP_WITH_VOIDS shell relationships.
//!
//! The immutable body checker checks manifold uses and geometric carriers. It does
//! not establish shell containment, so unsupported relationships are refused here.
//! No mesh is constructed; all decisions use plane halfspaces and exact predicates.

use std::collections::{HashMap, HashSet};

use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{Curve, Surface};
use gensolid_math::{Interval, Point2, Point3, Sign, Vec3, predicates::orient2d};
use gensolid_topo::{BodyId, EdgeId, FaceId, Model, Sense, ShellId};

struct Halfspace {
    origin: Point3,
    normal: Vec3,
}

struct Polyhedron {
    planes: Vec<Halfspace>,
    vertices: Vec<Point3>,
}

fn coordinate(value: f64, reference: f64) -> Interval {
    Interval::point(value) - Interval::point(reference)
}

fn point_vector(point: Point3, reference: Point3) -> [Interval; 3] {
    [
        coordinate(point.x, reference.x),
        coordinate(point.y, reference.y),
        coordinate(point.z, reference.z),
    ]
}

fn dot(a: [Interval; 3], b: [Interval; 3]) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

impl Halfspace {
    fn side(&self, point: Point3) -> Interval {
        dot(
            point_vector(point, self.origin),
            [
                Interval::point(self.normal.x),
                Interval::point(self.normal.y),
                Interval::point(self.normal.z),
            ],
        )
    }
}

/// Called only after the immutable checker established all carrier/manifold checks.
pub(crate) fn validate(model: &Model, body: BodyId) -> Result<()> {
    let regions = &model.body(body)?.regions;
    if regions.len() <= 2 {
        return Ok(());
    }
    let mut boundaries = None;
    for id in regions {
        let region = model.region(*id)?;
        if region.solid {
            if boundaries.is_some() {
                return Err(Error::unimplemented(
                    "STEP void certification needs one connected material region",
                ));
            }
            boundaries = Some(region.shells.as_slice());
        }
    }
    let boundaries = boundaries
        .ok_or_else(|| Error::invalid_argument("STEP void body has no material region"))?;
    let precision = model.precision();
    let outer_shell = boundaries
        .first()
        .ok_or_else(|| Error::internal("void body has no outer shell"))?;
    let outer = polyhedron(model, *outer_shell, false)?;
    let cavities = boundaries
        .iter()
        .skip(1)
        .map(|shell| polyhedron(model, *shell, true))
        .collect::<Result<Vec<_>>>()?;
    for cavity in &cavities {
        for point in &cavity.vertices {
            for plane in &outer.planes {
                let side = plane.side(*point);
                // A full linear-precision margin keeps a closed cavity distinct
                // from the outer boundary, including uncertain rounded placements.
                if (side + Interval::point(precision.linear())).sign() == Some(Sign::Negative) {
                    continue;
                }
                if (side - Interval::point(precision.linear())).sign() == Some(Sign::Positive) {
                    return Err(Error::invalid_argument(
                        "STEP cavity lies outside its convex outer shell",
                    ));
                }
                return Err(Error::unimplemented(
                    "STEP cavity containment cannot be certified at session precision",
                ));
            }
        }
    }
    for (index, a) in cavities.iter().enumerate() {
        for b in &cavities[index + 1..] {
            if separated(a, b, precision) || separated(b, a, precision) {
                continue;
            }
            return Err(Error::unimplemented(
                "STEP cavities have no certified separating face plane",
            ));
        }
    }
    Ok(())
}

fn separated(a: &Polyhedron, b: &Polyhedron, precision: &Precision) -> bool {
    a.planes.iter().any(|plane| {
        b.vertices.iter().all(|point| {
            (plane.side(*point) - Interval::point(precision.linear())).sign()
                == Some(Sign::Positive)
        })
    })
}

fn polyhedron(model: &Model, shell: ShellId, cavity: bool) -> Result<Polyhedron> {
    let face_ids = &model.shell(shell)?.faces;
    let mut vertices = Vec::new();
    let mut planes = Vec::new();
    let mut polygons = Vec::new();
    let mut uses: HashMap<EdgeId, Vec<(FaceId, Sense)>> = HashMap::new();
    for face_id in face_ids {
        let face = model.face(*face_id)?;
        if face.back != Some(shell) || face.loops.len() != 1 {
            return Err(Error::unimplemented(
                "STEP void certification needs one polygon per face bounding the material shell",
            ));
        }
        let surface = face
            .surface
            .ok_or_else(|| Error::invalid_argument("STEP void face has no surface"))?;
        let Surface::Plane(plane) = model.geom().surface(surface)? else {
            return Err(Error::unimplemented(
                "STEP void containment for curved faces is not implemented",
            ));
        };
        let loop_id = face
            .loops
            .first()
            .ok_or_else(|| Error::invalid_argument("STEP void face has no loop"))?;
        let mut polygon = Vec::new();
        for fin in &model.loop_(*loop_id)?.fins {
            let fin_value = model.fin(*fin)?;
            uses.entry(fin_value.edge)
                .or_default()
                .push((*face_id, fin_value.sense));
            let edge = model.edge(fin_value.edge)?;
            let curve = edge
                .curve
                .ok_or_else(|| Error::invalid_argument("STEP void edge has no curve"))?;
            if !matches!(model.geom().curve(curve)?, Curve::Line(_)) {
                return Err(Error::unimplemented(
                    "STEP void certification needs straight polygon edges",
                ));
            }
            let (start, _) = model.fin_vertices(*fin)?;
            let start = start.ok_or_else(|| {
                Error::unimplemented("STEP void certification needs polygon vertices")
            })?;
            polygon.push(model.vertex(start)?.point);
        }
        if polygon.len() < 3 {
            return Err(Error::invalid_argument(
                "STEP void polygon has fewer than three vertices",
            ));
        }
        let projected: Vec<_> = polygon
            .iter()
            .map(|point| {
                let local = plane.frame().to_local(*point);
                Point2::new(local.x, local.y)
            })
            .collect();
        let expected = if face.sense == Sense::Positive {
            Sign::Positive
        } else {
            Sign::Negative
        };
        let mut noncollinear = false;
        for (index, a) in projected.iter().enumerate() {
            let b = projected[(index + 1) % projected.len()];
            for point in &projected {
                match orient2d(*a, b, *point) {
                    Some(Sign::Zero) => (),
                    Some(sign) if sign == expected => noncollinear = true,
                    Some(_) => {
                        return Err(Error::unimplemented(
                            "STEP void certification needs convex face loops with consistent normal orientation",
                        ));
                    }
                    None => return Err(Error::degenerate("non-finite STEP polygon projection")),
                }
            }
        }
        if !noncollinear {
            return Err(Error::degenerate("STEP void polygon is collinear"));
        }
        let face_normal = if face.sense == Sense::Positive {
            plane.frame().z()
        } else {
            -plane.frame().z()
        };
        planes.push(Halfspace {
            origin: plane.frame().origin(),
            normal: if cavity { -face_normal } else { face_normal },
        });
        vertices.extend(polygon.iter().copied());
        polygons.push(polygon);
    }
    let mut adjacency: HashMap<FaceId, Vec<FaceId>> = HashMap::new();
    for edge_uses in uses.values() {
        let [(a, sense_a), (b, sense_b)] = edge_uses.as_slice() else {
            return Err(Error::invalid_argument(
                "STEP void boundary shell is not a closed manifold",
            ));
        };
        if a == b || sense_a == sense_b {
            return Err(Error::invalid_argument(
                "STEP void boundary has inconsistent manifold edge uses",
            ));
        }
        adjacency.entry(*a).or_default().push(*b);
        adjacency.entry(*b).or_default().push(*a);
    }
    let first_face = face_ids
        .first()
        .ok_or_else(|| Error::invalid_argument("STEP void shell is empty"))?;
    let mut seen = HashSet::new();
    let mut pending = vec![*first_face];
    while let Some(face) = pending.pop() {
        if seen.insert(face)
            && let Some(neighbors) = adjacency.get(&face)
        {
            pending.extend(neighbors);
        }
    }
    if seen.len() != face_ids.len() {
        return Err(Error::invalid_argument(
            "STEP closed shell has disconnected face components",
        ));
    }
    let precision = model.precision();
    let tolerance = precision.linear();
    let uncertainty = Interval::new(-tolerance, tolerance)
        .ok_or_else(|| Error::internal("invalid session precision interval"))?;
    for plane in &planes {
        for point in &vertices {
            let distance = plane.side(*point);
            // Input vertices and elementary placements are only required to agree
            // at session precision. Include that uncertainty before certifying a side.
            match (distance + uncertainty).sign() {
                Some(Sign::Negative | Sign::Zero) => (),
                Some(Sign::Positive) => {
                    return Err(Error::unimplemented(
                        "STEP void certification needs a convex, consistently oriented shell",
                    ));
                }
                None => {
                    // Explicit tolerance fallback for an interval straddling zero:
                    // tolerate a coincident boundary only within the session envelope.
                    if distance.hi() > tolerance {
                        return Err(Error::unimplemented(
                            "STEP shell convexity is unresolved at session precision",
                        ));
                    }
                }
            }
        }
    }
    let reference = vertices
        .first()
        .copied()
        .ok_or_else(|| Error::invalid_argument("STEP void shell has no vertices"))?;
    let mut signed_six_volume = Interval::point(0.0);
    for polygon in &polygons {
        let origin = polygon
            .first()
            .copied()
            .ok_or_else(|| Error::internal("empty certified polygon"))?;
        for pair in polygon[1..].windows(2) {
            signed_six_volume = signed_six_volume
                + dot(
                    point_vector(origin, reference),
                    cross(point_vector(pair[0], origin), point_vector(pair[1], origin)),
                );
        }
    }
    let expected = if cavity {
        Sign::Negative
    } else {
        Sign::Positive
    };
    match signed_six_volume.sign() {
        Some(sign) if sign == expected => (),
        Some(_) => {
            return Err(Error::invalid_argument(
                "STEP shell volume orientation does not bound a finite solid or cavity",
            ));
        }
        None => {
            return Err(Error::unimplemented(
                "STEP shell volume sign cannot be certified",
            ));
        }
    }
    Ok(Polyhedron { planes, vertices })
}
