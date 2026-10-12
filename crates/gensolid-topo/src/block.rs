//! Solid block primitive.

use std::collections::HashMap;

use gensolid_core::{Error, Result};
use gensolid_geom::{Curve, Line, Plane, Surface};
use gensolid_math::{Frame, Vec3};

use crate::entity::*;
use crate::model::Model;

/// Corner `i` has local coordinates `(i & 1, i & 2, i & 4)` scaled to the block.
/// Each face lists its corners anticlockwise about its outward normal.
const FACES: [(Vec3, [usize; 4]); 6] = [
    (Vec3::new(0.0, 0.0, -1.0), [0, 2, 3, 1]),
    (Vec3::new(0.0, 0.0, 1.0), [4, 5, 7, 6]),
    (Vec3::new(0.0, -1.0, 0.0), [0, 1, 5, 4]),
    (Vec3::new(0.0, 1.0, 0.0), [2, 6, 7, 3]),
    (Vec3::new(-1.0, 0.0, 0.0), [0, 4, 6, 2]),
    (Vec3::new(1.0, 0.0, 0.0), [1, 3, 7, 5]),
];

/// Creates a solid box of size `x × y × z`, centred on the origin of `frame` and aligned
/// with its axes.
pub fn make_block(model: &mut Model, x: f64, y: f64, z: f64, frame: &Frame) -> Result<BodyId> {
    let precision = *model.precision();
    for (name, value) in [("x", x), ("y", y), ("z", z)] {
        if !(value.is_finite() && value > precision.linear() && value < precision.size()) {
            return Err(Error::invalid_argument(format!(
                "block {name} must lie in ({}, {}), got {value}",
                precision.linear(),
                precision.size()
            )));
        }
    }
    let corner = |i: usize| {
        let pick = |bit: usize, size: f64| {
            if i & bit == 0 {
                -size / 2.0
            } else {
                size / 2.0
            }
        };
        frame.to_world(Vec3::new(pick(1, x), pick(2, y), pick(4, z)))
    };

    let mut raw = model.raw();
    let body = raw.insert_body(Body {
        kind: BodyKind::Solid,
        regions: Vec::new(),
    })?;
    let void = raw.insert_region(Region {
        body,
        solid: false,
        shells: Vec::new(),
    })?;
    let solid = raw.insert_region(Region {
        body,
        solid: true,
        shells: Vec::new(),
    })?;
    raw.body_mut(body)?.regions = vec![void, solid];
    let empty_shell = |region| Shell {
        region,
        faces: Vec::new(),
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    };
    let outer = raw.insert_shell(empty_shell(void))?;
    let inner = raw.insert_shell(empty_shell(solid))?;
    raw.region_mut(void)?.shells.push(outer);
    raw.region_mut(solid)?.shells.push(inner);

    let mut vertices = Vec::with_capacity(8);
    for i in 0..8 {
        vertices.push(raw.insert_vertex(Vertex {
            point: corner(i),
            tolerance: None,
        })?);
    }

    let mut edges: HashMap<(usize, usize), EdgeId> = HashMap::new();
    for (normal, corners) in FACES {
        let plane_frame = Frame::new(
            corner(corners[0]),
            frame.vector_to_world(normal),
            corner(corners[1]) - corner(corners[0]),
            &precision,
        )?;
        let surface = raw
            .geom_mut()
            .add_surface(Surface::Plane(Plane::new(plane_frame)))?;
        let face = raw.insert_face(Face {
            front: Some(outer),
            back: Some(inner),
            loops: Vec::new(),
            surface: Some(surface),
            sense: Sense::Positive,
        })?;
        let lp = raw.insert_loop(Loop {
            face,
            fins: Vec::new(),
        })?;
        for k in 0..4 {
            let (a, b) = (corners[k], corners[(k + 1) % 4]);
            let key = (a.min(b), a.max(b));
            let edge = match edges.get(&key) {
                Some(&edge) => edge,
                None => {
                    let line = Line::through(corner(key.0), corner(key.1), &precision)?;
                    let curve = raw.geom_mut().add_curve(Curve::Line(line))?;
                    let edge = raw.insert_edge(Edge {
                        fins: Vec::new(),
                        start: Some(vertices[key.0]),
                        end: Some(vertices[key.1]),
                        curve: Some(curve),
                        tolerance: None,
                    })?;
                    edges.insert(key, edge);
                    edge
                }
            };
            let sense = if a < b {
                Sense::Positive
            } else {
                Sense::Negative
            };
            let fin = raw.insert_fin(Fin {
                owner: lp,
                edge,
                sense,
            })?;
            raw.edge_mut(edge)?.fins.push(fin);
            raw.loop_mut(lp)?.fins.push(fin);
        }
        raw.face_mut(face)?.loops.push(lp);
        raw.shell_mut(outer)?.faces.push(face);
        raw.shell_mut(inner)?.faces.push(face);
    }
    Ok(body)
}

#[cfg(test)]
mod tests {
    use gensolid_core::Precision;

    use super::*;

    #[test]
    fn block_has_euler_characteristic_two() {
        let mut model = Model::new(Precision::DEFAULT);
        let body = make_block(&mut model, 1.0, 2.0, 3.0, &Frame::WORLD).unwrap();
        let f = model.body_faces(body).unwrap().len();
        let e = model.body_edges(body).unwrap().len();
        let v = model.body_vertices(body).unwrap().len();
        assert_eq!((f, e, v), (6, 12, 8));
        assert_eq!(v + f - e, 2);
    }

    #[test]
    fn rejects_degenerate_size() {
        let mut model = Model::new(Precision::DEFAULT);
        assert!(make_block(&mut model, 0.0, 1.0, 1.0, &Frame::WORLD).is_err());
        assert!(make_block(&mut model, 1.0, f64::NAN, 1.0, &Frame::WORLD).is_err());
    }
}
