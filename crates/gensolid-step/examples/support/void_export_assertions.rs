//! Frozen geometric and topological assertions for the planar void export contract.

use std::collections::HashSet;

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{Error, Result};
use gensolid_geom::{Curve, Surface};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{BodyId, BodyKind, Model, Sense, ShellId};

/// An analytic box, specified independently of the imported/exported topology.
#[derive(Clone, Copy)]
pub struct BoxSpec {
    pub frame: Frame,
    pub sizes: Vec3,
}

pub fn require(condition: bool, message: &str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(Error::invalid_argument(message))
    }
}

impl BoxSpec {
    fn corners(self) -> Vec<Point3> {
        let mut points = Vec::new();
        for x in [-0.5, 0.5] {
            for y in [-0.5, 0.5] {
                for z in [-0.5, 0.5] {
                    points.push(self.frame.to_world(Vec3::new(
                        self.sizes.x * x,
                        self.sizes.y * y,
                        self.sizes.z * z,
                    )));
                }
            }
        }
        points
    }

    fn volume(self) -> f64 {
        self.sizes.x * self.sizes.y * self.sizes.z
    }

    fn area(self) -> f64 {
        2.0 * (self.sizes.x * self.sizes.y
            + self.sizes.x * self.sizes.z
            + self.sizes.y * self.sizes.z)
    }
}

/// The independent OCCT fixture has outer dimensions 4x6x8 and a 2x3x4 cavity.
pub fn occt_boxes() -> Vec<BoxSpec> {
    vec![
        BoxSpec {
            frame: Frame::WORLD,
            sizes: Vec3::new(4.0, 6.0, 8.0),
        },
        BoxSpec {
            frame: Frame::WORLD,
            sizes: Vec3::new(2.0, 3.0, 4.0),
        },
    ]
}

fn shell_points(model: &Model, shell: ShellId) -> Result<Vec<Point3>> {
    let shell = model.shell(shell)?;
    require(
        shell.faces.len() == 6,
        "each analytic box must have six faces",
    )?;
    require(
        shell.wire_edges.is_empty() && shell.acorn_vertices.is_empty(),
        "unexpected wire or isolated vertex",
    )?;
    let mut vertices = HashSet::new();
    let mut edges = HashSet::new();
    for id in &shell.faces {
        let face = model.face(*id)?;
        require(face.loops.len() == 1, "box face must have one loop")?;
        let surface = face
            .surface
            .ok_or_else(|| Error::invalid_argument("missing face surface"))?;
        require(
            matches!(model.geom().surface(surface)?, Surface::Plane(_)),
            "box carrier must be planar",
        )?;
        let lp = model.loop_(face.loops[0])?;
        require(
            lp.face == *id && lp.fins.len() == 4,
            "box loop must have four fins and correct owner",
        )?;
        for fin in &lp.fins {
            let fin = model.fin(*fin)?;
            edges.insert(fin.edge);
            let edge = model.edge(fin.edge)?;
            require(edge.fins.len() == 2, "each box edge must have two uses")?;
            let curve = edge
                .curve
                .ok_or_else(|| Error::invalid_argument("missing edge curve"))?;
            require(
                matches!(model.geom().curve(curve)?, Curve::Line(_)),
                "box edge must be straight",
            )?;
            for vertex in [edge.start, edge.end] {
                vertices
                    .insert(vertex.ok_or_else(|| Error::invalid_argument("unbounded box edge"))?);
            }
        }
    }
    require(
        vertices.len() == 8 && edges.len() == 12,
        "box shell must have eight vertices and twelve edges",
    )?;
    vertices
        .into_iter()
        .map(|id| model.vertex(id).map(|v| v.point))
        .collect()
}

fn same_points(actual: &[Point3], expected: &[Point3], tolerance: f64) -> bool {
    if actual.len() != expected.len() {
        return false;
    }
    let mut used = vec![false; expected.len()];
    for point in actual {
        let Some(index) = expected
            .iter()
            .enumerate()
            .position(|(index, target)| !used[index] && (*point - *target).length() <= tolerance)
        else {
            return false;
        };
        used[index] = true;
    }
    true
}

/// Certifies exact shell membership, analytic coordinates, orientations and volume.
/// The first box is material's outer boundary; all subsequent boxes are cavities.
pub fn certify(model: &Model, body: BodyId, boxes: &[BoxSpec]) -> Result<()> {
    require(!boxes.is_empty(), "empty analytic contract")?;
    require(
        model.bodies().count() == 1,
        "contract requires exactly one body",
    )?;
    let value = model.body(body)?;
    require(value.kind == BodyKind::Solid, "body is not solid")?;
    require(
        value.regions.len() == boxes.len() + 1,
        "finite void region lost or duplicated",
    )?;
    let regions = &value.regions;
    require(
        !model.region(regions[0])?.solid && model.region(regions[1])?.solid,
        "infinite void/material layout changed",
    )?;
    for (index, id) in regions.iter().enumerate() {
        let region = model.region(*id)?;
        require(
            region.body == body && region.solid == (index == 1),
            "region ownership or material classification changed",
        )?;
        require(
            region.shells.len() == if index == 1 { boxes.len() } else { 1 },
            "region boundary components lost or duplicated",
        )?;
    }
    let material = model.region(regions[1])?;
    require(
        material.shells.iter().collect::<HashSet<_>>().len() == boxes.len(),
        "material repeats a boundary shell",
    )?;
    let mut all_faces = HashSet::new();
    let mut unmatched = boxes[1..].to_vec();
    for (index, shell_id) in material.shells.iter().enumerate() {
        let shell = model.shell(*shell_id)?;
        require(
            shell.region == regions[1],
            "material shell owner is inconsistent",
        )?;
        let points = shell_points(model, *shell_id)?;
        let expected = if index == 0 {
            boxes[0]
        } else {
            let matched = unmatched
                .iter()
                .position(|spec| same_points(&points, &spec.corners(), model.precision().linear()))
                .ok_or_else(|| Error::invalid_argument("cavity coordinates changed"))?;
            unmatched.remove(matched)
        };
        require(
            same_points(&points, &expected.corners(), model.precision().linear()),
            "analytic box coordinates changed",
        )?;
        let first = *shell
            .faces
            .first()
            .ok_or_else(|| Error::invalid_argument("empty box shell"))?;
        let paired = model
            .face(first)?
            .front
            .ok_or_else(|| Error::invalid_argument("missing paired void shell"))?;
        let paired_shell = model.shell(paired)?;
        let empty_region = paired_shell.region;
        require(
            empty_region != regions[1] && regions.contains(&empty_region),
            "paired shell has no void region",
        )?;
        require(
            (index == 0) == (empty_region == regions[0]),
            "outer/cavity shell classification changed",
        )?;
        require(
            model.region(empty_region)?.shells == vec![paired],
            "void shell backlink inconsistent",
        )?;
        require(
            paired_shell.faces.len() == shell.faces.len()
                && paired_shell.faces.iter().collect::<HashSet<_>>()
                    == shell.faces.iter().collect::<HashSet<_>>(),
            "shell sides do not share the same face set",
        )?;
        require(
            paired_shell.wire_edges.is_empty() && paired_shell.acorn_vertices.is_empty(),
            "unexpected void wire or vertex",
        )?;
        for face_id in &shell.faces {
            require(
                all_faces.insert(*face_id),
                "face duplicated across boundary components",
            )?;
            let face = model.face(*face_id)?;
            require(
                face.back == Some(*shell_id) && face.front == Some(paired),
                "face/shell side orientation changed",
            )?;
            let surface = face
                .surface
                .ok_or_else(|| Error::invalid_argument("missing plane"))?;
            let Surface::Plane(plane) = model.geom().surface(surface)? else {
                return Err(Error::invalid_argument("nonplanar face"));
            };
            let normal = if face.sense == Sense::Positive {
                plane.frame().z()
            } else {
                -plane.frame().z()
            };
            let signed = normal.dot(plane.frame().origin() - expected.frame.origin());
            require(
                if index == 0 {
                    signed > model.precision().linear()
                } else {
                    signed < -model.precision().linear()
                },
                "face normal does not point out of material",
            )?;
        }
    }
    require(unmatched.is_empty(), "analytic cavity omitted")?;
    require(
        model.body_faces(body)?.len() == boxes.len() * 6
            && model.body_edges(body)?.len() == boxes.len() * 12
            && model.body_vertices(body)?.len() == boxes.len() * 8,
        "global topology count changed",
    )?;
    let report = check_body(model, body, &CheckOptions::default())?;
    require(
        report.faults.is_empty() && report.unchecked == 0,
        "immutable checker found faults or skipped geometry",
    )?;
    let expected = boxes[0].volume() - boxes[1..].iter().map(|b| b.volume()).sum::<f64>();
    let tolerance = model.precision().linear() * boxes.iter().map(|b| b.area()).sum::<f64>();
    require(
        (gensolid_topo::body_volume(model, body)? - expected).abs() <= tolerance,
        "volume differs from independent analytic box difference",
    )
}
