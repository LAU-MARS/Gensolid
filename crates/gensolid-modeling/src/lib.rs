//! Checked B-rep modelling and mass properties. This crate never depends on meshes.

mod mass;

pub use gensolid_topo::EdgeSplit;
pub use mass::{MassProperties, body_ask_mass_properties};

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{Error, ErrorCode, Result};
use gensolid_math::{Frame, Point2, Point3, Vec2, Vec3};
use gensolid_topo::{BodyId, EdgeId, Model};

/// Require every enabled immutable body check to pass, with no skipped geometry.
pub fn body_check_complete(model: &Model, body: BodyId) -> Result<()> {
    let report = check_body(model, body, &CheckOptions::default())?;
    if !report.faults.is_empty() {
        return Err(Error::new(
            ErrorCode::InvalidArgument,
            format!(
                "body has {} checker faults: {:?}",
                report.faults.len(),
                report.faults.first()
            ),
        ));
    }
    if report.unchecked != 0 {
        return Err(Error::unimplemented(format!(
            "body has {} unchecked geometric conditions",
            report.unchecked
        )));
    }
    Ok(())
}

/// Build a solid block and require the immutable body checker to pass.
pub fn body_create_solid_block(
    model: &mut Model,
    x: f64,
    y: f64,
    z: f64,
    frame: &Frame,
) -> Result<BodyId> {
    let body = gensolid_topo::make_block(model, x, y, z, frame)?;
    body_check_complete(model, body)?;
    Ok(body)
}

/// Extrude a convex polygon from local Z = 0 to a positive height, and require
/// every enabled immutable body check to pass before returning the body.
pub fn body_create_solid_prism(
    model: &mut Model,
    profile: &[Point2],
    height: f64,
    frame: &Frame,
) -> Result<BodyId> {
    let body = gensolid_topo::make_prism(model, profile, height, frame)?;
    body_check_complete(model, body)?;
    Ok(body)
}

/// Sweep a convex local XY polygon by a world-coordinate displacement. Positive
/// and negative normal components are supported; the completed body is checked.
pub fn body_create_solid_linear_sweep(
    model: &mut Model,
    profile: &[Point2],
    displacement: Vec3,
    frame: &Frame,
) -> Result<BodyId> {
    let body = gensolid_topo::make_linear_sweep(model, profile, displacement, frame)?;
    body_check_complete(model, body)?;
    Ok(body)
}

/// Loft a convex local XY polygon to a uniformly scaled, translated copy in a
/// parallel plane, retaining planar side faces and checking the completed body.
pub fn body_create_solid_homothetic_loft(
    model: &mut Model,
    profile: &[Point2],
    height: f64,
    top_scale: f64,
    top_translation: Vec2,
    frame: &Frame,
) -> Result<BodyId> {
    let body = gensolid_topo::make_homothetic_loft(
        model,
        profile,
        height,
        top_scale,
        top_translation,
        frame,
    )?;
    body_check_complete(model, body)?;
    Ok(body)
}

/// Construct and check a convex polyhedron from world-coordinate vertices and
/// outward, counter-clockwise face index lists.
pub fn body_create_solid_convex_polyhedron(
    model: &mut Model,
    points: &[Point3],
    faces: &[Vec<usize>],
) -> Result<BodyId> {
    let body = gensolid_topo::make_convex_polyhedron(model, points, faces)?;
    body_check_complete(model, body)?;
    Ok(body)
}

/// Split a bounded line edge while preserving all co-edge and wire-shell uses.
/// Every body containing the shared edge must pass all enabled checks both before
/// and after mutation. An edge outside any body cannot bypass body verification.
pub fn edge_split_line(model: &mut Model, edge: EdgeId, fraction: f64) -> Result<EdgeSplit> {
    model.edge(edge)?;
    let mut bodies = Vec::new();
    for body in model.bodies() {
        if model.body_edges(body)?.contains(&edge) {
            bodies.push(body);
        }
    }
    if bodies.is_empty() {
        return Err(Error::invalid_argument("edge is not owned by any body"));
    }
    for &body in &bodies {
        body_check_complete(model, body)?;
    }
    let split = gensolid_topo::split_line_edge(model, edge, fraction)?;
    for body in bodies {
        body_check_complete(model, body)?;
    }
    Ok(split)
}

/// Revolve a convex meridian profile with exact rational circular carriers and check it.
pub fn body_create_solid_revolution(
    model: &mut Model,
    profile: &[Point2],
    angle: f64,
    frame: &Frame,
) -> Result<gensolid_topo::ParametricBody> {
    let result = gensolid_topo::make_revolution(model, profile, angle, frame)?;
    body_check_complete(model, result.body)?;
    Ok(result)
}
/// Loft multiple independently shaped convex sections and check the resulting B-rep.
pub fn body_create_solid_loft_sections(
    model: &mut Model,
    sections: &[gensolid_topo::LoftSection],
    frame: &Frame,
) -> Result<gensolid_topo::ParametricBody> {
    let result = gensolid_topo::make_loft_sections(model, sections, frame)?;
    body_check_complete(model, result.body)?;
    Ok(result)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use gensolid_core::Precision;
    use gensolid_topo::{Body, BodyKind, Region, Shell};
    use proptest::prelude::*;

    use super::*;

    fn rectangle(x: f64, y: f64, reversed: bool) -> Vec<Point2> {
        let mut points = vec![
            Point2::new(-x * 0.5, -y * 0.5),
            Point2::new(x * 0.5, -y * 0.5),
            Point2::new(x * 0.5, y * 0.5),
            Point2::new(-x * 0.5, y * 0.5),
        ];
        if reversed {
            points.reverse();
        }
        points
    }

    fn assert_mass(model: &Model, body: BodyId, volume: f64) {
        body_check_complete(model, body).unwrap();
        let mass = body_ask_mass_properties(model, body).unwrap();
        assert!((mass.volume - volume).abs() <= model.precision().linear() * mass.surface_area);
    }

    fn shared_wire_body(model: &mut Model, edge: EdgeId) -> BodyId {
        let mut raw = model.raw();
        let body = raw
            .insert_body(Body {
                kind: BodyKind::Wire,
                regions: Vec::new(),
            })
            .unwrap();
        let region = raw
            .insert_region(Region {
                body,
                solid: false,
                shells: Vec::new(),
            })
            .unwrap();
        let shell = raw
            .insert_shell(Shell {
                region,
                faces: Vec::new(),
                wire_edges: vec![edge],
                acorn_vertices: Vec::new(),
            })
            .unwrap();
        raw.body_mut(body).unwrap().regions.push(region);
        raw.region_mut(region).unwrap().shells.push(shell);
        body
    }

    #[test]
    fn shared_edge_checks_every_owner_and_preserves_every_wire_use() {
        let mut model = Model::new(Precision::DEFAULT);
        let solid = body_create_solid_block(&mut model, 1.0, 1.0, 1.0, &Frame::WORLD).unwrap();
        let edge = model.body_edges(solid).unwrap()[0];
        let wire = shared_wire_body(&mut model, edge);
        body_check_complete(&model, wire).unwrap();
        let split = edge_split_line(&mut model, edge, 0.5).unwrap();
        for body in [solid, wire] {
            body_check_complete(&model, body).unwrap();
            assert!(model.body_edges(body).unwrap().contains(&split.new_edge));
        }
        let shell = model
            .region(model.body(wire).unwrap().regions[0])
            .unwrap()
            .shells[0];
        let region = model.shell(shell).unwrap().region;
        model.raw().region_mut(region).unwrap().body = solid;
        let end_before = model.edge(split.new_edge).unwrap().end;
        assert_eq!(
            edge_split_line(&mut model, split.new_edge, 0.5)
                .unwrap_err()
                .code(),
            ErrorCode::InvalidArgument
        );
        assert_eq!(model.edge(split.new_edge).unwrap().end, end_before);
    }

    #[test]
    fn near_precision_constructors_and_split_keep_complete_checks() {
        let p = Precision::DEFAULT;
        let mut model = Model::new(p);
        let length = 16.0 * p.linear();
        let profile = rectangle(length, length, false);
        let prism = body_create_solid_prism(&mut model, &profile, length, &Frame::WORLD).unwrap();
        let sweep = body_create_solid_linear_sweep(
            &mut model,
            &profile,
            Vec3::new(length, 0.0, -length),
            &Frame::WORLD,
        )
        .unwrap();
        let loft = body_create_solid_homothetic_loft(
            &mut model,
            &profile,
            length,
            0.5,
            Vec2::new(length, 0.0),
            &Frame::WORLD,
        )
        .unwrap();
        for body in [prism, sweep, loft] {
            body_check_complete(&model, body).unwrap();
        }
        let edge = model.body_edges(prism).unwrap()[0];
        edge_split_line(&mut model, edge, 0.125).unwrap();
        body_check_complete(&model, prism).unwrap();
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]
        #[test]
        fn checked_constructors_retain_known_volumes_in_random_frames(
            x in 0.05f64..5.0, y in 0.05f64..5.0, height in 0.05f64..5.0,
            ox in -10.0f64..10.0, oy in -10.0f64..10.0, oz in -10.0f64..10.0,
            ax in -1.0f64..1.0, ay in -1.0f64..1.0, az in 0.1f64..1.0,
            scale in 0.25f64..2.0, reversed in any::<bool>(), negative_sweep in any::<bool>(),
        ) {
            let p = Precision::DEFAULT;
            let frame = Frame::from_axis(Point3::new(ox, oy, oz), Vec3::new(ax, ay, az), &p).unwrap();
            let profile = rectangle(x, y, reversed);
            let mut model = Model::new(p);
            let prism = body_create_solid_prism(&mut model, &profile, height, &frame).unwrap();
            assert_mass(&model, prism, x * y * height);
            let sign = if negative_sweep { -1.0 } else { 1.0 };
            let sweep = body_create_solid_linear_sweep(&mut model, &profile,
                frame.z() * (sign * height) + frame.x() * x + frame.y() * y, &frame).unwrap();
            assert_mass(&model, sweep, x * y * height);
            let loft = body_create_solid_homothetic_loft(&mut model, &profile, height, scale, Vec2::new(x, -y), &frame).unwrap();
            assert_mass(&model, loft, x * y * height * (1.0 + scale + scale * scale) / 3.0);
            let tetrahedron = [frame.origin(), frame.origin() + frame.x() * x,
                frame.origin() + frame.y() * y, frame.origin() + frame.z() * height];
            let faces = [vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 3]];
            let tetrahedron = body_create_solid_convex_polyhedron(&mut model, &tetrahedron, &faces).unwrap();
            assert_mass(&model, tetrahedron, x * y * height / 6.0);
        }

        #[test]
        fn checked_split_keeps_mass_geometry_and_body_incidence(
            x in 0.05f64..5.0, y in 0.05f64..5.0, z in 0.05f64..5.0,
            axis in (0.1f64..1.0, -1.0f64..1.0, -1.0f64..1.0),
            fraction in 0.05f64..0.95, edge_index in 0usize..12,
        ) {
            let p = Precision::DEFAULT;
            let frame = Frame::from_axis(Point3::new(2.0, -3.0, 1.0), Vec3::new(axis.0, axis.1, axis.2), &p).unwrap();
            let mut model = Model::new(p);
            let body = body_create_solid_block(&mut model, x, y, z, &frame).unwrap();
            let before = body_ask_mass_properties(&model, body).unwrap();
            let edge = model.body_edges(body).unwrap()[edge_index];
            let start = model.vertex(model.edge(edge).unwrap().start.unwrap()).unwrap().point;
            let end = model.vertex(model.edge(edge).unwrap().end.unwrap()).unwrap().point;
            let split = edge_split_line(&mut model, edge, fraction).unwrap();
            prop_assert_eq!(model.body_edges(body).unwrap().len(), 13);
            prop_assert_eq!(model.body_vertices(body).unwrap().len(), 9);
            prop_assert!(model.vertex(split.vertex).unwrap().point.distance(start + (end - start) * fraction) <= p.linear());
            let after = body_ask_mass_properties(&model, body).unwrap();
            prop_assert!((after.volume - before.volume).abs() <= p.linear() * after.surface_area);
        }
    }
}
