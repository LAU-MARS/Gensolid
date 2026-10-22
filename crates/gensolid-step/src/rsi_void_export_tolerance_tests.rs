//! Export precision regressions: local input tolerances are not STEP metadata.

#![cfg(test)]

use std::collections::{HashMap, HashSet};

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision, Tolerance};
use gensolid_geom::{Curve, Line, Plane, Surface};
use gensolid_math::{Frame, Point3, Vec3};
use gensolid_topo::{BodyId, Model, VertexId};
use proptest::prelude::*;

use crate::parse::{Instance, Value, bool_of, name_skipped, nums, parse_document, ref_of, refs};
use crate::{StepStatus, read_solids, validity_report, write_planar_solid};

fn fixture() -> (Model, BodyId) {
    let model = read_solids(
        include_str!("../tests/fixtures/hollow_block.step"),
        Precision::DEFAULT,
    )
    .unwrap();
    let body = model.bodies().next().unwrap();
    (model, body)
}

fn outer_vertices(model: &Model, body: BodyId) -> HashSet<VertexId> {
    let empty = model.body(body).unwrap().regions[0];
    let shell = model.region(empty).unwrap().shells[0];
    let mut vertices = HashSet::new();
    for face in &model.shell(shell).unwrap().faces {
        for loop_id in &model.face(*face).unwrap().loops {
            for fin in &model.loop_(*loop_id).unwrap().fins {
                let (a, b) = model.fin_vertices(*fin).unwrap();
                vertices.extend([a.unwrap(), b.unwrap()]);
            }
        }
    }
    assert_eq!(vertices.len(), 8);
    vertices
}

#[test]
fn local_tolerance_counterexample_is_explicitly_refused() {
    let (mut model, body) = fixture();
    let p = *model.precision();
    let d = 16.0 * p.linear();
    for vertex in outer_vertices(&model, body) {
        let mut raw = model.raw();
        let vertex = raw.vertex_mut(vertex).unwrap();
        vertex.point.x -= vertex.point.x.signum() * d;
        vertex.point.y -= vertex.point.y.signum() * d;
        vertex.point.z -= vertex.point.z.signum() * d;
        vertex.tolerance = Some(Tolerance::new(4.0 * d, &p).unwrap());
    }
    let report = check_body(&model, body, &CheckOptions::default()).unwrap();
    println!("source check_body: {report:?}");
    assert!(report.is_valid());
    assert_eq!(report.unchecked, 0);
    let relationship = crate::voids::validate(&model, body);
    println!("source voids::validate: {relationship:?}");
    assert!(relationship.is_ok());
    let output = write_planar_solid(&model, body);
    match &output {
        Ok(text) => {
            println!("writer: Ok ({} bytes)", text.len());
            println!("output validity: {:?}", validity_report(text, p));
            println!("output reimport: {:?}", read_solids(text, p).err());
        }
        Err(error) => println!("writer: {error:?}; no output to reimport"),
    }
    let error = output.expect_err("incompatible retained planes must not be exported");
    assert_eq!(error.code(), ErrorCode::Unimplemented);
}

fn checked(model: &Model, body: BodyId) {
    let report = check_body(model, body, &CheckOptions::default()).unwrap();
    assert!(report.is_valid(), "{report:?}");
    assert_eq!(report.unchecked, 0);
    crate::voids::validate(model, body).unwrap();
}

fn metadata(model: &mut Model, body: BodyId, radius: f64) {
    let tolerance = Some(Tolerance::new(radius, model.precision()).unwrap());
    for vertex in model.body_vertices(body).unwrap() {
        model.raw().vertex_mut(vertex).unwrap().tolerance = tolerance;
    }
    for edge in model.body_edges(body).unwrap() {
        model.raw().edge_mut(edge).unwrap().tolerance = tolerance;
    }
}

fn transformed_fixture(scales: Vec3, frame: Frame) -> (Model, BodyId) {
    let (mut model, body) = fixture();
    let p = *model.precision();
    let transform = |point: Point3| {
        frame.to_world(Vec3::new(
            point.x * scales.x,
            point.y * scales.y,
            point.z * scales.z,
        ))
    };
    for vertex in model.body_vertices(body).unwrap() {
        let point = transform(model.vertex(vertex).unwrap().point);
        model.raw().vertex_mut(vertex).unwrap().point = point;
    }
    for edge in model.body_edges(body).unwrap() {
        let value = model.edge(edge).unwrap();
        let line = Line::through(
            model.vertex(value.start.unwrap()).unwrap().point,
            model.vertex(value.end.unwrap()).unwrap().point,
            &p,
        )
        .unwrap();
        let curve = model.raw().geom_mut().add_curve(Curve::Line(line)).unwrap();
        model.raw().edge_mut(edge).unwrap().curve = Some(curve);
    }
    for face in model.body_faces(body).unwrap() {
        let Surface::Plane(plane) = model
            .geom()
            .surface(model.face(face).unwrap().surface.unwrap())
            .unwrap()
        else {
            panic!("fixture must be planar");
        };
        let placement = Frame::new(
            transform(plane.frame().origin()),
            frame.vector_to_world(plane.frame().z()),
            frame.vector_to_world(plane.frame().x()),
            &p,
        )
        .unwrap();
        let surface = model
            .raw()
            .geom_mut()
            .add_surface(Surface::Plane(Plane::new(placement)))
            .unwrap();
        model.raw().face_mut(face).unwrap().surface = Some(surface);
    }
    checked(&model, body);
    (model, body)
}

fn shrink_outer(model: &mut Model, body: BodyId, frame: Frame, d: f64) {
    for vertex in outer_vertices(model, body) {
        let local = frame.to_local(model.vertex(vertex).unwrap().point);
        let point = frame.to_world(Vec3::new(
            local.x - local.x.signum() * d,
            local.y - local.y.signum() * d,
            local.z - local.z.signum() * d,
        ));
        model.raw().vertex_mut(vertex).unwrap().point = point;
    }
}

type Entities = HashMap<u32, Instance>;

fn args<'a>(data: &'a Entities, id: u32, kind: &str) -> &'a [Value] {
    assert_eq!(data[&id].type_name, kind);
    name_skipped(&data[&id].args).unwrap()
}

fn triple(data: &Entities, id: u32, kind: &str) -> Vec3 {
    let values = nums(&args(data, id, kind)[0]).unwrap();
    assert_eq!(values.len(), 3);
    Vec3::new(values[0], values[1], values[2])
}

fn point(data: &Entities, id: u32) -> Point3 {
    let value = triple(data, id, "CARTESIAN_POINT");
    Point3::new(value.x, value.y, value.z)
}

fn vertex(data: &Entities, id: u32) -> Point3 {
    point(data, ref_of(&args(data, id, "VERTEX_POINT")[0]).unwrap())
}

// Inspect raw exported coordinates and carriers. This does not use the STEP
// reader, the writer's precision guard, or source-model geometry queries.
fn raw_shell(data: &Entities, shell: u32, center: Point3, p: Precision) -> f64 {
    let faces = refs(&args(data, shell, "CLOSED_SHELL")[0]).unwrap();
    assert_eq!(faces.len(), 6);
    let mut volume = 0.0;
    for face in faces {
        let face = args(data, face, "ADVANCED_FACE");
        let plane = args(data, ref_of(&face[1]).unwrap(), "PLANE");
        let placement = args(data, ref_of(&plane[0]).unwrap(), "AXIS2_PLACEMENT_3D");
        let origin = point(data, ref_of(&placement[0]).unwrap());
        let axis = triple(data, ref_of(&placement[1]).unwrap(), "DIRECTION");
        assert!((axis.length() - 1.0).abs() <= p.angular());
        let normal = axis
            * if bool_of(&face[2]).unwrap() {
                1.0
            } else {
                -1.0
            };
        let bounds = refs(&face[0]).unwrap();
        assert_eq!(bounds.len(), 1);
        let bound = args(data, bounds[0], "FACE_OUTER_BOUND");
        assert!(bool_of(&bound[1]).unwrap());
        let fins = refs(&args(data, ref_of(&bound[0]).unwrap(), "EDGE_LOOP")[0]).unwrap();
        assert_eq!(fins.len(), 4);
        let mut cycle = Vec::new();
        for fin in fins {
            let fin = args(data, fin, "ORIENTED_EDGE");
            let edge = args(data, ref_of(&fin[2]).unwrap(), "EDGE_CURVE");
            assert!(bool_of(&edge[3]).unwrap());
            let a = ref_of(&edge[0]).unwrap();
            let b = ref_of(&edge[1]).unwrap();
            let line = args(data, ref_of(&edge[2]).unwrap(), "LINE");
            let line_origin = point(data, ref_of(&line[0]).unwrap());
            let vector = args(data, ref_of(&line[1]).unwrap(), "VECTOR");
            let direction = triple(data, ref_of(&vector[0]).unwrap(), "DIRECTION");
            assert!((direction.length() - 1.0).abs() <= p.angular());
            for v in [a, b] {
                let q = vertex(data, v);
                // A straight segment's signed plane distance is affine, so its
                // maximum absolute distance is attained at an endpoint.
                assert!((q - origin).dot(axis).abs() <= p.default_local());
                assert!((q - line_origin).cross(direction).length() <= p.linear());
            }
            cycle.push(if bool_of(&fin[3]).unwrap() {
                (a, b)
            } else {
                (b, a)
            });
        }
        let mut area = Vec3::ZERO;
        for i in 0..cycle.len() {
            let next = (i + 1) % cycle.len();
            assert_eq!(cycle[i].1, cycle[next].0);
            area += (vertex(data, cycle[i].0) - center).cross(vertex(data, cycle[next].0) - center)
                * 0.5;
        }
        let q = vertex(data, cycle[0].0);
        assert!(area.dot(normal) > p.linear() * p.linear());
        assert!(
            (q - center).dot(normal) > p.linear(),
            "base shell must face outward"
        );
        volume += (q - center).dot(area) / 3.0;
    }
    volume
}

fn verify_success(model: &Model, body: BodyId, frame: Frame, scales: Vec3, d: f64) {
    checked(model, body);
    let p = *model.precision();
    let text = write_planar_solid(model, body).unwrap();
    let data = parse_document(&text).unwrap();
    let exported: Vec<_> = data
        .iter()
        .filter(|(_, e)| e.type_name == "VERTEX_POINT")
        .map(|(id, _)| vertex(&data, *id))
        .collect();
    assert_eq!(exported.len(), 16);
    for v in model.body_vertices(body).unwrap() {
        let source = model.vertex(v).unwrap().point;
        assert_eq!(
            exported
                .iter()
                .filter(|q| q.distance(source) <= p.solver())
                .count(),
            1
        );
    }
    let brep = data
        .values()
        .find(|e| e.type_name == "BREP_WITH_VOIDS")
        .unwrap();
    let brep = name_skipped(&brep.args).unwrap();
    let outer = ref_of(&brep[0]).unwrap();
    let uses = refs(&brep[1]).unwrap();
    assert_eq!(uses.len(), 1);
    let cavity = args(&data, uses[0], "ORIENTED_CLOSED_SHELL");
    assert!(!bool_of(&cavity[2]).unwrap());
    let cavity_shell = ref_of(&cavity[1]).unwrap();
    assert_ne!(outer, cavity_shell);
    let actual = raw_shell(&data, outer, frame.origin(), p)
        - raw_shell(&data, cavity_shell, frame.origin(), p);
    let dimensions = Vec3::new(
        4.0 * scales.x - 2.0 * d,
        6.0 * scales.y - 2.0 * d,
        8.0 * scales.z - 2.0 * d,
    );
    let expected =
        dimensions.x * dimensions.y * dimensions.z - 24.0 * scales.x * scales.y * scales.z;
    let budget = p.linear()
        * (dimensions.x * dimensions.y + dimensions.y * dimensions.z + dimensions.z * dimensions.x)
        * 8.0;
    assert!(
        (actual - expected).abs() <= budget,
        "raw boundary integral {actual}, analytic {expected}"
    );
    let report = validity_report(&text, p);
    assert_eq!(report.status, StepStatus::Valid, "{report:?}");
    let imported = read_solids(&text, p).unwrap();
    let imported_body = imported.bodies().next().unwrap();
    checked(&imported, imported_body);
    assert_eq!(imported.body(imported_body).unwrap().regions.len(), 3);
}

#[test]
fn exact_geometry_with_large_tolerance_metadata_remains_exportable() {
    let (mut model, body) = fixture();
    let radius = model.precision().size() * 0.25;
    metadata(&mut model, body, radius);
    verify_success(&model, body, Frame::WORLD, Vec3::new(1.0, 1.0, 1.0), 0.0);
    assert!(
        model.body_vertices(body).unwrap().iter().all(|v| model
            .vertex(*v)
            .unwrap()
            .tolerance
            .is_some())
    );
}

#[test]
fn regenerated_lines_can_remove_tolerated_input_curve_displacement() {
    let (mut model, body) = fixture();
    let p = *model.precision();
    let edge_id = model.body_edges(body).unwrap()[0];
    let Curve::Line(line) = *model
        .geom()
        .curve(model.edge(edge_id).unwrap().curve.unwrap())
        .unwrap()
    else {
        panic!("fixture must have line edges");
    };
    let d = 16.0 * p.linear();
    let displaced = Line::new(line.origin() + Vec3::new(d, d, d), line.direction(), &p).unwrap();
    let curve = model
        .raw()
        .geom_mut()
        .add_curve(Curve::Line(displaced))
        .unwrap();
    model.raw().edge_mut(edge_id).unwrap().curve = Some(curve);
    assert!(
        !check_body(&model, body, &CheckOptions::default())
            .unwrap()
            .is_valid()
    );
    model.raw().edge_mut(edge_id).unwrap().tolerance = Some(Tolerance::new(4.0 * d, &p).unwrap());
    verify_success(&model, body, Frame::WORLD, Vec3::new(1.0, 1.0, 1.0), 0.0);
}

#[test]
fn retained_planes_set_the_output_precision_boundary() {
    let p = Precision::DEFAULT;
    for factor in [0.0, 0.25, 0.9, 1.1, 2.0, 32.0] {
        let (mut model, body) = fixture();
        let d = factor * p.default_local();
        metadata(&mut model, body, p.linear() * 128.0);
        shrink_outer(&mut model, body, Frame::WORLD, d);
        checked(&model, body);
        if factor <= 0.9 {
            verify_success(&model, body, Frame::WORLD, Vec3::new(1.0, 1.0, 1.0), d);
        } else {
            assert_eq!(
                write_planar_solid(&model, body).unwrap_err().code(),
                ErrorCode::Unimplemented
            );
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn local_tolerances_do_not_hide_export_errors_in_random_frames(
        factors in (0.5f64..2.0, 0.5f64..2.0, 0.5f64..2.0),
        translation in (-20.0f64..20.0, -20.0f64..20.0, -20.0f64..20.0),
        axis in (-2.0f64..2.0, -2.0f64..2.0),
        roll in -3.0f64..3.0,
        exponent in -5i32..=1,
        tiny in any::<bool>(),
        radius in 128.0f64..4096.0,
        boundary in prop::sample::select(vec![0.0, 0.25, 0.9, 1.0, 1.1, 2.0, 32.0]),
    ) {
        let p = Precision::DEFAULT;
        let scale = if tiny { p.linear() * 256.0 } else { 10.0f64.powi(exponent) };
        let scales = Vec3::new(factors.0 * scale, factors.1 * scale, factors.2 * scale);
        let basis = Frame::from_axis(Point3::new(translation.0, translation.1, translation.2), Vec3::new(axis.0, axis.1, 1.0), &p).unwrap();
        let frame = Frame::new(basis.origin(), basis.z(), basis.x() * roll.cos() + basis.y() * roll.sin(), &p).unwrap();
        let (mut model, body) = transformed_fixture(scales, frame);
        let d = boundary * p.default_local();
        metadata(&mut model, body, radius * p.default_local());
        shrink_outer(&mut model, body, frame, d);
        checked(&model, body);
        if boundary <= 0.9 {
            verify_success(&model, body, frame, scales, d);
        } else if boundary >= 1.1 {
            prop_assert_eq!(write_planar_solid(&model, body).unwrap_err().code(), ErrorCode::Unimplemented);
        } else {
            // At exactly the default tolerance, floating-point placement can
            // land on either side. Every successful file must still independently
            // satisfy raw carrier checks and complete reimport validation.
            match write_planar_solid(&model, body) {
                Ok(_) => verify_success(&model, body, frame, scales, d),
                Err(error) => prop_assert_eq!(error.code(), ErrorCode::Unimplemented),
            }
        }
    }
}
