//! Independent invariants and immutable-checker tests for topology construction.

use gensolid_check::{CheckOptions, check_body};
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::{Curve, Line};
use gensolid_math::{Frame, Point2, Point3, Vec2, Vec3};
use proptest::prelude::*;

use gensolid_topo::*;

fn profile(count: usize, x: f64, y: f64, angle: f64, jitter: &[f64]) -> Vec<Point2> {
    (0..count)
        .map(|i| {
            let theta = angle + std::f64::consts::TAU * (i as f64 + jitter[i] * 0.2) / count as f64;
            Point2::new(x * theta.cos(), y * theta.sin())
        })
        .collect()
}

fn area(profile: &[Point2]) -> f64 {
    profile
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let q = profile[(i + 1) % profile.len()];
            p.x * q.y - q.x * p.y
        })
        .sum::<f64>()
        .abs()
        * 0.5
}

fn checked(model: &Model, body: BodyId) -> gensolid_core::Result<()> {
    let report = check_body(model, body, &CheckOptions::default())?;
    assert!(report.is_valid(), "{report:?}");
    assert_eq!(report.unchecked, 0);
    Ok(())
}

#[test]
fn tetrahedron_builder_shares_topology_and_has_analytic_volume() {
    let mut model = Model::new(Precision::DEFAULT);
    let points = [
        Point3::ORIGIN,
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
        Point3::new(0.0, 0.0, 4.0),
    ];
    let faces = vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 3]];
    let body = make_convex_polyhedron(&mut model, &points, &faces).unwrap();
    checked(&model, body).unwrap();
    assert_eq!(model.body_faces(body).unwrap().len(), 4);
    assert_eq!(model.body_edges(body).unwrap().len(), 6);
    assert_eq!(model.body_vertices(body).unwrap().len(), 4);
    assert!((body_volume(&model, body).unwrap() - 4.0).abs() <= model.precision().linear() * 26.0);
    for edge in model.body_edges(body).unwrap() {
        let edge = model.edge(edge).unwrap();
        assert_eq!(edge.fins.len(), 2);
        assert_ne!(
            model.fin(edge.fins[0]).unwrap().sense,
            model.fin(edge.fins[1]).unwrap().sense
        );
    }
}

#[test]
fn bad_polyhedron_inputs_do_not_publish_partial_bodies_or_geometry() {
    let points = [
        Point3::ORIGIN,
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
        Point3::new(0.0, 0.0, 4.0),
    ];
    for faces in [
        vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 9]],
        vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 3, 2]],
        vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![1, 2, 1]],
        vec![vec![0, 2, 1], vec![0, 1, 3], vec![0, 3, 2], vec![0, 2, 1]],
    ] {
        let mut model = Model::new(Precision::DEFAULT);
        assert!(make_convex_polyhedron(&mut model, &points, &faces).is_err());
        assert_eq!(model.bodies().count(), 0);
        assert_eq!(model.geom().curve_count(), 0);
        assert_eq!(model.geom().surface_count(), 0);
    }
    let mut model = Model::new(Precision::DEFAULT);
    let concave = [
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(0.5, 0.5),
        Point2::new(0.0, 2.0),
    ];
    assert!(make_prism(&mut model, &concave, 2.0, &Frame::WORLD).is_err());
    assert_eq!(model.bodies().count(), 0);
}

#[test]
fn shared_edge_split_preserves_fin_order_carrier_and_bounds() {
    let mut model = Model::new(Precision::DEFAULT);
    let body = make_block(&mut model, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
    let edge = model.body_edges(body).unwrap()[0];
    let before = model.edge(edge).unwrap().clone();
    let old_fins = before.fins.clone();
    let bounds = edge_parameter_bounds(&model, edge).unwrap();
    let old_volume = body_volume(&model, body).unwrap();
    let curves = model.geom().curve_count();
    let result = split_line_edge(&mut model, edge, 0.375).unwrap();
    checked(&model, body).unwrap();
    assert_eq!(result.original_edge, edge);
    assert_eq!(model.edge(edge).unwrap().start, before.start);
    assert_eq!(model.edge(edge).unwrap().end, Some(result.vertex));
    assert_eq!(model.edge(result.new_edge).unwrap().end, before.end);
    assert_eq!(model.edge(result.new_edge).unwrap().curve, before.curve);
    assert_eq!(model.edge(edge).unwrap().fins, old_fins);
    assert_eq!(model.geom().curve_count(), curves);
    assert_eq!(model.body_edges(body).unwrap().len(), 13);
    assert_eq!(model.body_vertices(body).unwrap().len(), 9);
    let first = edge_parameter_bounds(&model, edge).unwrap();
    let second = edge_parameter_bounds(&model, result.new_edge).unwrap();
    let tolerance = model.precision().solver();
    assert!((first.lo - bounds.lo).abs() <= tolerance);
    assert!((first.hi - second.lo).abs() <= tolerance);
    assert!((second.hi - bounds.hi).abs() <= tolerance);
    assert!((body_volume(&model, body).unwrap() - old_volume).abs() <= model.precision().linear());
    for fin in old_fins {
        let fin_value = model.fin(fin).unwrap();
        let fins = &model.loop_(fin_value.owner).unwrap().fins;
        let position = fins.iter().position(|id| *id == fin).unwrap();
        let neighbor = if fin_value.sense == Sense::Positive {
            fins[(position + 1) % fins.len()]
        } else {
            fins[(position + fins.len() - 1) % fins.len()]
        };
        assert_eq!(model.fin(neighbor).unwrap().edge, result.new_edge);
        assert_eq!(model.fin(neighbor).unwrap().sense, fin_value.sense);
    }
}

#[test]
fn rejected_splits_do_not_change_existing_topology() {
    let mut model = Model::new(Precision::DEFAULT);
    let body = make_block(&mut model, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
    let edge = model.body_edges(body).unwrap()[0];
    let before = model.edge(edge).unwrap().clone();
    let vertices = model.body_vertices(body).unwrap();
    for fraction in [
        0.0,
        1.0,
        -1.0,
        f64::NAN,
        Precision::DEFAULT.linear() * 0.125,
    ] {
        assert!(split_line_edge(&mut model, edge, fraction).is_err());
        assert_eq!(model.edge(edge).unwrap(), &before);
        assert_eq!(model.body_vertices(body).unwrap(), vertices);
        checked(&model, body).unwrap();
    }
    model.raw().edge_mut(edge).unwrap().start = None;
    assert_eq!(
        split_line_edge(&mut model, edge, 0.5).err().unwrap().code(),
        ErrorCode::Unimplemented
    );
}

#[test]
fn wire_split_updates_its_shell_edge_sequence() {
    let precision = Precision::DEFAULT;
    let mut model = Model::new(precision);
    let (body, shell, edge) = {
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
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            })
            .unwrap();
        raw.body_mut(body).unwrap().regions.push(region);
        raw.region_mut(region).unwrap().shells.push(shell);
        let start = raw
            .insert_vertex(Vertex {
                point: Point3::ORIGIN,
                tolerance: None,
            })
            .unwrap();
        let end = raw
            .insert_vertex(Vertex {
                point: Point3::new(3.0, 0.0, 0.0),
                tolerance: None,
            })
            .unwrap();
        let curve = raw
            .geom_mut()
            .add_curve(Curve::Line(
                Line::new(Point3::ORIGIN, Vec3::X, &precision).unwrap(),
            ))
            .unwrap();
        let edge = raw
            .insert_edge(Edge {
                fins: Vec::new(),
                start: Some(start),
                end: Some(end),
                curve: Some(curve),
                tolerance: None,
            })
            .unwrap();
        raw.shell_mut(shell).unwrap().wire_edges.push(edge);
        (body, shell, edge)
    };
    checked(&model, body).unwrap();
    let split = split_line_edge(&mut model, edge, 0.5).unwrap();
    checked(&model, body).unwrap();
    assert_eq!(
        model.shell(shell).unwrap().wire_edges,
        vec![edge, split.new_edge]
    );
    assert_eq!(model.body_vertices(body).unwrap().len(), 3);
}

#[test]
fn redundant_collinear_profile_corners_preserve_the_extruded_solid() {
    let polygon = [
        Point2::new(-1.0, -1.0),
        Point2::new(0.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 1.0),
    ];
    let precision = Precision::DEFAULT;
    let frame = Frame::from_axis(
        Point3::new(1.0, 2.0, 3.0),
        Vec3::new(1.0, 2.0, 3.0),
        &precision,
    )
    .unwrap();
    let mut model = Model::new(precision);
    let body = make_prism(&mut model, &polygon, 3.0, &frame).unwrap();
    checked(&model, body).unwrap();
    assert_eq!(model.body_vertices(body).unwrap().len(), 8);
    assert!((body_volume(&model, body).unwrap() - 12.0).abs() <= precision.linear() * 32.0);
}

#[test]
fn failed_sweeps_and_lofts_preserve_an_existing_checked_scene() {
    let precision = Precision::DEFAULT;
    let mut model = Model::new(precision);
    let existing = make_block(&mut model, 2.0, 3.0, 4.0, &Frame::WORLD).unwrap();
    let polygon = [
        Point2::new(-1.0, -1.0),
        Point2::new(1.0, -1.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 1.0),
    ];
    let curves = model.geom().curve_count();
    let surfaces = model.geom().surface_count();
    assert!(make_linear_sweep(&mut model, &polygon, Vec3::X, &Frame::WORLD).is_err());
    assert!(
        make_linear_sweep(
            &mut model,
            &polygon,
            Vec3::new(0.0, 0.0, precision.linear() * 0.5),
            &Frame::WORLD
        )
        .is_err()
    );
    for scale in [0.0, -1.0, f64::NAN, precision.linear() * 0.125] {
        assert!(
            make_homothetic_loft(
                &mut model,
                &polygon,
                2.0,
                scale,
                Vec2::new(0.0, 0.0),
                &Frame::WORLD
            )
            .is_err()
        );
    }
    assert!(make_prism(&mut model, &polygon, f64::INFINITY, &Frame::WORLD).is_err());
    assert!(
        make_prism(
            &mut model,
            &[
                Point2::new(f64::NAN, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(0.0, 1.0)
            ],
            2.0,
            &Frame::WORLD
        )
        .is_err()
    );
    assert_eq!(model.bodies().collect::<Vec<_>>(), vec![existing]);
    assert_eq!(model.geom().curve_count(), curves);
    assert_eq!(model.geom().surface_count(), surfaces);
    checked(&model, existing).unwrap();
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn polygon_extrusion_sweep_and_loft_pass_the_immutable_checker(
        count in 3usize..=10,
        jitter in prop::collection::vec(-1.0f64..1.0, 10),
        factors in (1.0f64..3.0, 1.0f64..3.0, 1.0f64..3.0),
        exponent in -3i32..=1,
        translation in (-20.0f64..20.0, -20.0f64..20.0, -20.0f64..20.0),
        axis in (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0),
        theta in -3.0f64..3.0,
        loft_scale in 0.75f64..2.0,
        shear in (-0.5f64..0.5, -0.5f64..0.5),
        clockwise in any::<bool>(),
        negative in any::<bool>(),
        near_degenerate in any::<bool>(),
        fraction in 0.4f64..0.6,
    ) {
        let precision = Precision::DEFAULT;
        let axis = Vec3::new(axis.0, axis.1, axis.2);
        prop_assume!(axis.length() > precision.linear());
        let frame = Frame::from_axis(Point3::new(translation.0, translation.1, translation.2), axis, &precision).unwrap();
        let scale = if near_degenerate { precision.linear() * 8.0 } else { 10.0f64.powi(exponent) };
        let height = scale * factors.2;
        let mut polygon = profile(count, scale * factors.0, scale * factors.1, theta, &jitter);
        if clockwise { polygon.reverse(); }
        let profile_area = area(&polygon);
        let displacement = frame.vector_to_world(Vec3::new(scale * shear.0, scale * shear.1, if negative { -height } else { height }));
        let mut model = Model::new(precision);
        let prism = make_prism(&mut model, &polygon, height, &frame).unwrap();
        let sweep = make_linear_sweep(&mut model, &polygon, displacement, &frame).unwrap();
        let loft = make_homothetic_loft(&mut model, &polygon, height, loft_scale, Vec2::new(scale * shear.0, scale * shear.1), &frame).unwrap();
        for (body, expected) in [(prism, profile_area * height), (sweep, profile_area * height), (loft, profile_area * height * (1.0 + loft_scale + loft_scale * loft_scale) / 3.0)] {
            checked(&model, body).unwrap();
            let faces = model.body_faces(body).unwrap().len();
            let edges = model.body_edges(body).unwrap().len();
            let vertices = model.body_vertices(body).unwrap().len();
            prop_assert_eq!((faces, edges, vertices), (count + 2, count * 3, count * 2));
            prop_assert_eq!(vertices + faces - edges, 2);
            let extent = scale * 10.0;
            prop_assert!((body_volume(&model, body).unwrap() - expected).abs() <= precision.linear() * extent * extent * faces as f64);
        }
        let edge = model.body_edges(prism).unwrap()[0];
        let old_volume = body_volume(&model, prism).unwrap();
        split_line_edge(&mut model, edge, fraction).unwrap();
        checked(&model, prism).unwrap();
        prop_assert!((body_volume(&model, prism).unwrap() - old_volume).abs() <= precision.linear() * scale * scale * count as f64);
    }
}
