//! Independent equations and adversarial seams for the third Rust migration round.
//! Sampling below is a regression oracle, never the production certificate.
#![allow(clippy::unwrap_used)]

use gensolid_api::{Session, SessionOptions, SheetCutOptions, SurfaceConstruction};
use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_geom::{
    BCurve, BSurface,
    backend::{GeometryBackend, SurfaceIntersectionCandidate},
    cylindrical_offset::{certify_cylindrical_graph, certify_offset, offset_cylindrical},
    nurbs_curve_ops::arc,
    rust_surface::{RustGeometryBackend, extrude_nurbs},
};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn cylinder(frame: Frame, radius: f64, angle: f64) -> BSurface {
    let p = Precision::DEFAULT;
    let arc = arc(Point3::ORIGIN, radius, [-angle / 2., angle / 2.], &p).unwrap();
    let profile = BCurve::new(
        2,
        arc.knots()
            .iter()
            .map(|&t| -2. + 7. * ((t + angle / 2.) / angle))
            .collect(),
        arc.points()
            .iter()
            .map(|q| frame.to_world(q.to_vec()))
            .collect(),
        arc.weights().map(<[f64]>::to_vec),
    )
    .unwrap();
    extrude_nurbs(&profile, frame.z() * 2., &p).unwrap()
}

fn cutter(frame: Frame, slope: f64) -> BSurface {
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-8., -8.), (-8., 8.), (8., -8.), (8., 8.)]
            .into_iter()
            .map(|(x, y)| frame.to_world(Vec3::new(x, y, 0.9 + slope * y)))
            .collect(),
        None,
    )
    .unwrap()
}

fn changed(source: &BSurface, points: Vec<Point3>) -> BSurface {
    let (pu, pv) = source.degrees();
    let (nu, nv) = source.counts();
    BSurface::new(
        pu,
        pv,
        source.knots().0.to_vec(),
        source.knots().1.to_vec(),
        nu,
        nv,
        points,
        source.weights().map(<[f64]>::to_vec),
    )
    .unwrap()
}

/// Inspect the published shared edge and both actual fin parameter curves.
fn inspect_cut(session: &Session, frame: Frame, radius: f64, slope: f64) {
    let p = Precision::DEFAULT;
    let model = session.model();
    let body = model.bodies().last().unwrap();
    let faces = model.body_faces(body).unwrap();
    assert_eq!(faces.len(), 2);
    let edge = model
        .body_edges(body)
        .unwrap()
        .into_iter()
        .find(|&id| {
            let e = model.edge(id).unwrap();
            e.fins.len() == 2
        })
        .unwrap();
    let data = model.edge(edge).unwrap();
    let curve = model.geom().curve(data.curve.unwrap()).unwrap();
    let r = curve.range();
    assert_eq!([r.lo, r.hi].map(f64::to_bits), [-2., 5.].map(f64::to_bits));
    for i in 0..=31 {
        let t = r.lo + (r.hi - r.lo) * (i as f64 / 31.);
        let point = curve.eval(t).point;
        let q = frame.to_local(point);
        assert!((q.x.hypot(q.y) - radius).abs() <= p.linear());
        assert!((q.z - 0.9 - slope * q.y).abs() <= p.linear());
        for &fin in &data.fins {
            let geometry = model.fin_geometry(fin).unwrap().unwrap();
            let uv_range = geometry.pcurve.range();
            let fraction =
                (t - geometry.edge_range.lo) / (geometry.edge_range.hi - geometry.edge_range.lo);
            let uv = geometry
                .pcurve
                .eval(uv_range.lo + (uv_range.hi - uv_range.lo) * fraction)
                .point;
            // Profile parameter is affine even if the shared edge runs backwards.
            let start = geometry.pcurve.eval(uv_range.lo).point.x;
            let end = geometry.pcurve.eval(uv_range.hi).point.x;
            assert!((uv.x - (start + (end - start) * fraction)).abs() <= p.solver());
            let face = model.loop_(model.fin(fin).unwrap().owner).unwrap().face;
            let surface = model.face(face).unwrap().surface.unwrap();
            assert!(
                model
                    .geom()
                    .eval_surface(surface, uv.x, uv.y)
                    .unwrap()
                    .point
                    .distance(point)
                    <= p.linear()
            );
        }
    }
}

fn combined(frame: Frame, radius: f64, distance: f64, backend: &dyn GeometryBackend) {
    let p = Precision::DEFAULT;
    let native = RustGeometryBackend::new(p);
    let source = cylinder(frame, radius, 1.2);
    let mut session = Session::start(SessionOptions::default());
    let original = session.body_create_sheet_surface(&source).unwrap();
    let old_faces = session.body_ask_faces(original).unwrap();
    let offset = session
        .body_offset_sheet_backend(original, &native, distance)
        .unwrap();
    let target = session.body_ask_faces(offset).unwrap()[0];
    let tool = session
        .body_create_sheet_surface(&cutter(frame, 0.12))
        .unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let result = session
        .body_split_face_by_sheet_with_options(
            offset,
            target,
            tool_face,
            Some(backend),
            SheetCutOptions {
                allow_rational_oblique: true,
            },
        )
        .unwrap();
    assert!(session.body_check(result.body).unwrap().is_empty());
    assert!(session.body_check(original).unwrap().is_empty());
    assert!(session.body_check(offset).unwrap().is_empty());
    assert_eq!(session.body_ask_faces(original).unwrap(), old_faces);
    inspect_cut(&session, frame, radius + distance, 0.12);
}

#[test]
fn offset_then_oblique_cut_publishes_an_analytic_cylinder_section_and_pcurves() {
    combined(
        Frame::WORLD,
        2.,
        0.25,
        &RustGeometryBackend::new(Precision::DEFAULT),
    );
}

struct Misleading;
impl GeometryBackend for Misleading {
    fn precision(&self) -> Precision {
        Precision::DEFAULT
    }
    fn construct_surface(&self, _: &SurfaceConstruction) -> Result<BSurface> {
        Err(Error::unimplemented("intersection-only adversary"))
    }
    fn offset_surface(&self, _: &BSurface, _: f64) -> Result<BSurface> {
        Err(Error::unimplemented("intersection-only adversary"))
    }
    fn intersect_surfaces(
        &self,
        first: &BSurface,
        second: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        let mut candidates =
            RustGeometryBackend::new(self.precision()).intersect_surfaces(first, second)?;
        for candidate in &mut candidates {
            let wrong = BCurve::new(
                1,
                vec![-2., -2., 5., 5.],
                vec![Point3::new(20., 30., 40.), Point3::new(25., 35., 45.)],
                None,
            )?;
            candidate.curve = wrong.clone();
            candidate.first_uv = wrong.clone();
            candidate.second_uv = wrong;
            candidate.reported_deviation = 0.;
        }
        Ok(candidates)
    }
}

#[test]
fn false_provider_geometry_cannot_replace_the_independently_constructed_cut() {
    // The kernel may ignore proposals and construct its own certified section.
    // A provider's claimed zero deviation must never enter the published edge.
    combined(Frame::WORLD, 2., 0.25, &Misleading);
}

#[test]
fn multispan_offset_cuts_isoparametrically_and_oblique_failure_preserves_tags() {
    let p = Precision::DEFAULT;
    let backend = RustGeometryBackend::new(p);
    let source = cylinder(Frame::WORLD, 2., 2.4);
    assert!(source.counts().0 > 3);
    let mut session = Session::start(SessionOptions::default());
    let original = session.body_create_sheet_surface(&source).unwrap();
    let offset = session
        .body_offset_sheet_backend(original, &backend, 0.25)
        .unwrap();
    let target = session.body_ask_faces(offset).unwrap()[0];
    let tool = session
        .body_create_sheet_surface(&cutter(Frame::WORLD, 0.12))
        .unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let count = session.model().bodies().count();
    let mut untouched = session.clone();
    let result = session.body_split_face_by_sheet_with_options(
        offset,
        target,
        tool_face,
        Some(&backend),
        SheetCutOptions {
            allow_rational_oblique: true,
        },
    );
    assert_eq!(result.unwrap_err().code(), ErrorCode::Unimplemented);
    assert_eq!(session.model().bodies().count(), count);
    let expected = untouched
        .body_create_sheet_surface(&cutter(Frame::WORLD, 0.))
        .unwrap();
    let flat = session
        .body_create_sheet_surface(&cutter(Frame::WORLD, 0.))
        .unwrap();
    assert_eq!(flat, expected);
    let flat_face = session.body_ask_faces(flat).unwrap()[0];
    let split = session
        .body_split_face_by_sheet(offset, target, flat_face, Some(&backend))
        .unwrap();
    assert!(session.body_check(split.body).unwrap().is_empty());
    assert!(session.body_check(offset).unwrap().is_empty());
    inspect_cut(&session, Frame::WORLD, 2.25, 0.);
}

/// Both pieces are exact circular arcs with a shared point and tangent, but
/// their radii and centers differ. Per-piece recognition alone is insufficient.
fn changing_radius() -> BSurface {
    let p = Precision::DEFAULT;
    let first = arc(Point3::ORIGIN, 2., [-0.8, 0.], &p).unwrap();
    let second = arc(Point3::new(0.4, 0., 0.), 1.6, [0., 0.8], &p).unwrap();
    let points = first
        .points()
        .iter()
        .chain(second.points().iter().skip(1))
        .copied()
        .collect();
    let weights = first
        .weights()
        .unwrap()
        .iter()
        .chain(second.weights().unwrap().iter().skip(1))
        .copied()
        .collect();
    let curve = BCurve::new(
        2,
        vec![-2., -2., -2., 1., 1., 5., 5., 5.],
        points,
        Some(weights),
    )
    .unwrap();
    extrude_nurbs(&curve, Vec3::Z * 2., &p).unwrap()
}

#[test]
fn tangent_continuous_pieces_of_different_cylinders_do_not_become_one_cylinder() {
    let p = Precision::DEFAULT;
    let source = changing_radius();
    let error = offset_cylindrical(&source, 0.25, &p).unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
}

#[test]
fn span_local_frames_cannot_certify_a_wide_globally_folded_projection() {
    let p = Precision::DEFAULT;
    let source = cylinder(Frame::WORLD, 2., 3.6);
    assert!(source.counts().0 > 3);
    assert!(certify_cylindrical_graph(&source, &p).is_err());
    let mut session = Session::start(SessionOptions::default());
    assert!(session.body_create_sheet_surface(&source).is_err());
    assert_eq!(session.model().bodies().count(), 0);
}

#[test]
fn shared_seam_corruption_and_axial_drift_do_not_pass_offset_certification() {
    let p = Precision::DEFAULT;
    let source = cylinder(Frame::WORLD, 2., 2.4);
    let output = offset_cylindrical(&source, 0.25, &p).unwrap();
    let mut seam_points = output.points().to_vec();
    seam_points[4].z += p.linear() * 128.;
    seam_points[5].z -= p.linear() * 128.;
    assert!(certify_offset(&source, &changed(&output, seam_points), 0.25, &p).is_err());
    let mut drift_points = source.points().to_vec();
    let last = drift_points.len() - 1;
    drift_points[last].x += 0.125;
    assert!(offset_cylindrical(&changed(&source, drift_points), 0.25, &p).is_err());
}

#[test]
fn nonuniform_simple_knot_refinement_preserves_the_cylinder_and_offset() {
    let p = Precision::DEFAULT;
    let source = cylinder(Frame::WORLD, 2., 2.4);
    let profile = source.iso_curve(false, source.ranges().1.lo).unwrap();
    // These simple knots have different neighboring spans; the existing seam
    // remains double. A fixed local-to-global derivative scale is incorrect.
    let refined = profile.insert_knots(&[-1.3, -0.2, 3.25]).unwrap();
    let source = extrude_nurbs(&refined, Vec3::Z * 2., &p).unwrap();
    certify_cylindrical_graph(&source, &p).unwrap();
    let output = offset_cylindrical(&source, 0.25, &p).unwrap();
    assert_eq!(source.knots(), output.knots());
    assert_eq!(source.weights(), output.weights());
    for &u in source
        .knots()
        .0
        .iter()
        .chain([-1.8, -0.7, 0.2, 2., 4.5].iter())
    {
        for v in [0., 0.31, 1.] {
            let before = source.eval(u, v).point;
            let after = output.eval(u, v).point;
            assert!((before.x.hypot(before.y) - 2.).abs() <= p.linear());
            assert!((after.x.hypot(after.y) - 2.25).abs() <= p.linear());
            assert!((after.z - before.z).abs() <= p.linear());
        }
    }
    let mut session = Session::start(SessionOptions::default());
    let input = session.body_create_sheet_surface(&source).unwrap();
    let output = session
        .body_offset_sheet_backend(input, &RustGeometryBackend::new(p), 0.25)
        .unwrap();
    assert!(session.body_check(input).unwrap().is_empty());
    assert!(session.body_check(output).unwrap().is_empty());
}

#[test]
fn near_tangent_rational_cut_cannot_trade_plane_residual_for_spatial_accuracy() {
    let p = Precision::DEFAULT;
    let epsilon = p.linear() * 10_000.;
    let profile = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![
            Point3::ORIGIN,
            Point3::new(epsilon * 0.1, 0., 1.5),
            Point3::new(epsilon * 0.4, 0., 3.),
        ],
        Some(vec![1., 0.7, 1.]),
    )
    .unwrap();
    let source = extrude_nurbs(&profile, Vec3::Y * 2., &p).unwrap();
    // x + epsilon*y = epsilon crosses every ruling strictly inside (0,2).
    // However a solver-sized plane perturbation gives a much larger spatial
    // displacement along the ruling. Merely checking signed distance is unsafe.
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (5., -5.), (5., 5.)]
            .into_iter()
            .map(|(y, z)| Point3::new(epsilon * (1. - y), y, z))
            .collect(),
        None,
    )
    .unwrap();
    let mut session = Session::start(SessionOptions::default());
    let body = session.body_create_sheet_surface(&source).unwrap();
    let face = session.body_ask_faces(body).unwrap()[0];
    let tool = session.body_create_sheet_surface(&plane).unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let count = session.model().bodies().count();
    let result = session.body_split_face_by_sheet_with_options(
        body,
        face,
        tool_face,
        None,
        SheetCutOptions {
            allow_rational_oblique: true,
        },
    );
    assert_eq!(result.unwrap_err().code(), ErrorCode::Unimplemented);
    assert_eq!(session.model().bodies().count(), count);
    assert!(session.body_check(body).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn checked_offset_and_oblique_cut_compose_in_random_frames(
        radius in 1.0f64..3., distance in -0.25f64..0.4,
        tx in -3.0f64..3., tilt in -0.5f64..0.5,
    ) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(tx,1.,-2.),Vec3::new(tilt,0.3,1.),&p).unwrap();
        combined(frame,radius,distance,&RustGeometryBackend::new(p));
    }
}
