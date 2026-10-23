//! Pure Rust geometry must enter the same checked Session pipeline as SMLib.
//! These are independent integration tests; the original SMLib tests stay intact.
#![allow(clippy::unwrap_used)]

use gensolid_api::{Session, SessionOptions, SurfaceConstruction, Tag};
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::{BCurve, BSurface, rust_surface::RustGeometryBackend};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn profile(frame: Frame, width: f64, bend: f64) -> BCurve {
    BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        [
            Vec3::ZERO,
            Vec3::new(width / 2., 0., bend),
            Vec3::new(width, 0., 0.),
        ]
        .into_iter()
        .map(|q| frame.to_world(q))
        .collect(),
        None,
    )
    .unwrap()
}

fn plane(points: Vec<Point3>) -> BSurface {
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        points,
        None,
    )
    .unwrap()
}

fn cutter(frame: Frame, width: f64, height: f64, at: f64, tilt: f64) -> BSurface {
    plane(
        [
            (-1., -2.),
            (-1., height + 2.),
            (width + 1., -2.),
            (width + 1., height + 2.),
        ]
        .into_iter()
        .map(|(x, z)| frame.to_world(Vec3::new(x, at + tilt * x + tilt * z, z)))
        .collect(),
    )
}

#[derive(Debug, PartialEq, Eq)]
struct BodySnapshot {
    faces: Vec<Tag>,
    edges: Vec<Tag>,
    vertices: Vec<(Tag, [u64; 3])>,
}

fn snapshot(s: &mut Session, body: Tag) -> BodySnapshot {
    BodySnapshot {
        faces: s.body_ask_faces(body).unwrap(),
        edges: s.body_ask_edges(body).unwrap(),
        vertices: s
            .body_ask_vertices(body)
            .unwrap()
            .into_iter()
            .map(|v| (v, s.vertex_ask_point(v).unwrap().map(f64::to_bits)))
            .collect(),
    }
}

#[test]
fn rust_surface_enters_session_and_intersection_modifies_topology() {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    let mut session = Session::start(SessionOptions::default());
    let request =
        SurfaceConstruction::Extrude(profile(Frame::WORLD, 3., 0.4), Vec3::new(0., 2., 0.));
    let sheet = session
        .body_create_sheet_backend(&backend, &request)
        .unwrap();
    assert!(session.body_check(sheet).unwrap().is_empty());
    let body = session
        .body_create_solid_backend_prism(&backend, &request, 1.)
        .unwrap();
    let before = snapshot(&mut session, body);
    let target = session.body_ask_faces(body).unwrap()[1];
    let tool = session
        .body_create_sheet_surface(&cutter(Frame::WORLD, 3., 1., 0.7, 0.))
        .unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let split = session
        .body_split_face_by_sheet(body, target, tool_face, Some(&backend))
        .unwrap();
    assert_eq!(session.body_ask_faces(split.body).unwrap().len(), 7);
    assert!(session.body_check(split.body).unwrap().is_empty());
    assert!(session.body_check(body).unwrap().is_empty());
    assert_eq!(snapshot(&mut session, body), before);
}

#[test]
fn rust_planar_offset_enters_kernel_with_certified_distance() {
    let mut s = Session::start(SessionOptions::default());
    let surface = plane(vec![
        Point3::new(0., 0., 0.),
        Point3::new(0., 2., 0.),
        Point3::new(3., 0., 0.),
        Point3::new(3., 2., 0.),
    ]);
    let input = s.body_create_sheet_surface(&surface).unwrap();
    let before = snapshot(&mut s, input);
    let offset = {
        let backend = RustGeometryBackend::new(Precision::DEFAULT);
        s.body_offset_sheet_backend(input, &backend, 0.5).unwrap()
    };
    // Kernel geometry remains owned after the provider's lifetime ends.
    assert!(s.body_check(offset).unwrap().is_empty());
    for v in s.body_ask_vertices(offset).unwrap() {
        assert!((s.vertex_ask_point(v).unwrap()[2] - 0.5).abs() < Precision::DEFAULT.linear());
    }
    assert_eq!(snapshot(&mut s, input), before);
}

#[test]
fn rust_oblique_section_creates_non_isoparametric_trim_loops() {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    let mut s = Session::start(SessionOptions::default());
    let body = s
        .body_create_solid_backend_prism(
            &backend,
            &SurfaceConstruction::Extrude(profile(Frame::WORLD, 3., 0.4), Vec3::new(0., 2., 0.)),
            1.,
        )
        .unwrap();
    let before = snapshot(&mut s, body);
    let target = s.body_ask_faces(body).unwrap()[1];
    let tool = s
        .body_create_sheet_surface(&plane(
            [(-5., -5.), (-5., 5.), (8., -5.), (8., 5.)]
                .into_iter()
                .map(|(x, z)| Point3::new(x, 0.7 + 0.2 * x + 0.2 * z, z))
                .collect(),
        ))
        .unwrap();
    let tool_face = s.body_ask_faces(tool).unwrap()[0];
    let split = s
        .body_split_face_by_sheet(body, target, tool_face, Some(&backend))
        .unwrap();
    assert_eq!(s.body_ask_faces(split.body).unwrap().len(), 7);
    assert!(s.body_check(split.body).unwrap().is_empty());
    assert!(s.body_check(body).unwrap().is_empty());
    assert_eq!(snapshot(&mut s, body), before);
}

#[test]
fn unsupported_curved_sweep_and_curved_offset_leave_session_unchanged() {
    let p = Precision::DEFAULT;
    let backend = RustGeometryBackend::new(p);
    let mut s = Session::start(SessionOptions::default());
    let source = profile(Frame::WORLD, 3., 0.4);
    let body = s
        .body_create_sheet_backend(
            &backend,
            &SurfaceConstruction::Extrude(source.clone(), Vec3::new(0., 2., 0.)),
        )
        .unwrap();
    let before = snapshot(&mut s, body);
    let count = s.model().bodies().count();
    let error = s
        .body_create_sheet_backend(
            &backend,
            &SurfaceConstruction::Sweep(source.clone(), source),
        )
        .unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert_eq!(s.model().bodies().count(), count);
    let error = s
        .body_offset_sheet_backend(body, &backend, 0.5)
        .unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert_eq!(s.model().bodies().count(), count);
    assert_eq!(snapshot(&mut s, body), before);
    assert!(s.body_check(body).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn rust_extrusions_and_bounded_cuts_in_random_frames(
        width in 1.0f64..5., depth in 1.0f64..5., height in 0.2f64..2.,
        bend in -0.1f64..0.1, ox in -10.0f64..10., ax in -0.5f64..0.5,
        ay in -0.5f64..0.5, fraction in 0.3f64..0.7, tilt in -0.02f64..0.02,
    ) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(ox, 1., 2.), Vec3::new(ax, ay, 1.), &p).unwrap();
        let backend = RustGeometryBackend::new(p);
        let mut s = Session::start(SessionOptions::default());
        let body = s.body_create_solid_backend_prism(&backend,
            &SurfaceConstruction::Extrude(profile(frame, width, bend), frame.y()*depth), height).unwrap();
        let before = snapshot(&mut s, body);
        let face = s.body_ask_faces(body).unwrap()[1];
        let tool = s.body_create_sheet_surface(&cutter(frame, width, height, depth*fraction, tilt)).unwrap();
        let tool_face = s.body_ask_faces(tool).unwrap()[0];
        let split = s.body_split_face_by_sheet(body, face, tool_face, Some(&backend)).unwrap();
        prop_assert_eq!(s.body_ask_faces(split.body).unwrap().len(), 7);
        prop_assert!(s.body_check(split.body).unwrap().is_empty());
        prop_assert!(s.body_check(body).unwrap().is_empty());
        prop_assert_eq!(snapshot(&mut s, body), before);
    }

    #[test]
    fn rust_planar_offsets_have_expected_signed_distance_in_random_frames(
        width in 0.5f64..5., depth in 0.5f64..5., distance in -2.0f64..2.,
        ox in -10.0f64..10., ax in -0.5f64..0.5, ay in -0.5f64..0.5,
    ) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(ox, 1., 2.), Vec3::new(ax, ay, 1.), &p).unwrap();
        let surface = plane([Vec3::ZERO, Vec3::new(0.,depth,0.), Vec3::new(width,0.,0.), Vec3::new(width,depth,0.)]
            .into_iter().map(|q| frame.to_world(q)).collect());
        let backend = RustGeometryBackend::new(p);
        let mut s = Session::start(SessionOptions::default());
        let body = s.body_create_sheet_surface(&surface).unwrap();
        let before = snapshot(&mut s, body);
        let offset = s.body_offset_sheet_backend(body, &backend, distance).unwrap();
        prop_assert!(s.body_check(offset).unwrap().is_empty());
        for tag in s.body_ask_vertices(offset).unwrap() {
            let [x,y,z] = s.vertex_ask_point(tag).unwrap();
            prop_assert!((frame.to_local(Point3::new(x,y,z)).z - distance).abs() <= p.linear());
        }
        prop_assert_eq!(snapshot(&mut s, body), before);
    }

    #[test]
    fn near_degenerate_extrusions_are_refused_without_partial_bodies(factor in 0.01f64..0.99) {
        let p = Precision::DEFAULT;
        let backend = RustGeometryBackend::new(p);
        let mut s = Session::start(SessionOptions::default());
        let existing = s.body_create_solid_block(2., 3., 4., None).unwrap();
        let before = snapshot(&mut s, existing);
        let count = s.model().bodies().count();
        let request = SurfaceConstruction::Extrude(profile(Frame::WORLD, 3., 0.4), Vec3::Y*(p.linear()*factor));
        prop_assert!(s.body_create_sheet_backend(&backend, &request).is_err());
        prop_assert_eq!(s.model().bodies().count(), count);
        prop_assert_eq!(snapshot(&mut s, existing), before);
        prop_assert!(s.body_check(existing).unwrap().is_empty());
    }
}
