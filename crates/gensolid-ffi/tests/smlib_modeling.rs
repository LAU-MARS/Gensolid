//! Real SMLib-to-kernel integration: owned surfaces and certified bounded face cuts.
#![cfg(feature = "smlib")]
use gensolid::smlib::Backend;
use gensolid_api::{Session, SessionOptions, SurfaceConstruction};
use gensolid_core::Precision;
use gensolid_geom::{BCurve, BSurface};
use gensolid_math::{Point3, Vec3};

#[test]
fn backend_surface_enters_session_and_intersection_modifies_topology() {
    let backend = Backend::new(Precision::DEFAULT).unwrap();
    let mut session = Session::start(SessionOptions::default());
    let curve = BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.4),
            Point3::new(3.0, 0.0, 0.0),
        ],
        None,
    )
    .unwrap();
    let request = SurfaceConstruction::Extrude(curve, Vec3::new(0.0, 2.0, 0.0));
    let sheet = session
        .body_create_sheet_backend(&backend, &request)
        .unwrap();
    assert!(session.body_check(sheet).unwrap().is_empty());
    let body = session
        .body_create_solid_backend_prism(&backend, &request, 1.0)
        .unwrap();
    let target = session.body_ask_faces(body).unwrap()[1];
    let cutter = BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        2,
        2,
        vec![
            Point3::new(-1.0, 0.7, -2.0),
            Point3::new(-1.0, 0.7, 3.0),
            Point3::new(4.0, 0.7, -2.0),
            Point3::new(4.0, 0.7, 3.0),
        ],
        None,
    )
    .unwrap();
    let tool = session.body_create_sheet_surface(&cutter).unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let split = session
        .body_split_face_by_sheet(body, target, tool_face, Some(&backend))
        .unwrap();
    assert_eq!(session.body_ask_faces(split.body).unwrap().len(), 7);
    drop(backend);
    assert!(session.body_check(split.body).unwrap().is_empty());
    assert!(session.body_check(body).unwrap().is_empty());
}

#[test]
fn native_planar_offset_enters_kernel_with_certified_distance() {
    let backend = Backend::new(Precision::DEFAULT).unwrap();
    let mut s = Session::start(SessionOptions::default());
    let plane = BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        2,
        2,
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 2.0, 0.0),
        ],
        None,
    )
    .unwrap();
    let input = s.body_create_sheet_surface(&plane).unwrap();
    let offset = s.body_offset_sheet_backend(input, &backend, 0.5).unwrap();
    drop(backend);
    assert!(s.body_check(offset).unwrap().is_empty());
    let b = s.model().bodies().last().unwrap();
    for v in s.model().body_vertices(b).unwrap() {
        assert!((s.model().vertex(v).unwrap().point.z - 0.5).abs() < Precision::DEFAULT.linear());
    }
}

#[test]
fn native_oblique_section_creates_non_isoparametric_trim_loops() {
    let backend = Backend::new(Precision::DEFAULT).unwrap();
    let mut s = Session::start(SessionOptions::default());
    let profile = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(1.5, 0., 0.4),
            Point3::new(3., 0., 0.),
        ],
        None,
    )
    .unwrap();
    let body = s
        .body_create_solid_backend_prism(
            &backend,
            &SurfaceConstruction::Extrude(profile, Vec3::new(0., 2., 0.)),
            1.,
        )
        .unwrap();
    let target = s.body_ask_faces(body).unwrap()[1];
    let cutter = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (8., -5.), (8., 5.)]
            .into_iter()
            .map(|(x, z)| Point3::new(x, 0.7 + 0.2 * x + 0.2 * z, z))
            .collect(),
        None,
    )
    .unwrap();
    let tool = s.body_create_sheet_surface(&cutter).unwrap();
    let tool_face = s.body_ask_faces(tool).unwrap()[0];
    let split = s
        .body_split_face_by_sheet(body, target, tool_face, Some(&backend))
        .unwrap();
    assert_eq!(s.body_ask_faces(split.body).unwrap().len(), 7);
    drop(backend);
    assert!(s.body_check(split.body).unwrap().is_empty());
    assert!(s.body_check(body).unwrap().is_empty());
}
