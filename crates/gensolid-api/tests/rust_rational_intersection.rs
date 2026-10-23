//! Rational spline candidates must pass independent modelling checks and roll back failures.
#![allow(clippy::unwrap_used)]
use gensolid_api::{Session, SessionOptions, SurfaceConstruction};
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::{BCurve, BSurface, rust_surface::RustGeometryBackend};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn profile(frame: Frame, width: f64, weight: f64, multispan: bool) -> BCurve {
    let (knots, points, weights) = if multispan {
        (
            vec![0., 0., 0., 0.5, 0.5, 1., 1., 1.],
            vec![(0., 0.), (0.25, 0.02), (0.5, 0.), (0.75, -0.02), (1., 0.)],
            vec![1., weight, 1., weight, 1.],
        )
    } else {
        (
            vec![0., 0., 0., 1., 1., 1.],
            vec![(0., 0.), (0.5, 0.04), (1., 0.)],
            vec![1., weight, 1.],
        )
    };
    BCurve::new(
        2,
        knots,
        points
            .into_iter()
            .map(|(x, z)| frame.to_world(Vec3::new(x * width, 0., z)))
            .collect(),
        Some(weights),
    )
    .unwrap()
}

fn tool(frame: Frame, width: f64, depth: f64, fraction: f64, slope: f64) -> BSurface {
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        [(-5., -5.), (-5., 5.), (width + 5., -5.), (width + 5., 5.)]
            .into_iter()
            .map(|(x, z)| frame.to_world(Vec3::new(x, depth * fraction + slope * x, z)))
            .collect(),
        None,
    )
    .unwrap()
}

fn cut(frame: Frame, width: f64, depth: f64, weight: f64, fraction: f64, multispan: bool) {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    let mut s = Session::start(SessionOptions::default());
    let body = s
        .body_create_sheet_backend(
            &backend,
            &SurfaceConstruction::Extrude(
                profile(frame, width, weight, multispan),
                frame.y() * depth,
            ),
        )
        .unwrap();
    let old_faces = s.body_ask_faces(body).unwrap();
    let old_edges = s.body_ask_edges(body).unwrap();
    let target = old_faces[0];
    let cutter = s
        .body_create_sheet_surface(&tool(frame, width, depth, fraction, 0.))
        .unwrap();
    let tool_face = s.body_ask_faces(cutter).unwrap()[0];
    let split = s
        .body_split_face_by_sheet(body, target, tool_face, Some(&backend))
        .unwrap();
    assert_eq!(s.body_ask_faces(split.body).unwrap().len(), 2);
    assert!(s.body_check(split.body).unwrap().is_empty());
    assert!(s.body_check(body).unwrap().is_empty());
    assert_eq!(s.body_ask_faces(body).unwrap(), old_faces);
    assert_eq!(s.body_ask_edges(body).unwrap(), old_edges);
}

#[test]
fn rational_single_and_multispan_intersections_split_checked_session_sheets() {
    for multispan in [false, true] {
        cut(Frame::WORLD, 3., 2., 0.8, 0.4, multispan);
    }
}

#[test]
fn rational_oblique_geometry_is_not_mistaken_for_a_modelling_certificate() {
    let backend = RustGeometryBackend::new(Precision::DEFAULT);
    let mut s = Session::start(SessionOptions::default());
    let body = s
        .body_create_sheet_backend(
            &backend,
            &SurfaceConstruction::Extrude(profile(Frame::WORLD, 3., 0.8, false), Vec3::Y * 2.),
        )
        .unwrap();
    let faces = s.body_ask_faces(body).unwrap();
    let edges = s.body_ask_edges(body).unwrap();
    let cutter = s
        .body_create_sheet_surface(&tool(Frame::WORLD, 3., 2., 0.4, 0.1))
        .unwrap();
    let tool_face = s.body_ask_faces(cutter).unwrap()[0];
    let count = s.model().bodies().count();
    let error = s
        .body_split_face_by_sheet(body, faces[0], tool_face, Some(&backend))
        .unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert_eq!(s.model().bodies().count(), count);
    assert_eq!(s.body_ask_faces(body).unwrap(), faces);
    assert_eq!(s.body_ask_edges(body).unwrap(), edges);
    assert!(s.body_check(body).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]
    #[test]
    fn rational_isoparametric_cuts_survive_random_sizes_and_frames(
        width in 2.0f64..4.,depth in 1.5f64..3.,weight in 0.7f64..1.3,
        fraction in 0.25f64..0.75,ox in -5.0f64..5.,ax in -0.3f64..0.3,ay in -0.3f64..0.3,
        multispan in any::<bool>(),
    ) {
        let frame=Frame::from_axis(Point3::new(ox,1.,2.),Vec3::new(ax,ay,1.),&Precision::DEFAULT).unwrap();
        cut(frame,width,depth,weight,fraction,multispan);
    }

    #[test]
    fn cuts_near_the_precision_boundary_leave_the_source_unchanged(factor in 0.01f64..0.9) {
        let p=Precision::DEFAULT;
        let backend=RustGeometryBackend::new(p);
        let mut s=Session::start(SessionOptions::default());
        let body=s.body_create_sheet_backend(&backend,&SurfaceConstruction::Extrude(profile(Frame::WORLD,3.,0.8,false),Vec3::Y*2.)).unwrap();
        let faces=s.body_ask_faces(body).unwrap();
        let edges=s.body_ask_edges(body).unwrap();
        let cutter=s.body_create_sheet_surface(&tool(Frame::WORLD,3.,2.,p.linear()*factor/2.,0.)).unwrap();
        let tool_face=s.body_ask_faces(cutter).unwrap()[0];
        let count=s.model().bodies().count();
        prop_assert!(s.body_split_face_by_sheet(body,faces[0],tool_face,Some(&backend)).is_err());
        prop_assert_eq!(s.model().bodies().count(),count);
        prop_assert_eq!(s.body_ask_faces(body).unwrap(),faces);
        prop_assert_eq!(s.body_ask_edges(body).unwrap(),edges);
        prop_assert!(s.body_check(body).unwrap().is_empty());
    }
}
