//! Fixed-orientation curved-path sweeps enter unchanged complete body validation.
#![allow(clippy::unwrap_used)]
use gensolid_api::{RustGeometryBackend, Session, SessionOptions, SurfaceConstruction, Tag};
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::BCurve;
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn curves(frame: Frame, width: f64, depth: f64, bend: f64, w: f64) -> (BCurve, BCurve) {
    let profile = BCurve::new(
        2,
        vec![2., 2., 2., 5., 5., 5.],
        [
            Vec3::ZERO,
            Vec3::new(width / 2., 0., bend),
            Vec3::new(width, 0., 0.),
        ]
        .map(|q| frame.to_world(q))
        .to_vec(),
        Some(vec![1., w, 1.]),
    )
    .unwrap();
    let path = BCurve::new(
        3,
        vec![-1., -1., -1., -1., 2., 2., 2., 2.],
        [
            Vec3::ZERO,
            Vec3::new(depth / 30., depth / 3., bend / 2.),
            Vec3::new(-depth / 30., 2. * depth / 3., bend),
            Vec3::new(0., depth, 0.),
        ]
        .map(|q| frame.to_world(q))
        .to_vec(),
        Some(vec![1., w, 1., 1.]),
    )
    .unwrap();
    (profile, path)
}

#[derive(Debug, PartialEq, Eq)]
struct Snapshot {
    faces: Vec<Tag>,
    edges: Vec<Tag>,
    vertices: Vec<(Tag, [u64; 3])>,
}
fn snapshot(s: &mut Session, body: Tag) -> Snapshot {
    Snapshot {
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
fn curved_spatial_path_creates_owned_checked_sheet_and_solid() {
    let p = Precision::DEFAULT;
    let (profile, path) = curves(Frame::WORLD, 3., 2., 0.2, 1.1);
    let request = SurfaceConstruction::TranslationSweep { profile, path };
    let mut session = Session::start(SessionOptions::default());
    let (sheet, solid) = {
        let backend = RustGeometryBackend::new(p);
        (
            session
                .body_create_sheet_backend(&backend, &request)
                .unwrap(),
            session
                .body_create_solid_backend_prism(&backend, &request, 0.5)
                .unwrap(),
        )
    };
    assert!(session.body_check(sheet).unwrap().is_empty());
    assert!(session.body_check(solid).unwrap().is_empty());
    assert_eq!(session.body_ask_faces(sheet).unwrap().len(), 1);
    assert_eq!(session.body_ask_faces(solid).unwrap().len(), 6);
    assert_eq!(session.body_ask_edges(solid).unwrap().len(), 12);
    assert_eq!(session.body_ask_vertices(solid).unwrap().len(), 8);
}

#[test]
fn singular_and_discontinuous_sweeps_do_not_publish_partial_state() {
    let p = Precision::DEFAULT;
    let backend = RustGeometryBackend::new(p);
    let mut session = Session::start(SessionOptions::default());
    let existing = session.body_create_solid_block(2., 3., 4., None).unwrap();
    let before = snapshot(&mut session, existing);
    let count = session.model().bodies().count();
    let (profile, _) = curves(Frame::WORLD, 3., 2., 0.2, 1.);
    let legacy = SurfaceConstruction::Sweep(profile.clone(), profile.clone());
    assert_eq!(
        session
            .body_create_sheet_backend(&backend, &legacy)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let singular = SurfaceConstruction::TranslationSweep {
        profile: profile.clone(),
        path: profile.clone(),
    };
    assert!(
        session
            .body_create_sheet_backend(&backend, &singular)
            .is_err()
    );
    assert!(
        session
            .body_create_solid_backend_prism(&backend, &singular, 0.5)
            .is_err()
    );
    let discontinuous = BCurve::new(
        1,
        vec![0., 0., 0.5, 0.5, 1., 1.],
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(0., 1., 0.),
            Point3::new(0., 2., 0.),
            Point3::new(0., 3., 0.),
        ],
        None,
    )
    .unwrap();
    let request = SurfaceConstruction::TranslationSweep {
        profile,
        path: discontinuous,
    };
    assert!(
        session
            .body_create_sheet_backend(&backend, &request)
            .is_err()
    );
    assert_eq!(session.model().bodies().count(), count);
    assert_eq!(snapshot(&mut session, existing), before);
    assert!(session.body_check(existing).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]
    #[test]
    fn random_curved_path_prisms_pass_the_full_checker(
        width in 1.0f64..4.,depth in 1.0f64..4.,height in 0.2f64..1.,
        bend in 0.02f64..0.12,weight in 0.9f64..1.1,offset in -5.0f64..5.,tilt in -0.3f64..0.3
    ) {
        let p=Precision::DEFAULT;let frame=Frame::from_axis(Point3::new(offset,1.,2.),Vec3::new(tilt,0.2,1.),&p).unwrap();
        let (profile,path)=curves(frame,width,depth,bend,weight);
        let request=SurfaceConstruction::TranslationSweep {profile,path};
        let backend=RustGeometryBackend::new(p);let mut s=Session::start(SessionOptions::default());
        let existing=s.body_create_solid_block(2.,2.,2.,None).unwrap();let before=snapshot(&mut s,existing);
        let sheet=s.body_create_sheet_backend(&backend,&request).unwrap();
        let solid=s.body_create_solid_backend_prism(&backend,&request,height).unwrap();
        prop_assert!(s.body_check(sheet).unwrap().is_empty());prop_assert!(s.body_check(solid).unwrap().is_empty());
        prop_assert_eq!(s.body_ask_faces(solid).unwrap().len(),6);prop_assert_eq!(snapshot(&mut s,existing),before);
    }

    #[test]
    fn near_degenerate_paths_leave_existing_bodies_unchanged(factor in 0.01f64..0.99) {
        let p=Precision::DEFAULT;let backend=RustGeometryBackend::new(p);let mut s=Session::start(SessionOptions::default());
        let existing=s.body_create_solid_block(1.,2.,3.,None).unwrap();let before=snapshot(&mut s,existing);let count=s.model().bodies().count();
        let (profile,_)=curves(Frame::WORLD,3.,2.,0.2,1.);
        let path=BCurve::new(2,vec![0.,0.,0.,1.,1.,1.],vec![Point3::ORIGIN,Point3::new(0.,p.linear()*factor/2.,0.),Point3::new(0.,p.linear()*factor,0.)],None).unwrap();
        let request=SurfaceConstruction::TranslationSweep {profile,path};
        prop_assert_eq!(s.body_create_sheet_backend(&backend,&request).unwrap_err().code(),ErrorCode::Degenerate);
        prop_assert!(s.body_create_solid_backend_prism(&backend,&request,0.5).is_err());
        prop_assert_eq!(s.model().bodies().count(),count);prop_assert_eq!(snapshot(&mut s,existing),before);
    }
}
