//! End-to-end validation of the pure-Rust `SolidcoreBackend` inside checked `Session`
//! B-rep operations: incompatible-knot lofts, curved sweeps, planar/cylindrical offsets,
//! and single-span / multi-span rational oblique cuts.
#![allow(clippy::unwrap_used)]

use gensolid_api::{Session, SessionOptions, SheetCutOptions, SurfaceConstruction};
use gensolid_core::Precision;
use gensolid_geom::{BCurve, BSurface, solidcore_backend::SolidcoreBackend};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn plane(points: Vec<Point3>) -> BSurface {
    BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        2,
        2,
        points,
        None,
    )
    .unwrap()
}

#[test]
fn solidcore_incompatible_loft_and_curved_sweep_publish_checked_bodies() {
    let p = Precision::DEFAULT;
    let backend = SolidcoreBackend::new(p).unwrap();
    let mut session = Session::start(SessionOptions { precision: p });

    // Degree-1 two-point line lofted to a degree-2 three-point arch (incompatible knots/degrees).
    let base = BCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)],
        None,
    )
    .unwrap();
    let top = BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(1.5, 2.0, 0.35),
            Point3::new(3.0, 2.0, 0.0),
        ],
        None,
    )
    .unwrap();
    let loft_req = SurfaceConstruction::Loft(vec![base.clone(), top]);
    let loft_sheet = session
        .body_create_sheet_backend(&backend, &loft_req)
        .unwrap();
    assert!(session.body_check(loft_sheet).unwrap().is_empty());

    let loft_prism = session
        .body_create_solid_backend_prism(&backend, &loft_req, 1.0)
        .unwrap();
    assert_eq!(session.body_ask_faces(loft_prism).unwrap().len(), 6);
    assert!(session.body_check(loft_prism).unwrap().is_empty());

    // Moving-frame curved path sweep (`SurfaceConstruction::Sweep`), previously unimplemented
    // in `RustGeometryBackend`. In `solidcore`, the profile lies in local YZ (`y, z` cross-section
    // coordinates) and is placed along the trajectory's moving frame (`x` = start tangent).
    let sweep_profile = BCurve::new(
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![Point3::new(0.0, -1.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
        None,
    )
    .unwrap();
    let curved_path = BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.2, 0.15, 0.25),
            Point3::new(2.4, 0.0, 0.0),
        ],
        None,
    )
    .unwrap();
    let sweep_sheet = session
        .body_create_sheet_backend(
            &backend,
            &SurfaceConstruction::Sweep(sweep_profile, curved_path),
        )
        .unwrap();
    assert!(session.body_check(sweep_sheet).unwrap().is_empty());
}

#[test]
fn solidcore_cylindrical_offset_and_multispan_rational_oblique_cut_pass_checker() {
    let p = Precision::DEFAULT;
    let backend = SolidcoreBackend::new(p).unwrap();
    let mut session = Session::start(SessionOptions { precision: p });

    // Single-span rational quadratic circular arc extruded and cut obliquely.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let arc = BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 2.0),
            Point3::new(0.0, 0.0, 2.0),
        ],
        Some(vec![1.0, w, 1.0]),
    )
    .unwrap();
    let cyl_surface = BSurface::new(
        2,
        1,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        3,
        2,
        vec![
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
            Point3::new(2.0, 0.0, 2.0),
            Point3::new(2.0, 3.0, 2.0),
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(0.0, 3.0, 2.0),
        ],
        Some(vec![1.0, 1.0, w, w, 1.0, 1.0]),
    )
    .unwrap();
    let sheet = session.body_create_sheet_surface(&cyl_surface).unwrap();
    let offset_sheet = session
        .body_offset_sheet_backend(sheet, &backend, 0.25)
        .unwrap();
    assert!(session.body_check(offset_sheet).unwrap().is_empty());

    let offset_face = session.body_ask_faces(offset_sheet).unwrap()[0];
    let tool = session
        .body_create_sheet_surface(&plane(
            [(-2.0, -2.0), (-2.0, 5.0), (5.0, -2.0), (5.0, 5.0)]
                .into_iter()
                .map(|(x, z)| Point3::new(x, 1.4 + 0.1 * x + 0.08 * z, z))
                .collect(),
        ))
        .unwrap();
    let tool_face = session.body_ask_faces(tool).unwrap()[0];
    let split = session
        .body_split_face_by_sheet_with_options(
            offset_sheet,
            offset_face,
            tool_face,
            Some(&backend),
            SheetCutOptions {
                allow_rational_oblique: true,
            },
        )
        .unwrap();
    assert_eq!(session.body_ask_faces(split.body).unwrap().len(), 2);
    assert!(session.body_check(split.body).unwrap().is_empty());

    // Multi-span rational cylinder cut via SolidcoreBackend.
    let multi_arc = backend
        .import_curve(&arc)
        .unwrap()
        .insert_knots(&[0.5])
        .unwrap()
        .nurbs()
        .unwrap();
    let multi_cyl = BSurface::new(
        2,
        1,
        multi_arc.knots().to_vec(),
        vec![0.0, 0.0, 1.0, 1.0],
        multi_arc.points().len(),
        2,
        multi_arc
            .points()
            .iter()
            .flat_map(|&pt| [pt, pt + Vec3::new(0.0, 3.0, 0.0)])
            .collect(),
        multi_arc
            .weights()
            .map(|ws| ws.iter().flat_map(|&wt| [wt, wt]).collect()),
    )
    .unwrap();
    let multi_sheet = session.body_create_sheet_surface(&multi_cyl).unwrap();
    let multi_face = session.body_ask_faces(multi_sheet).unwrap()[0];
    let multi_split = session
        .body_split_face_by_sheet_multispan(multi_sheet, multi_face, tool_face, Some(&backend))
        .unwrap();
    assert_eq!(session.body_ask_faces(multi_split.body).unwrap().len(), 2);
    assert!(session.body_check(multi_split.body).unwrap().is_empty());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    #[test]
    fn solidcore_backend_extruded_prisms_and_cuts_pass_in_random_frames(
        width in 1.0f64..6.0,
        depth in 1.2f64..5.0,
        height in 0.6f64..3.0,
        arch in 0.05f64..0.35,
        cut_frac in 0.3f64..0.7,
        angle in -2.5f64..2.5,
        shift in -8.0f64..8.0,
    ) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(
            Point3::new(shift, -shift * 0.5, shift * 0.25),
            Vec3::new(angle.cos(), angle.sin(), 1.0),
            &p,
        )?;
        let backend = SolidcoreBackend::new(p)?;
        let mut session = Session::start(SessionOptions { precision: p });
        let prof = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                frame.to_world(Vec3::new(0.0, 0.0, 0.0)),
                frame.to_world(Vec3::new(width * 0.5, 0.0, arch * width)),
                frame.to_world(Vec3::new(width, 0.0, 0.0)),
            ],
            None,
        )?;
        let req = SurfaceConstruction::Extrude(prof, frame.y() * depth);
        let sheet = session.body_create_sheet_backend(&backend, &req)?;
        prop_assert!(session.body_check(sheet)?.is_empty());
        let sheet_face = session.body_ask_faces(sheet)?[0];

        let y_at = depth * cut_frac;
        let pad = width + height + 2.0;
        let cutter = plane(
            [(-pad, -pad), (-pad, pad), (pad, -pad), (pad, pad)]
                .into_iter()
                .map(|(x, z)| frame.to_world(Vec3::new(x, y_at, z)))
                .collect(),
        );
        let tool = session.body_create_sheet_surface(&cutter)?;
        let tool_face = session.body_ask_faces(tool)?[0];
        let cut = session.body_split_face_by_sheet(sheet, sheet_face, tool_face, Some(&backend))?;
        prop_assert_eq!(session.body_ask_faces(cut.body)?.len(), 2);
        prop_assert!(session.body_check(cut.body)?.is_empty());
    }
}
