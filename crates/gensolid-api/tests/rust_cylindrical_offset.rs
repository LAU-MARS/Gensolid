//! Cylindrical offsets pass the real Session/provider trust boundary.
#![allow(clippy::unwrap_used)]

use gensolid_api::{Session, SessionOptions, SurfaceConstruction};
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

fn cylinder(frame: Frame, radius: f64, height: f64, angle: f64, reversed: bool) -> BSurface {
    let p = Precision::DEFAULT;
    let circular = arc(Point3::ORIGIN, radius, [-angle / 2.0, angle / 2.0], &p).unwrap();
    let mut profile = BCurve::new(
        circular.degree(),
        circular.knots().to_vec(),
        circular
            .points()
            .iter()
            .map(|q| frame.to_world(q.to_vec()))
            .collect(),
        circular.weights().map(<[f64]>::to_vec),
    )
    .unwrap();
    if reversed {
        profile = profile.reversed().unwrap();
    }
    extrude_nurbs(&profile, frame.z() * height, &p).unwrap()
}

fn changed(surface: &BSurface, points: Vec<Point3>, weights: Option<Vec<f64>>) -> Result<BSurface> {
    let (uk, vk) = surface.knots();
    let (pu, pv) = surface.degrees();
    let (nu, nv) = surface.counts();
    BSurface::new(pu, pv, uk.to_vec(), vk.to_vec(), nu, nv, points, weights)
}

#[test]
fn small_rotated_circular_patch_has_a_complete_graph_certificate() {
    let p = Precision::DEFAULT;
    let frame = Frame::from_axis(Point3::ORIGIN, Vec3::new(0.0, 1.0, 1.0), &p).unwrap();
    let source = cylinder(frame, 0.2, 0.3, 0.3, false);
    // This exact case exhausted the existing general interval graph budget.
    // The specialized proof satisfies the same complete-domain inequalities.
    certify_cylindrical_graph(&source, &p).unwrap();
    let mut session = Session::start(SessionOptions::default());
    let input = session.body_create_sheet_surface(&source).unwrap();
    let output = session
        .body_offset_sheet_backend(input, &RustGeometryBackend::new(p), 0.0)
        .unwrap();
    assert!(session.body_check(input).unwrap().is_empty());
    assert!(session.body_check(output).unwrap().is_empty());
}

#[test]
fn positive_negative_zero_and_inward_orientation_offsets_enter_session() {
    let p = Precision::DEFAULT;
    for reversed in [false, true] {
        let source = cylinder(Frame::WORLD, 2.0, 1.5, 1.0, reversed);
        let original = source.clone();
        for distance in [-0.5, 0.0, 0.75] {
            let mut session = Session::start(SessionOptions::default());
            let body = session.body_create_sheet_surface(&source).unwrap();
            let offset = session
                .body_offset_sheet_backend(body, &RustGeometryBackend::new(p), distance)
                .unwrap();
            assert!(session.body_check(body).unwrap().is_empty());
            assert!(session.body_check(offset).unwrap().is_empty());
            let radius = 2.0 + if reversed { -distance } else { distance };
            for tag in session.body_ask_vertices(offset).unwrap() {
                let [x, y, _] = session.vertex_ask_point(tag).unwrap();
                assert!((x.hypot(y) - radius).abs() <= p.linear());
            }
            assert_eq!(source, original);
        }
    }
}

#[derive(Clone, Copy)]
enum Attack {
    Unchanged,
    WrongDistance,
    InteriorTwist,
    Basis,
    CrossAxis,
}
struct Incorrect(Attack);
impl GeometryBackend for Incorrect {
    fn precision(&self) -> Precision {
        Precision::DEFAULT
    }
    fn construct_surface(&self, _: &SurfaceConstruction) -> Result<BSurface> {
        Err(Error::unimplemented("test provider only offsets"))
    }
    fn intersect_surfaces(
        &self,
        _: &BSurface,
        _: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        Err(Error::unimplemented("test provider only offsets"))
    }
    fn offset_surface(&self, source: &BSurface, distance: f64) -> Result<BSurface> {
        let p = Precision::DEFAULT;
        match self.0 {
            Attack::Unchanged => Ok(source.clone()),
            Attack::WrongDistance => offset_cylindrical(source, distance * 0.5, &p),
            Attack::CrossAxis => changed(
                source,
                source
                    .points()
                    .iter()
                    .map(|q| Point3::new(-q.x, -q.y, q.z))
                    .collect(),
                source.weights().map(<[f64]>::to_vec),
            ),
            Attack::InteriorTwist => {
                let output = offset_cylindrical(source, distance, &p)?;
                let mut points = output.points().to_vec();
                // All four corners and the center are unchanged. A samples-only
                // certificate using those positions would miss this alteration.
                points[2].z += p.linear() * 64.0;
                points[3].z -= p.linear() * 64.0;
                changed(&output, points, output.weights().map(<[f64]>::to_vec))
            }
            Attack::Basis => {
                let output = offset_cylindrical(source, distance, &p)?;
                let mut weights = output.weights().unwrap().to_vec();
                weights[2] *= 2.0;
                changed(&output, output.points().to_vec(), Some(weights))
            }
        }
    }
}

#[test]
fn actual_provider_output_is_certified_and_rejection_does_not_publish() {
    let source = cylinder(Frame::WORLD, 2.0, 1.0, 1.0, false);
    for attack in [
        Attack::Unchanged,
        Attack::WrongDistance,
        Attack::InteriorTwist,
        Attack::Basis,
        Attack::CrossAxis,
    ] {
        let mut session = Session::start(SessionOptions::default());
        let body = session.body_create_sheet_surface(&source).unwrap();
        let vertices = session.body_ask_vertices(body).unwrap();
        let before = vertices
            .iter()
            .map(|&v| session.vertex_ask_point(v).unwrap().map(f64::to_bits))
            .collect::<Vec<_>>();
        let count = session.model().bodies().count();
        let distance = if matches!(attack, Attack::CrossAxis) {
            -4.0
        } else {
            0.5
        };
        let error = session
            .body_offset_sheet_backend(body, &Incorrect(attack), distance)
            .unwrap_err();
        assert_eq!(error.code(), ErrorCode::Unimplemented);
        assert_eq!(session.model().bodies().count(), count);
        let after = vertices
            .iter()
            .map(|&v| session.vertex_ask_point(v).unwrap().map(f64::to_bits))
            .collect::<Vec<_>>();
        assert_eq!(before, after);
        assert!(session.body_check(body).unwrap().is_empty());
    }
}

#[test]
fn noncircular_extrusion_and_collapsed_radius_remain_explicit() {
    let p = Precision::DEFAULT;
    let source = cylinder(Frame::WORLD, 2.0, 1.0, 1.0, false);
    for distance in [-2.0, -3.0, -2.0 + p.linear() * 0.5] {
        assert_eq!(
            offset_cylindrical(&source, distance, &p)
                .unwrap_err()
                .code(),
            ErrorCode::Degenerate
        );
    }
    let mut points = source.points().to_vec();
    for point in &mut points {
        point.x *= 1.5;
    }
    let ellipse = changed(&source, points, source.weights().map(<[f64]>::to_vec)).unwrap();
    assert_eq!(
        offset_cylindrical(&ellipse, 0.25, &p).unwrap_err().code(),
        ErrorCode::Unimplemented
    );
    assert!(offset_cylindrical(&source, f64::NAN, &p).is_err());
    assert!(certify_offset(&source, &source, f64::INFINITY, &p).is_err());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]
    #[test]
    fn circular_offsets_are_checked_in_random_frames_scales_and_orientations(
        radius in 0.2f64..6.0,
        height in 0.3f64..4.0,
        angle in 0.3f64..1.3,
        tilt in -3.0f64..3.0,
        translation in -5.0f64..5.0,
        fraction in -0.6f64..0.8,
        reversed in any::<bool>(),
    ) {
        let p = Precision::DEFAULT;
        let frame = Frame::from_axis(Point3::new(translation,-translation,translation),Vec3::new(tilt.sin(),tilt.cos(),1.0),&p)?;
        let source = cylinder(frame,radius,height,angle,reversed);
        let mut session = Session::start(SessionOptions::default());
        let input = session.body_create_sheet_surface(&source).map_err(|e|TestCaseError::fail(format!("source: {e}")))?;
        let distance = fraction*radius;
        let output = session.body_offset_sheet_backend(input,&RustGeometryBackend::new(p),distance).map_err(|e|TestCaseError::fail(format!("offset: {e}")))?;
        prop_assert!(session.body_check(input)?.is_empty());
        prop_assert!(session.body_check(output)?.is_empty());
        let expected = radius+if reversed {-distance} else {distance};
        for vertex in session.body_ask_vertices(output)? {
            let [x,y,z] = session.vertex_ask_point(vertex)?;
            let point = frame.to_local(Point3::new(x,y,z));
            prop_assert!((point.x.hypot(point.y)-expected).abs() <= p.linear());
        }
    }
    #[test]
    fn near_collapse_does_not_modify_session(radius in 0.2f64..6.0, factor in 0.1f64..0.8) {
        let p = Precision::DEFAULT;
        let source = cylinder(Frame::WORLD,radius,1.0,1.0,false);
        let mut session = Session::start(SessionOptions::default());
        let input = session.body_create_sheet_surface(&source)?;
        let count = session.model().bodies().count();
        prop_assert!(session.body_offset_sheet_backend(input,&RustGeometryBackend::new(p),-radius+p.linear()*factor).is_err());
        prop_assert_eq!(session.model().bodies().count(),count);
        prop_assert!(session.body_check(input)?.is_empty());
    }
}
