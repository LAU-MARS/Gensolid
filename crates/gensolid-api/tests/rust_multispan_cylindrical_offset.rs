//! Multiple circular spans and refined knots cross the real Session boundary.
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

fn joined_profile(radius: f64, angles: &[f64], parameters: &[f64]) -> BCurve {
    let mut points = Vec::new();
    let mut weights = Vec::new();
    let mut knots = vec![parameters[0]; 3];
    for i in 0..angles.len() - 1 {
        let piece = arc(
            Point3::ORIGIN,
            radius,
            [angles[i], angles[i + 1]],
            &Precision::DEFAULT,
        )
        .unwrap();
        assert_eq!(piece.points().len(), 3);
        points.extend_from_slice(&piece.points()[usize::from(i > 0)..]);
        weights.extend_from_slice(&piece.weights().unwrap()[usize::from(i > 0)..]);
        knots.extend(std::iter::repeat_n(
            parameters[i + 1],
            if i + 2 == angles.len() { 3 } else { 2 },
        ));
    }
    BCurve::new(2, knots, points, Some(weights)).unwrap()
}

fn placed(profile: &BCurve, frame: Frame, height: f64, reverse: bool) -> BSurface {
    let mut profile = BCurve::new(
        2,
        profile.knots().to_vec(),
        profile
            .points()
            .iter()
            .map(|p| frame.to_world(p.to_vec()))
            .collect(),
        profile.weights().map(<[f64]>::to_vec),
    )
    .unwrap();
    if reverse {
        profile = profile.reversed().unwrap();
    }
    let s = extrude_nurbs(&profile, frame.z() * height, &Precision::DEFAULT).unwrap();
    // V also need not be [0,1]. Its single Bernstein span already represents
    // differentiation in the full normalized V parameter.
    BSurface::new(
        2,
        1,
        s.knots().0.to_vec(),
        vec![-3.0, -3.0, 2.0, 2.0],
        s.counts().0,
        2,
        s.points().to_vec(),
        s.weights().map(<[f64]>::to_vec),
    )
    .unwrap()
}

fn changed(source: &BSurface, points: Vec<Point3>, weights: Option<Vec<f64>>) -> BSurface {
    BSurface::new(
        2,
        1,
        source.knots().0.to_vec(),
        source.knots().1.to_vec(),
        source.counts().0,
        2,
        points,
        weights,
    )
    .unwrap()
}

fn verify_session(source: &BSurface, frame: Frame, radius: f64, distance: f64, reverse: bool) {
    let p = Precision::DEFAULT;
    let original = source.clone();
    let candidate = offset_cylindrical(source, distance, &p).unwrap();
    assert_eq!(candidate.knots(), source.knots());
    assert_eq!(candidate.weights(), source.weights());
    certify_offset(source, &candidate, distance, &p).unwrap();
    let mut session = Session::start(SessionOptions::default());
    let input = session.body_create_sheet_surface(source).unwrap();
    let before = session
        .body_ask_vertices(input)
        .unwrap()
        .into_iter()
        .map(|v| session.vertex_ask_point(v).unwrap().map(f64::to_bits))
        .collect::<Vec<_>>();
    let output = session
        .body_offset_sheet_backend(input, &RustGeometryBackend::new(p), distance)
        .unwrap();
    assert!(session.body_check(input).unwrap().is_empty());
    assert!(session.body_check(output).unwrap().is_empty());
    let expected = radius + if reverse { -distance } else { distance };
    for vertex in session.body_ask_vertices(output).unwrap() {
        let [x, y, z] = session.vertex_ask_point(vertex).unwrap();
        let local = frame.to_local(Point3::new(x, y, z));
        assert!((local.x.hypot(local.y) - expected).abs() <= p.linear());
    }
    let after = session
        .body_ask_vertices(input)
        .unwrap()
        .into_iter()
        .map(|v| session.vertex_ask_point(v).unwrap().map(f64::to_bits))
        .collect::<Vec<_>>();
    assert_eq!(before, after);
    assert_eq!(*source, original);
}

#[test]
fn unequal_arc_angles_and_nonuniform_domains_preserve_the_actual_basis() {
    let p = Precision::DEFAULT;
    let frame =
        Frame::from_axis(Point3::new(2.0, -1.0, 3.0), Vec3::new(1.0, 2.0, 3.0), &p).unwrap();
    let profile = joined_profile(2.0, &[-1.0, -0.4, 0.3, 1.0], &[5.0, 6.0, 8.0, 11.0]);
    for reverse in [false, true] {
        let source = placed(&profile, frame, 1.5, reverse);
        certify_cylindrical_graph(&source, &p).unwrap();
        for distance in [-0.5, 0.0, 0.75] {
            verify_session(&source, frame, 2.0, distance, reverse);
        }
    }
}

#[test]
fn non_midpoint_refinement_with_unequal_local_endpoint_weights_is_supported() {
    let p = Precision::DEFAULT;
    let original = arc(Point3::ORIGIN, 2.0, [0.0, 1.2], &p).unwrap();
    let once = original.insert_knots(&[0.36]).unwrap();
    let first = once.trim([0.0, 0.36]).unwrap();
    assert_ne!(
        first.weights().unwrap()[0].to_bits(),
        first.weights().unwrap()[2].to_bits()
    );
    for profile in [once.clone(), once.insert_knots(&[0.36, 0.84]).unwrap()] {
        let source = placed(&profile, Frame::WORLD, 1.25, false);
        verify_session(&source, Frame::WORLD, 2.0, 0.35, false);
    }
}

struct LaterSpanAttack(bool);
impl GeometryBackend for LaterSpanAttack {
    fn precision(&self) -> Precision {
        Precision::DEFAULT
    }
    fn construct_surface(&self, _: &SurfaceConstruction) -> Result<BSurface> {
        Err(Error::unimplemented("offset-only test provider"))
    }
    fn intersect_surfaces(
        &self,
        _: &BSurface,
        _: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>> {
        Err(Error::unimplemented("offset-only test provider"))
    }
    fn offset_surface(&self, source: &BSurface, distance: f64) -> Result<BSurface> {
        let output = offset_cylindrical(source, distance, &Precision::DEFAULT)?;
        let mut points = output.points().to_vec();
        let i = points.len() - 4; // Last span's interior controls; all corners unchanged.
        if self.0 {
            points[i].z += Precision::DEFAULT.linear() * 64.0;
            points[i + 1].z -= Precision::DEFAULT.linear() * 64.0;
        } else {
            points[i] = source.points()[i];
            points[i + 1] = source.points()[i + 1];
        }
        Ok(changed(
            &output,
            points,
            output.weights().map(<[f64]>::to_vec),
        ))
    }
}

#[test]
fn actual_last_span_provider_corruption_is_rejected_without_publication() {
    let profile = joined_profile(2.0, &[-0.7, 0.0, 0.7], &[0.0, 0.4, 1.0]);
    let source = placed(&profile, Frame::WORLD, 1.0, false);
    for twist in [false, true] {
        let mut session = Session::start(SessionOptions::default());
        let input = session.body_create_sheet_surface(&source).unwrap();
        let count = session.model().bodies().count();
        let vertices = session.body_ask_vertices(input).unwrap();
        let before = vertices
            .iter()
            .map(|&v| session.vertex_ask_point(v).unwrap())
            .collect::<Vec<_>>();
        let error = session
            .body_offset_sheet_backend(input, &LaterSpanAttack(twist), 0.3)
            .unwrap_err();
        assert_eq!(error.code(), ErrorCode::Unimplemented);
        assert_eq!(count, session.model().bodies().count());
        assert_eq!(
            before,
            vertices
                .iter()
                .map(|&v| session.vertex_ask_point(v).unwrap())
                .collect::<Vec<_>>()
        );
        assert!(session.body_check(input).unwrap().is_empty());
    }
}

#[test]
fn individually_circular_spans_must_share_a_cylinder_and_direction() {
    let p = Precision::DEFAULT;
    let profile = joined_profile(2.0, &[-0.6, 0.0, 0.6], &[0.0, 0.5, 1.0]);
    let source = placed(&profile, Frame::WORLD, 1.0, false);
    // The second span is an exact circle of radius 1.5, tangent-continuous at
    // the shared point (2,0). It must not pass by certifying each center alone.
    let mut points = source.points().to_vec();
    for point in &mut points[6..] {
        point.x = 0.5 + point.x * 0.75;
        point.y *= 0.75;
    }
    let different_radius = changed(&source, points, source.weights().map(<[f64]>::to_vec));
    assert_eq!(
        offset_cylindrical(&different_radius, 0.25, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let mut session = Session::start(SessionOptions::default());
    let input = session
        .body_create_sheet_surface(&different_radius)
        .unwrap();
    let count = session.model().bodies().count();
    assert!(
        session
            .body_offset_sheet_backend(input, &RustGeometryBackend::new(p), 0.25)
            .is_err()
    );
    assert_eq!(count, session.model().bodies().count());
    assert!(session.body_check(input).unwrap().is_empty());
    for flip in [false, true] {
        let mut points = source.points().to_vec();
        for i in [7, 9] {
            if flip {
                points[i].z = -1.0;
            } else {
                points[i].x += 0.25;
            }
        }
        let different_axis = changed(&source, points, source.weights().map(<[f64]>::to_vec));
        assert_eq!(
            offset_cylindrical(&different_axis, 0.25, &p)
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
    }
    let mut points = source.points().to_vec();
    points[4].x += p.linear() * 64.0;
    points[5].x += p.linear() * 64.0;
    assert!(
        offset_cylindrical(
            &changed(&source, points, source.weights().map(<[f64]>::to_vec)),
            0.25,
            &p
        )
        .is_err()
    );
}

#[test]
fn local_graphs_do_not_admit_a_closed_circle_or_a_global_fold() {
    let p = Precision::DEFAULT;
    for angle in [3.5, std::f64::consts::TAU] {
        let profile = arc(Point3::ORIGIN, 2.0, [-angle / 2.0, angle / 2.0], &p).unwrap();
        let source = placed(&profile, Frame::WORLD, 1.0, false);
        assert!(certify_cylindrical_graph(&source, &p).is_err());
        assert!(offset_cylindrical(&source, 0.2, &p).is_err());
        let mut session = Session::start(SessionOptions::default());
        assert!(session.body_create_sheet_surface(&source).is_err());
        assert_eq!(session.model().bodies().count(), 0);
    }
}

#[test]
fn disconnected_full_multiplicity_seam_remains_unsupported() {
    let p = Precision::DEFAULT;
    let profile = joined_profile(2.0, &[-0.6, 0.0, 0.6], &[0.0, 0.5, 1.0])
        .insert_knots(&[0.5])
        .unwrap();
    let source = placed(&profile, Frame::WORLD, 1.0, false);
    assert_eq!(
        offset_cylindrical(&source, 0.2, &p).unwrap_err().code(),
        ErrorCode::Unimplemented
    );
    let mut session = Session::start(SessionOptions::default());
    assert!(session.body_create_sheet_surface(&source).is_err());
    assert_eq!(session.model().bodies().count(), 0);
}

#[test]
fn common_weight_scales_are_invariant_and_unbounded_denominators_reject_explicitly() {
    let p = Precision::DEFAULT;
    let profile = joined_profile(2.0, &[-0.6, 0.0, 0.6], &[0.0, 0.3, 1.0]);
    let source = placed(&profile, Frame::WORLD, 1.0, false);
    for factor in [f64::MIN_POSITIVE.sqrt(), f64::MAX.sqrt()] {
        let weights = source
            .weights()
            .unwrap()
            .iter()
            .map(|w| w * factor)
            .collect();
        let scaled = changed(&source, source.points().to_vec(), Some(weights));
        verify_session(&scaled, Frame::WORLD, 2.0, 0.3, false);
    }
    let mut weights = source.weights().unwrap().to_vec();
    weights[2] = f64::from_bits(1);
    weights[3] = f64::from_bits(1);
    let unbounded = changed(&source, source.points().to_vec(), Some(weights));
    let error = offset_cylindrical(&unbounded, 0.3, &p).unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert!(error.to_string().contains("denominator"));
    // A common scale can itself make the reciprocal unrepresentable; it must
    // be an explicit refusal, never a certificate based on an infinite bound.
    let weights = source
        .weights()
        .unwrap()
        .iter()
        .map(|_| f64::from_bits(1))
        .collect();
    let unbounded = changed(&source, source.points().to_vec(), Some(weights));
    let error = offset_cylindrical(&unbounded, 0.3, &p).unwrap_err();
    assert_eq!(error.code(), ErrorCode::Unimplemented);
    assert!(error.to_string().contains("unbounded"));
}

#[test]
fn a_small_positive_remaining_radius_is_validated_as_a_real_sheet() {
    let p = Precision::DEFAULT;
    let profile = joined_profile(2.0, &[-0.6, 0.0, 0.6], &[0.0, 0.5, 1.0]);
    let source = placed(&profile, Frame::WORLD, 1.0, false);
    verify_session(&source, Frame::WORLD, 2.0, -2.0 + p.linear() * 64.0, false);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]
    #[test]
    fn many_spans_are_checked_in_random_frames_scales_and_knot_domains(
        radius in 0.3f64..5.0,height in 0.4f64..3.0,angle in 0.7f64..1.4,
        count in 2usize..5,tilt in -3.0f64..3.0,translation in -4.0f64..4.0,
        start in -5.0f64..5.0,width in 0.5f64..8.0,
        fraction in -0.5f64..0.7,reverse in any::<bool>(),
    ) {
        let p=Precision::DEFAULT;
        let angles=(0..=count).map(|i|angle*(i as f64/count as f64-0.5)).collect::<Vec<_>>();
        let knots=(0..=count).map(|i|start+width*(i as f64/count as f64).powi(2)).collect::<Vec<_>>();
        let frame=Frame::from_axis(Point3::new(translation,-translation,translation),Vec3::new(tilt.sin(),tilt.cos(),1.0),&p)?;
        let profile=joined_profile(radius,&angles,&knots);
        let source=placed(&profile,frame,height,reverse);
        verify_session(&source,frame,radius,fraction*radius,reverse);
    }
    #[test]
    fn refined_simple_knots_and_near_collapse_preserve_session(
        radius in 0.3f64..5.0,angle in 0.7f64..1.4,
        split in 0.25f64..0.75,tilt in -2.0f64..2.0,factor in 0.1f64..0.7,
    ) {
        let p=Precision::DEFAULT;
        let frame=Frame::from_axis(Point3::new(tilt,-tilt,tilt),Vec3::new(tilt,1.0,1.0),&p)?;
        let profile=arc(Point3::ORIGIN,radius,[0.0,angle],&p)?.insert_knots(&[angle*split])?;
        let source=placed(&profile,frame,1.0,false);
        verify_session(&source,frame,radius,radius*0.25,false);
        let mut session=Session::start(SessionOptions::default());
        let input=session.body_create_sheet_surface(&source)?;
        let count=session.model().bodies().count();
        prop_assert!(session.body_offset_sheet_backend(input,&RustGeometryBackend::new(p),-radius+p.linear()*factor).is_err());
        prop_assert_eq!(session.model().bodies().count(),count);
        prop_assert!(session.body_check(input)?.is_empty());
    }
}
