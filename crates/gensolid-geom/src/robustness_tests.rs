//! Geometry regressions with independently known points, derivatives and projections.

use std::f64::consts::FRAC_1_SQRT_2;

use gensolid_core::Precision;
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

use crate::{BCurve, BSurface, Curve, GeomStore, ProjectionKind, Surface};

fn frame_strategy() -> impl Strategy<Value = Frame> {
    (
        -100.0f64..100.0,
        -100.0f64..100.0,
        -100.0f64..100.0,
        -1.0f64..1.0,
        -1.0f64..1.0,
        0.1f64..1.0,
    )
        .prop_map(|(x, y, z, ax, ay, az)| {
            Frame::from_axis(
                Point3::new(x, y, z),
                Vec3::new(ax, ay, az),
                &Precision::DEFAULT,
            )
            .unwrap()
        })
}

fn quarter_circle(frame: Frame, radius: f64, knot_scale: f64, weight_scale: f64) -> BCurve {
    BCurve::new(
        2,
        vec![0.0, 0.0, 0.0, knot_scale, knot_scale, knot_scale],
        vec![
            frame.to_world(Vec3::new(radius, 0.0, 0.0)),
            frame.to_world(Vec3::new(radius, radius, 0.0)),
            frame.to_world(Vec3::new(0.0, radius, 0.0)),
        ],
        Some(vec![
            weight_scale,
            FRAC_1_SQRT_2 * weight_scale,
            weight_scale,
        ]),
    )
    .unwrap()
}

#[test]
fn full_multiplicity_knot_retains_both_one_sided_values() {
    let p = Precision::DEFAULT;
    let curve = BCurve::new(
        1,
        vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
        ],
        None,
    )
    .unwrap();
    assert!(
        curve
            .eval(0.5f64.next_down())
            .point
            .distance(Point3::new(1.0, 0.0, 0.0))
            <= p.linear()
    );
    assert!(curve.eval(0.5).point.distance(Point3::new(3.0, 0.0, 0.0)) <= p.linear());
    for t in [0.0, 0.5f64.next_down(), 0.5, 1.0] {
        assert!((curve.eval(t).d1 - Vec3::new(2.0, 0.0, 0.0)).length() <= p.linear());
    }
    let mut store = GeomStore::new();
    let id = store.add_curve(Curve::BCurve(curve)).unwrap();
    for x in [0.25, 0.99, 3.01, 3.75] {
        let projected = store
            .project_to_curve(id, Point3::new(x, 0.5, 0.0), &p)
            .unwrap();
        assert!((projected.distance - 0.5).abs() <= p.linear());
        assert!((projected.point.x - x).abs() <= p.linear());
    }
}

proptest! {
    #[test]
    fn rational_curve_is_invariant_under_weight_and_parameter_scales(
        frame in frame_strategy(),
        log_radius in -5.0f64..2.0,
        log_knot_scale in -6.0f64..3.0,
        log_weight_scale in -300.0f64..300.0,
        t in 0.01f64..0.99,
    ) {
        let p = Precision::DEFAULT;
        let radius = 10.0f64.powf(log_radius);
        let knot_scale = 10.0f64.powf(log_knot_scale);
        let curve = quarter_circle(frame, radius, knot_scale, 10.0f64.powf(log_weight_scale));
        let base = quarter_circle(frame, radius, 1.0, 1.0);
        let value = curve.eval(t * knot_scale);
        let expected = base.eval(t);
        prop_assert!(value.point.distance(expected.point) <= p.linear());
        prop_assert!((value.d1 * knot_scale - expected.d1).length() <= p.linear() * (1.0 + expected.d1.length()));
        prop_assert!((value.point.distance(frame.origin()) - radius).abs() <= p.linear());
        let mut store = GeomStore::new();
        let id = store.add_curve(Curve::BCurve(curve)).unwrap();
        let query = expected.point + frame.z() * radius;
        let projected = store.project_to_curve(id, query, &p).unwrap();
        prop_assert_eq!(projected.kind, ProjectionKind::NumericalCandidate);
        prop_assert!(projected.point.distance(expected.point) <= p.linear());
        prop_assert!((projected.distance - radius).abs() <= p.linear());
    }

    #[test]
    fn rational_surface_derivatives_and_projection_work_in_random_frames(
        frame in frame_strategy(),
        log_scale in -5.0f64..2.0,
        u in 0.05f64..0.95,
        v in 0.05f64..0.95,
        weights in prop::array::uniform4(0.2f64..5.0),
        offset_factor in -0.5f64..0.5,
    ) {
        let p = Precision::DEFAULT;
        let size = 10.0f64.powf(log_scale);
        let surface = BSurface::new(1, 1, vec![0.0, 0.0, 1.0, 1.0], vec![0.0, 0.0, 1.0, 1.0], 2, 2,
            vec![frame.origin(), frame.origin() + frame.y() * size,
                 frame.origin() + frame.x() * size, frame.origin() + (frame.x() + frame.y()) * size],
            Some(weights.to_vec())).unwrap();
        let value = surface.eval(u, v);
        let h = p.linear().sqrt();
        let du = (surface.eval(u + h, v).point - surface.eval(u - h, v).point) / (2.0 * h);
        let dv = (surface.eval(u, v + h).point - surface.eval(u, v - h).point) / (2.0 * h);
        prop_assert!((du - value.du).length() <= h * (1.0 + value.du.length()));
        prop_assert!((dv - value.dv).length() <= h * (1.0 + value.dv.length()));
        let offset = size * offset_factor;
        let query = value.point + frame.z() * offset;
        let mut store = GeomStore::new();
        let id = store.add_surface(Surface::BSurface(surface)).unwrap();
        let projected = store.project_to_surface(id, query, &p).unwrap();
        prop_assert_eq!(projected.kind, ProjectionKind::NumericalCandidate);
        prop_assert!(projected.point.distance(value.point) <= p.linear());
        prop_assert!((projected.distance - offset.abs()).abs() <= p.linear());
    }
}
