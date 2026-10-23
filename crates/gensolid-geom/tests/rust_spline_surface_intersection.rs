//! Exact rational paired curves, with independently evaluated source surfaces.
#![allow(clippy::unwrap_used)]
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::{
    BCurve, BSurface, rust_surface::extrude_nurbs,
    rust_surface_intersection::intersect_spline_extrusion_plane,
};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn profile(frame: Frame, multispan: bool, weight: f64) -> BCurve {
    let (knots, points, weights) = if multispan {
        (
            vec![2., 2., 2., 3.5, 3.5, 5., 5., 5.],
            vec![(0., 0.), (0.75, 0.4), (1.5, 0.2), (2.25, -0.2), (3., 0.)],
            vec![1., weight, 1., weight * 1.1, 1.],
        )
    } else {
        (
            vec![2., 2., 2., 5., 5., 5.],
            vec![(0., 0.), (1.5, 0.4), (3., 0.)],
            vec![1., weight, 1.],
        )
    };
    BCurve::new(
        2,
        knots,
        points
            .into_iter()
            .map(|(x, z)| frame.to_world(Vec3::new(x, 0., z)))
            .collect(),
        Some(weights),
    )
    .unwrap()
}

fn target(frame: Frame, multispan: bool, weight: f64, transpose: bool) -> BSurface {
    let s = extrude_nurbs(
        &profile(frame, multispan, weight),
        frame.y() * 2.,
        &Precision::DEFAULT,
    )
    .unwrap();
    let count = s.counts().0;
    if transpose {
        let indices: Vec<_> = (0..2)
            .flat_map(|j| (0..count).map(move |i| i * 2 + j))
            .collect();
        BSurface::new(
            1,
            2,
            vec![3., 3., 7., 7.],
            s.knots().0.to_vec(),
            2,
            count,
            indices.iter().map(|&i| s.points()[i]).collect(),
            Some(indices.iter().map(|&i| s.weights().unwrap()[i]).collect()),
        )
        .unwrap()
    } else {
        BSurface::new(
            2,
            1,
            s.knots().0.to_vec(),
            vec![3., 3., 7., 7.],
            count,
            2,
            s.points().to_vec(),
            s.weights().map(<[f64]>::to_vec),
        )
        .unwrap()
    }
}

fn plane(frame: Frame, offset: f64, slope: f64, bounded: bool) -> BSurface {
    let xs = if bounded { [-5., 8.] } else { [1., 2.] };
    BSurface::new(
        1,
        1,
        vec![-2., -2., 4., 4.],
        vec![8., 8., 12., 12.],
        2,
        2,
        xs.into_iter()
            .flat_map(|x| {
                [-5., 5.].map(move |z| frame.to_world(Vec3::new(x, offset + slope * (x + z), z)))
            })
            .collect(),
        None,
    )
    .unwrap()
}

fn assert_paired(a: &BSurface, b: &BSurface, spans: usize) {
    let p = Precision::DEFAULT;
    let hits = intersect_spline_extrusion_plane(a, b, &p).unwrap();
    assert_eq!(hits.len(), spans);
    for hit in hits {
        assert!(hit.reported_deviation <= p.linear());
        let r = hit.curve.range();
        assert_eq!(hit.first_uv.range(), r);
        assert_eq!(hit.second_uv.range(), r);
        for i in 0..=80 {
            let t = r.lo + (r.hi - r.lo) * (i as f64 / 80.);
            let xyz = hit.curve.eval(t).point;
            let uv = hit.first_uv.eval(t).point;
            let ab = hit.second_uv.eval(t).point;
            assert!(
                xyz.distance(a.eval(uv.x, uv.y).point) <= p.linear(),
                "first UV mismatch at {t}: {uv:?}"
            );
            assert!(
                xyz.distance(b.eval(ab.x, ab.y).point) <= p.linear(),
                "second UV mismatch at {t}"
            );
        }
    }
}

#[test]
fn rational_and_multiple_spans_have_paired_uv_in_both_orders_and_axes() {
    for multispan in [false, true] {
        for transpose in [false, true] {
            for slope in [0., 0.2] {
                let a = target(Frame::WORLD, multispan, 0.5, transpose);
                let b = plane(Frame::WORLD, 0.7, slope, true);
                assert_paired(&a, &b, 1 + usize::from(multispan));
                assert_paired(&b, &a, 1 + usize::from(multispan));
            }
        }
    }
}

#[test]
fn rational_profile_uv_reproduces_parameter_not_rational_greville_average() {
    let a = target(Frame::WORLD, false, 0.5, false);
    let b = plane(Frame::WORLD, 0.7, 0.2, true);
    let hit = intersect_spline_extrusion_plane(&a, &b, &Precision::DEFAULT)
        .unwrap()
        .remove(0);
    assert_eq!(hit.curve.degree(), 2);
    assert_eq!(hit.first_uv.degree(), 3);
    for i in 0..=64 {
        let t = 2. + 3. * (i as f64 / 64.);
        assert!((hit.first_uv.eval(t).point.x - t).abs() <= Precision::DEFAULT.solver());
    }
    let wrong = BCurve::new(
        2,
        vec![2., 2., 2., 5., 5., 5.],
        vec![
            Point3::new(2., 0., 0.),
            Point3::new(3.5, 0., 0.),
            Point3::new(5., 0., 0.),
        ],
        Some(vec![1., 0.5, 1.]),
    )
    .unwrap();
    assert!((wrong.eval(2.75).point.x - 2.75).abs() > Precision::DEFAULT.linear());
}

#[test]
fn disjoint_is_proved_and_partial_overlap_and_ruling_weight_changes_are_refused() {
    let p = Precision::DEFAULT;
    let a = target(Frame::WORLD, true, 0.5, false);
    assert!(
        intersect_spline_extrusion_plane(&a, &plane(Frame::WORLD, 4., 0., true), &p)
            .unwrap()
            .is_empty()
    );
    for b in [
        plane(Frame::WORLD, 0.7, 0., false),
        plane(Frame::WORLD, -0.1, 0.2, true),
    ] {
        assert_eq!(
            intersect_spline_extrusion_plane(&a, &b, &p)
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
    }
    assert_eq!(
        intersect_spline_extrusion_plane(&a, &a, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let mut weights = a.weights().unwrap().to_vec();
    weights[1] *= 2.;
    let changed = BSurface::new(
        2,
        1,
        a.knots().0.to_vec(),
        a.knots().1.to_vec(),
        a.counts().0,
        2,
        a.points().to_vec(),
        Some(weights),
    )
    .unwrap();
    assert_eq!(
        intersect_spline_extrusion_plane(&changed, &plane(Frame::WORLD, 0.7, 0., true), &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}

#[test]
fn large_translated_uv_domains_are_refused_even_if_control_uv_is_representable() {
    let p = Precision::DEFAULT;
    let a = target(Frame::WORLD, false, 0.5, false);
    let b = plane(Frame::WORLD, 0.5, 0., true);
    let lo = 1_000_000_000_000_000.;
    for (u, v) in [
        (
            vec![lo, lo, lo, lo + 1., lo + 1., lo + 1.],
            a.knots().1.to_vec(),
        ),
        (a.knots().0.to_vec(), vec![lo, lo, lo + 1., lo + 1.]),
    ] {
        let shifted = BSurface::new(
            2,
            1,
            u,
            v,
            3,
            2,
            a.points().to_vec(),
            a.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        assert_eq!(
            intersect_spline_extrusion_plane(&shifted, &b, &p)
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
    }
    let tool = BSurface::new(
        1,
        1,
        vec![lo, lo, lo + 1., lo + 1.],
        vec![lo, lo, lo + 1., lo + 1.],
        2,
        2,
        b.points().to_vec(),
        None,
    )
    .unwrap();
    assert_eq!(
        intersect_spline_extrusion_plane(&a, &tool, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}

#[test]
fn discontinuous_profiles_and_unresolved_weight_conditioning_are_refused() {
    let p = Precision::DEFAULT;
    let b = plane(Frame::WORLD, 0.7, 0., true);
    let broken = BCurve::new(
        2,
        vec![0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1.],
        [
            (0., 0.),
            (0.5, 0.1),
            (1., 0.),
            (2., 0.),
            (2.5, 0.1),
            (3., 0.),
        ]
        .into_iter()
        .map(|(x, z)| Point3::new(x, 0., z))
        .collect(),
        Some(vec![1., 0.8, 1., 1., 0.8, 1.]),
    )
    .unwrap();
    let a = extrude_nurbs(&broken, Vec3::Y * 2., &p).unwrap();
    assert_eq!(
        intersect_spline_extrusion_plane(&a, &b, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let a = target(Frame::WORLD, false, 1_000_000_000_000., false);
    assert_eq!(
        intersect_spline_extrusion_plane(&a, &b, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]
    #[test]
    fn rational_section_equations_hold_in_random_frames(
        weight in 0.4f64..2., ox in -10.0f64..10., ax in -0.4f64..0.4,
        ay in -0.4f64..0.4, slope in -0.05f64..0.15, multispan in any::<bool>(), transpose in any::<bool>(),
    ) {
        let p=Precision::DEFAULT;
        let frame=Frame::from_axis(Point3::new(ox,1.,2.),Vec3::new(ax,ay,1.),&p).unwrap();
        let a=target(frame,multispan,weight,transpose);
        let b=plane(frame,0.8,slope,true);
        assert_paired(&a,&b,1+usize::from(multispan));
        assert_paired(&b,&a,1+usize::from(multispan));
    }
}
