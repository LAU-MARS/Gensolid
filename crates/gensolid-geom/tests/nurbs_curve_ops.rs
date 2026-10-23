//! Independent geometry identities and randomized exact-edit checks.
use gensolid_core::{ErrorCode, Precision, Result};
use gensolid_geom::BCurve;
use gensolid_geom::nurbs_curve_ops::{arc, circle, ellipse, line_segment};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

#[test]
fn rational_second_derivatives_match_the_closed_form() -> Result<()> {
    let p = Precision::DEFAULT;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let c = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![point(1., 0., 0.), point(1., 1., 0.), point(0., 1., 0.)],
        Some(vec![1., w, 1.]),
    )?;
    let e = c.eval2(0.5)?;
    assert!(e.point.distance(point(w, w, 0.)) <= p.linear());
    let d = 2. / (1. + w);
    assert!((e.d1 - Vec3::new(-d, d, 0.)).length() <= p.linear());
    assert!((e.d2 - Vec3::new(-2. * d * d * w, -2. * d * d * w, 0.)).length() <= p.linear());
    assert!(c.eval2(f64::NAN).is_err());
    assert!(c.eval2(f64::INFINITY).is_err());
    let shifted = BCurve::new(
        c.degree(),
        c.knots().to_vec(),
        c.points()
            .iter()
            .map(|&q| q + Vec3::new(p.size(), -p.size(), p.size()))
            .collect(),
        Some(
            vec![1., w, 1.]
                .into_iter()
                .map(|v| v * f64::MAX / 4.)
                .collect(),
        ),
    )?;
    let translated = shifted.eval2(0.5)?;
    assert!((translated.d2 - e.d2).length() <= p.linear());
    Ok(())
}

#[test]
#[allow(clippy::float_cmp)] // Knot endpoints must retain the supplied parameter bits.
fn repeated_knots_full_breaks_and_trims_preserve_parameter_values() -> Result<()> {
    let p = Precision::DEFAULT;
    let source = BCurve::new(
        2,
        vec![0., 0., 0., 0.5, 0.5, 1., 1., 1.],
        vec![
            point(0., 0., 0.),
            point(1., 2., 0.),
            point(2., 0., 0.),
            point(3., -2., 0.),
            point(4., 0., 0.),
        ],
        Some(vec![1., 2., 1., 0.5, 1.]),
    )?;
    let broken = source.insert_knots(&[0.5])?;
    let raised = broken.elevate_degree(5)?;
    let trimmed = raised.trim([0.2, 0.8])?;
    assert_eq!(trimmed.range().lo, 0.2);
    assert_eq!(trimmed.range().hi, 0.8);
    for t in [0.2, 0.3, 0.499, 0.5, 0.501, 0.7, 0.8] {
        let a = source.eval2(t)?;
        let b = trimmed.eval2(t)?;
        assert!(a.point.distance(b.point) <= p.linear());
        assert!((a.d1 - b.d1).length() <= p.linear());
        assert!((a.d2 - b.d2).length() <= p.linear());
    }
    assert!(broken.insert_knots(&[0.5]).is_err());
    assert!(source.trim([0.8, 0.2]).is_err());
    assert!(source.trim([f64::NAN, 0.8]).is_err());
    assert!(source.insert_knots(&[f64::NAN]).is_err());
    assert!(source.insert_knots(&[0.]).is_err());
    assert!(source.elevate_degree(1).is_err());
    assert_eq!(
        source.elevate_degree(65).err().map(|e| e.code()),
        Some(ErrorCode::CapacityExceeded)
    );
    Ok(())
}

#[test]
fn conics_satisfy_their_implicit_equations_and_endpoint_contracts() -> Result<()> {
    let p = Precision::DEFAULT;
    let center = point(3., -2., 1.);
    for radii in [[2., 2.], [3., 0.5]] {
        let c = ellipse(center, radii, &p)?;
        for i in 0..=120 {
            let t = std::f64::consts::TAU * i as f64 / 120.;
            let q = c.eval2(t)?.point - center;
            assert!(
                ((q.x / radii[0]).powi(2) + (q.y / radii[1]).powi(2) - 1.).abs() <= p.angular()
            );
            assert!(q.z.abs() <= p.linear());
        }
        assert!(
            c.eval(0.)
                .point
                .distance(c.eval(std::f64::consts::TAU).point)
                <= p.linear()
        );
    }
    let angles = [-0.7, 2.3];
    let c = arc(center, 2., angles, &p)?;
    for t in angles {
        let expected = center + Vec3::new(2. * t.cos(), 2. * t.sin(), 0.);
        assert!(c.eval(t).point.distance(expected) <= p.linear());
    }
    assert!(circle(center, p.linear() / 2., &p).is_err());
    assert!(line_segment(center, center + Vec3::new(p.linear() / 2., 0., 0.), &p).is_err());
    assert!(arc(center, 1., [0., p.angular() / 2.], &p).is_err());
    assert!(arc(center, 1., [0., std::f64::consts::TAU * 2.], &p).is_err());
    assert!(ellipse(center, [1., f64::NAN], &p).is_err());
    let large = p.size() / p.angular();
    let large_angles = [large, large + 1.0];
    let c = arc(center, 2., large_angles, &p)?;
    for t in large_angles {
        let expected = center + Vec3::new(2. * t.cos(), 2. * t.sin(), 0.);
        assert!(c.eval(t).point.distance(expected) <= p.linear());
    }
    Ok(())
}

#[test]
fn unclamped_edits_fail_explicitly_and_do_not_change_source() -> Result<()> {
    let c = BCurve::new(
        2,
        vec![-2., -1., 0., 1., 2., 3.],
        vec![point(0., 0., 0.), point(1., 1., 0.), point(2., 0., 0.)],
        None,
    )?;
    let before = c.clone();
    assert_eq!(
        c.insert_knots(&[0.5]).err().map(|e| e.code()),
        Some(ErrorCode::Unimplemented)
    );
    assert_eq!(
        c.trim([0.1, 0.9]).err().map(|e| e.code()),
        Some(ErrorCode::Unimplemented)
    );
    assert_eq!(
        c.elevate_degree(3).err().map(|e| e.code()),
        Some(ErrorCode::Unimplemented)
    );
    assert_eq!(c, before);
    assert!(c.eval2(0.5)?.point.is_finite());
    Ok(())
}

#[test]
fn nonunit_domains_and_disconnected_spans_keep_their_parameterization() -> Result<()> {
    let p = Precision::DEFAULT;
    let source = BCurve::new(
        2,
        vec![-3., -3., -3., 1., 1., 1., 5., 5., 5.],
        vec![
            point(0., 0., 0.),
            point(1., 1., 0.),
            point(2., 0., 0.),
            point(4., 2., 0.),
            point(5., 3., 0.),
            point(6., 2., 0.),
        ],
        Some(vec![1., 2., 1., 1., 0.5, 1.]),
    )?;
    let edited = source.insert_knots(&[-1., 3.])?.elevate_degree(4)?;
    let trimmed = edited.trim([-2., 4.])?;
    for t in [-2., -1., 0., 1., 2., 3., 4.] {
        let a = source.eval2(t)?;
        let b = trimmed.eval2(t)?;
        assert!(a.point.distance(b.point) <= p.linear());
        assert!((a.d1 - b.d1).length() <= p.linear());
        assert!((a.d2 - b.d2).length() <= p.linear());
    }
    assert!(
        source
            .eval(1. - p.angular())
            .point
            .distance(source.eval(1.).point)
            > p.linear()
    );
    assert!(
        trimmed
            .eval(1. - p.angular())
            .point
            .distance(trimmed.eval(1.).point)
            > p.linear()
    );
    let right = source.trim([1., 5.])?;
    let left = source.trim([-3., 1.])?;
    assert!(right.eval(1.).point.distance(point(4., 2., 0.)) <= p.linear());
    assert!(left.eval(1.).point.distance(point(2., 0., 0.)) <= p.linear());
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]
    #[test]
    fn rational_edits_preserve_points_and_derivatives_in_random_frames(
        size in 0.1f64..20., angle in -3.0f64..3., offset in -20.0f64..20.,
        weight in 0.25f64..4., knot in 0.1f64..0.9, t in 0.2f64..0.8,
        scale in -200i32..200
    ) {
        let p=Precision::DEFAULT;
        let frame=Frame::from_axis(point(offset,-offset,offset),Vec3::new(angle.cos(),angle.sin(),1.),&p)?;
        let controls=[Vec3::ZERO,Vec3::new(size,size/3.,0.),Vec3::new(2.*size,0.,size/2.)];
        let source=BCurve::new(2,vec![0.,0.,0.,1.,1.,1.],controls.into_iter().map(|q|frame.to_world(q)).collect(),Some(vec![1.,weight,1.].into_iter().map(|w|w*2_f64.powi(scale)).collect()))?;
        let before=source.clone();
        let edited=source.insert_knots(&[knot,knot])?.elevate_degree(4)?.trim([0.2,0.8])?;
        let a=source.eval2(t)?;let b=edited.eval2(t)?;
        prop_assert!(a.point.distance(b.point)<=p.linear());
        prop_assert!((a.d1-b.d1).length()<=p.linear());
        prop_assert!((a.d2-b.d2).length()<=p.linear());
        prop_assert_eq!(source,before);
    }

    #[test]
    fn near_precision_lines_reject_degeneracy_without_poisoning_valid_curves(factor in 0.01f64..0.99, offset in -10.0f64..10.) {
        let p=Precision::DEFAULT;let a=point(offset,0.,0.);
        prop_assert!(line_segment(a,a+Vec3::new(p.linear()*factor,0.,0.),&p).is_err());
        let c=line_segment(a,a+Vec3::new(p.linear()*4.,0.,0.),&p)?;
        prop_assert!(c.eval2(0.5)?.point.is_finite());
        prop_assert!(c.eval2(0.5)?.d2.length()<=p.linear());
    }
}
