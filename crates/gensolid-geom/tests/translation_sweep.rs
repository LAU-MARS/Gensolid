//! Independent curve-sum and derivative identities for fixed-orientation sweeps.
use gensolid_core::{ErrorCode, Precision, Result};
use gensolid_geom::{BCurve, translation_sweep::sweep_translation};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn curve(
    degree: usize,
    knots: Vec<f64>,
    controls: &[[f64; 3]],
    weights: Vec<f64>,
) -> Result<BCurve> {
    BCurve::new(
        degree,
        knots,
        controls
            .iter()
            .map(|q| Point3::new(q[0], q[1], q[2]))
            .collect(),
        Some(weights),
    )
}

#[test]
fn unclamped_path_uses_its_evaluated_start_and_ignores_common_translation() -> Result<()> {
    let p = Precision::DEFAULT;
    let profile = curve(
        2,
        vec![2., 2., 2., 5., 5., 5.],
        &[[0., 0., 0.], [2., 0., 1.], [4., 0., 0.]],
        vec![1., 2., 1.],
    )?;
    let path = curve(
        2,
        vec![-2., -1., 0., 1., 2., 3., 4.],
        &[[0., 0., 0.], [1., 2., 1.], [-1., 4., 2.], [0., 6., 0.]],
        vec![1., 2., 1., 2.],
    )?;
    let anchor = path.eval(path.range().lo).point;
    assert!(anchor.distance(path.points()[0]) > p.linear());
    let swept = sweep_translation(&profile, &path, &p)?;
    assert_eq!(swept.knots(), (profile.knots(), path.knots()));
    for u in [2., 2.7, 3.6, 5.] {
        for v in [0., 0.2, 0.8, 1., 1.7, 2.] {
            let actual = swept.eval(u, v);
            let c = profile.eval(u);
            let q = path.eval(v);
            assert!(actual.point.distance(c.point + (q.point - anchor)) <= p.linear());
            assert!((actual.du - c.d1).length() <= p.linear());
            assert!((actual.dv - q.d1).length() <= p.linear());
        }
    }
    let shifted = BCurve::new(
        path.degree(),
        path.knots().to_vec(),
        path.points()
            .iter()
            .map(|&q| q + Vec3::new(2_f64.powi(50), -2_f64.powi(50), 2_f64.powi(50)))
            .collect(),
        path.weights().map(|w| w.to_vec()),
    )?;
    let shifted_surface = sweep_translation(&profile, &shifted, &p)?;
    for u in [2., 3., 5.] {
        for v in [0., 0.25, 0.75, 1.5, 2.] {
            assert!(
                swept
                    .eval(u, v)
                    .point
                    .distance(shifted_surface.eval(u, v).point)
                    <= p.linear()
            );
        }
    }
    Ok(())
}

#[test]
fn common_weight_scales_cancel_and_underflow_is_explicit() -> Result<()> {
    let p = Precision::DEFAULT;
    let make = |scale: f64| {
        curve(
            1,
            vec![0., 0., 1., 1.],
            &[[0., 0., 0.], [2., 0., 1.]],
            vec![scale, scale * 2.],
        )
    };
    let path = curve(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        &[[0., 0., 0.], [0., 2., 1.], [0., 4., 0.]],
        vec![1., 2., 1.],
    )?;
    let baseline = sweep_translation(&make(1.)?, &path, &p)?;
    for scale in [2_f64.powi(-700), 2_f64.powi(700)] {
        let scaled_path = BCurve::new(
            path.degree(),
            path.knots().to_vec(),
            path.points().to_vec(),
            path.weights()
                .map(|w| w.iter().map(|q| q * scale).collect()),
        )?;
        let actual = sweep_translation(&make(scale)?, &scaled_path, &p)?;
        for u in [0., 0.2, 0.7, 1.] {
            for v in [0., 0.3, 0.8, 1.] {
                assert!(actual.eval(u, v).point.distance(baseline.eval(u, v).point) <= p.linear());
            }
        }
    }
    let tiny = curve(
        1,
        vec![0., 0., 1., 1.],
        &[[0., 0., 0.], [2., 0., 1.]],
        vec![f64::MIN_POSITIVE, 1.],
    )?;
    assert_eq!(
        sweep_translation(&tiny, &tiny, &p).err().map(|e| e.code()),
        Some(ErrorCode::Degenerate)
    );
    let extreme = curve(
        1,
        vec![0., 0., 1., 1.],
        &[[0., 0., 0.], [2., 0., 1.]],
        vec![f64::MIN_POSITIVE, f64::MAX],
    )?;
    assert_eq!(
        sweep_translation(&extreme, &path, &p)
            .err()
            .map(|e| e.code()),
        Some(ErrorCode::Degenerate)
    );
    let subnormal = curve(
        1,
        vec![0., 0., 1., 1.],
        &[[0., 0., 0.], [2., 0., 1.]],
        vec![f64::MIN_POSITIVE / 2., 1.],
    )?;
    assert_eq!(
        sweep_translation(&subnormal, &path, &p)
            .err()
            .map(|e| e.code()),
        Some(ErrorCode::Degenerate)
    );
    let first = curve(
        1,
        vec![0., 0., 1., 1.],
        &[[0., 0., 0.], [2., 0., 1.]],
        vec![2_f64.powi(-600), 1.],
    )?;
    let second = curve(
        1,
        vec![0., 0., 1., 1.],
        &[[0., 0., 0.], [0., 2., 1.]],
        vec![2_f64.powi(-470), 1.],
    )?;
    assert_eq!(
        sweep_translation(&first, &second, &p)
            .err()
            .map(|e| e.code()),
        Some(ErrorCode::Degenerate)
    );
    Ok(())
}

#[test]
fn large_profile_hulls_refuse_interior_quantization_even_with_exact_controls() -> Result<()> {
    let p = Precision::DEFAULT;
    let base = 2_f64.powi(50);
    let profile = curve(
        1,
        vec![0., 0., 1., 1.],
        &[[base, 0., 0.], [base + 1., 0., 0.]],
        vec![1., 1.],
    )?;
    // 0.2 is lost/quantized in control addition; 0.5 gives exactly representable
    // controls but still cannot resolve arbitrary interior sweep evaluations.
    for dx in [0.2, 0.5] {
        let path = curve(
            1,
            vec![0., 0., 1., 1.],
            &[[0., 0., 0.], [dx, 1., 0.]],
            vec![1., 1.],
        )?;
        assert_eq!(
            sweep_translation(&profile, &path, &p)
                .err()
                .map(|e| e.code()),
            Some(ErrorCode::Unimplemented)
        );
    }
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]
    #[test]
    fn rational_spatial_path_preserves_sum_and_partials_in_random_frames(
        size in 0.5f64..10., tilt in -1.0f64..1., offset in -10.0f64..10.,
        w in 0.5f64..2., q in 0.5f64..2., u in 2.0f64..5., v in -1.0f64..2.
    ) {
        let p=Precision::DEFAULT;
        let frame=Frame::from_axis(Point3::new(offset,-offset,offset),Vec3::new(tilt,1.,2.),&p)?;
        let profile=BCurve::new(2,vec![2.,2.,2.,5.,5.,5.],
            [Vec3::ZERO,Vec3::new(size/2.,0.,size/8.),Vec3::new(size,0.,0.)].map(|x|frame.to_world(x)).to_vec(),Some(vec![1.,w,1.]))?;
        // Four non-coplanar controls define an actual spatial cubic path.
        let path=BCurve::new(3,vec![-1.,-1.,-1.,-1.,2.,2.,2.,2.],
            [Vec3::ZERO,Vec3::new(size/8.,size/3.,size/9.),Vec3::new(-size/8.,2.*size/3.,size/7.),Vec3::new(0.,size,0.)].map(|x|frame.to_world(x)).to_vec(),Some(vec![1.,q,w,1.]))?;
        let source_profile=profile.clone();let source_path=path.clone();
        let surface=sweep_translation(&profile,&path,&p)?;
        let actual=surface.eval(u,v);let a=profile.eval(u);let b=path.eval(v);
        prop_assert!(actual.point.distance(a.point+(b.point-path.eval(-1.).point))<=p.linear());
        prop_assert!((actual.du-a.d1).length()<=p.linear());
        prop_assert!((actual.dv-b.d1).length()<=p.linear());
        prop_assert_eq!(profile,source_profile);prop_assert_eq!(path,source_path);
    }

    #[test]
    fn subprecision_paths_fail_explicitly(factor in 0.01f64..0.99) {
        let p=Precision::DEFAULT;
        let profile=curve(1,vec![0.,0.,1.,1.],&[[0.,0.,0.],[2.,0.,0.]],vec![1.,1.])?;
        let path=curve(2,vec![0.,0.,0.,1.,1.,1.],&[[0.,0.,0.],[0.,p.linear()*factor/2.,0.],[0.,p.linear()*factor,0.]],vec![1.,1.,1.])?;
        prop_assert_eq!(sweep_translation(&profile,&path,&p).err().map(|e|e.code()),Some(ErrorCode::Degenerate));
    }
}
