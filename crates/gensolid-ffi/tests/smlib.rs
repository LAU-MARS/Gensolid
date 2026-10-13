//! Integration and property tests against the real native SMLib library.
#![cfg(feature = "smlib")]

use gensolid::smlib::{Backend, Direction, Parameterization};
use gensolid_core::{ErrorCode, Precision, Result};
use gensolid_geom::{BCurve, BSurface};
use gensolid_math::{Frame, Point3, Vec3};
use proptest::prelude::*;

fn point(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}
fn line(a: Point3, b: Point3) -> Result<BCurve> {
    BCurve::new(1, vec![0., 0., 1., 1.], vec![a, b], None)
}
fn plane(z: f64) -> Result<BSurface> {
    BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            point(0., 0., z),
            point(0., 1., z),
            point(1., 0., z),
            point(1., 1., z),
        ],
        None,
    )
}
fn close(a: Point3, b: Point3, p: Precision) {
    assert!(a.distance(b) <= p.linear(), "{a:?} != {b:?}");
}

#[test]
fn rational_roundtrip_evaluation_and_owned_lifetime() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let circle = b.circle(point(0., 0., 0.), 2.)?;
    let data = circle.nurbs()?;
    let copy = b.import_curve(&data)?;
    for t in [0.1, 0.3, 0.7, 0.9] {
        let [lo, hi] = copy.domain()?;
        let t = lo + t * (hi - lo);
        let evaluated = copy.evaluate(t)?;
        close(evaluated.point, data.eval(t).point, p);
        assert!((evaluated.d1 - data.eval(t).d1).length() <= p.linear());
    }
    let surface = b.import_surface(&plane(2.)?)?;
    let exported = surface.nurbs()?;
    assert_eq!(exported.knots(), plane(2.)?.knots());
    drop(b);
    close(
        surface.evaluate([0.25, 0.75])?.point,
        point(0.25, 0.75, 2.),
        p,
    );
    assert!(surface.evaluate([0.25, 0.75])?.duv.length() <= p.linear());
    close(
        copy.evaluate(copy.domain()?[0])?.point,
        circle.evaluate(circle.domain()?[0])?.point,
        p,
    );
    Ok(())
}
#[test]
fn intersections_keep_parameters_and_overlap_intervals() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let c = b.import_curve(&line(point(0., 0., 0.), point(1., 1., 0.))?)?;
    let other = b.import_curve(&line(point(0., 1., 0.), point(1., 0., 0.))?)?;
    let hits = c.intersect_curve(&other)?;
    assert_eq!(hits.len(), 1);
    close(
        c.evaluate(hits[0].start[0])?.point,
        other.evaluate(hits[0].start[1])?.point,
        p,
    );
    let offset = b.import_curve(&line(point(0., 2., 0.), point(1., 3., 0.))?)?;
    assert!(c.intersect_curve(&offset)?.is_empty());
    let overlaps = c.intersect_curve(&c)?;
    assert!(overlaps.iter().any(|s| s.native_kind == 2));
    let surface = b.import_surface(&plane(0.)?)?;
    let pierce = b.import_curve(&line(point(0.5, 0.5, -1.), point(0.5, 0.5, 1.))?)?;
    let hits = pierce.intersect_surface(&surface)?;
    assert_eq!(hits.len(), 1);
    let [t, u, v] = hits[0].start;
    close(
        pierce.evaluate(t)?.point,
        surface.evaluate([u, v])?.point,
        p,
    );
    let point_hits = surface.project_point(point(0.25, 0.75, 1.))?;
    assert!(!point_hits.is_empty());
    close(
        surface.evaluate(point_hits[0].start)?.point,
        point(0.25, 0.75, 0.),
        p,
    );
    let curve_hits = c.project_point(point(0.5, 0.5, 1.))?;
    assert!(!curve_hits.is_empty());
    close(
        c.evaluate(curve_hits[0].start[0])?.point,
        point(0.5, 0.5, 0.),
        p,
    );
    Ok(())
}
#[test]
fn surface_intersection_returns_paired_uv_curves() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let a = b.import_surface(&plane(0.)?)?;
    let q = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            point(0., 0.5, -1.),
            point(0., 0.5, 1.),
            point(1., 0.5, -1.),
            point(1., 0.5, 1.),
        ],
        None,
    )?;
    let other = b.import_surface(&q)?;
    let branches = a.intersect(&other)?;
    assert!(!branches.is_empty());
    for branch in branches {
        for fraction in [0., 0.3, 0.7, 1.] {
            let [lo, hi] = branch.curve.domain()?;
            let t = lo + (hi - lo) * fraction;
            let xyz = branch.curve.evaluate(t)?.point;
            let uv = branch.first_uv.evaluate(t)?.point;
            let st = branch.second_uv.evaluate(t)?.point;
            close(xyz, a.evaluate([uv.x, uv.y])?.point, p);
            close(xyz, other.evaluate([st.x, st.y])?.point, p);
        }
    }
    assert!(a.intersect(&b.import_surface(&plane(2.)?)?)?.is_empty());
    Ok(())
}
#[test]
fn copies_edits_and_surface_construction() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let original = b.import_curve(&line(point(0., 0., 0.), point(1., 0., 0.))?)?;
    let elevated = original.elevate_degree(3)?.insert_knots(&[0.4, 0.6])?;
    let before = original.nurbs()?;
    for t in [0., 0.25, 0.5, 1.] {
        close(original.evaluate(t)?.point, elevated.evaluate(t)?.point, p);
    }
    let reversed = elevated.reversed()?;
    close(
        reversed.evaluate(0.)?.point,
        original.evaluate(1.)?.point,
        p,
    );
    let trim = elevated.trim([0.2, 0.8])?;
    close(trim.evaluate(0.2)?.point, original.evaluate(0.2)?.point, p);
    let simplified = elevated.simplify()?;
    close(
        simplified.evaluate(0.5)?.point,
        original.evaluate(0.5)?.point,
        p,
    );
    assert_eq!(original.nurbs()?.points(), before.points());
    let s = b.extrude(&original, Vec3::new(0., 1., 0.))?;
    let [u0, u1, v0, v1] = s.domain()?;
    close(
        s.evaluate([(u0 + u1) / 2., (v0 + v1) / 2.])?.point,
        point(0.5, 0.5, 0.),
        p,
    );
    let (left, right) = s.split(Direction::U, (u0 + u1) / 2.)?;
    close(
        left.evaluate([(u0 + u1) / 2., v0])?.point,
        right.evaluate([(u0 + u1) / 2., v0])?.point,
        p,
    );
    let iso = s.iso_curve(Direction::U, (u0 + u1) / 2.)?;
    assert!(iso.domain()?[1] > iso.domain()?[0]);
    let raised = s
        .elevate_degree(Direction::U, 3)?
        .insert_knots(Direction::U, &[0.4])?;
    close(
        raised.evaluate([0.5, 0.5])?.point,
        s.evaluate([0.5, 0.5])?.point,
        p,
    );
    let flipped = s.reversed(Direction::V)?;
    close(
        flipped.evaluate([u0, v0])?.point,
        s.evaluate([u0, v1])?.point,
        p,
    );
    let cropped = s.trim([0.2, 0.8, 0.1, 0.9])?;
    close(
        cropped.evaluate([0.5, 0.5])?.point,
        s.evaluate([0.5, 0.5])?.point,
        p,
    );
    let second = b.import_curve(&line(point(0., 1., 0.), point(1., 1., 0.))?)?;
    close(
        b.ruled(&original, &second)?.evaluate([0.5, 0.5])?.point,
        point(0.5, 0.5, 0.),
        p,
    );
    let compatible = b.compatible_curves(&[original.clone(), second.clone()])?;
    assert_eq!(compatible.len(), 2);
    assert_eq!(b.loft(&compatible)?.nurbs()?.counts().0, 2);
    Ok(())
}
#[test]
fn projection_lift_offsets_and_factories() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let s = b.import_surface(&plane(0.)?)?;
    let uv = b.import_curve(&line(point(0.2, 0.3, 0.), point(0.8, 0.3, 0.))?)?;
    let lifted = uv.lift(&s)?;
    assert_eq!(lifted.len(), 1);
    close(lifted[0].evaluate(0.5)?.point, point(0.5, 0.3, 0.), p);
    let above = b.import_curve(&line(point(0.2, 0.3, 1.), point(0.8, 0.3, 1.))?)?;
    for projected in [above.drop_to(&s)?, above.project_to(&s)?] {
        assert!(!projected.uv.is_empty());
        assert_eq!(projected.uv.len(), projected.spatial.len());
        for curve in projected.spatial {
            let [a, z] = curve.domain()?;
            assert!(curve.evaluate((a + z) / 2.)?.point.z.abs() <= p.linear());
        }
    }
    let shifted = s.offset(0.5)?;
    close(shifted.evaluate([0.5, 0.5])?.point, point(0.5, 0.5, 0.5), p);
    let offset = uv.offset_in_plane(0.1, Vec3::new(0., 0., 1.))?;
    assert!(offset.nurbs()?.points().len() >= 2);
    assert_eq!(uv.lines_and_arcs()?.len(), 1);
    let rectangle = b.rectangle(point(0., 0., 0.), [2., 3.])?;
    assert_eq!(rectangle.len(), 4);
    assert_eq!(b.order_curves(&rectangle)?.len(), 4);
    assert_eq!(b.regular_polygon(point(0., 0., 0.), 5, 2.)?.len(), 5);
    assert!(
        b.ellipse(point(0., 0., 0.), [2., 1.])?
            .evaluate(0.)?
            .point
            .is_finite()
    );
    assert!(
        b.arc(point(0., 0., 0.), 1., [0., std::f64::consts::PI])?
            .nurbs()?
            .degree()
            >= 2
    );
    let input = [
        point(0., 0., 0.),
        point(1., 1., 0.),
        point(2., 1., 0.),
        point(3., 0., 0.),
    ];
    for curve in [
        b.control_curve(&input)?,
        b.interpolate_curve(&input, Parameterization::ChordLength)?,
        b.approximate_curve(&input)?,
    ] {
        assert!(curve.evaluate(curve.domain()?[0])?.point.is_finite());
    }
    Ok(())
}
#[test]
fn invalid_inputs_and_session_mismatch_return_errors() -> Result<()> {
    let p = Precision::DEFAULT;
    let a = Backend::new(p)?;
    let b = Backend::new(p)?;
    let ca = a.circle(point(0., 0., 0.), 1.)?;
    let cb = b.circle(point(0., 0., 0.), 1.)?;
    assert_eq!(
        ca.intersect_curve(&cb).err().map(|e| e.code()),
        Some(ErrorCode::InvalidArgument)
    );
    assert!(ca.evaluate(f64::NAN).is_err());
    assert!(ca.trim([1., 0.]).is_err());
    assert!(ca.insert_knots(&[0., 0., 0., 0.]).is_err());
    assert!(a.circle(point(0., 0., 0.), p.linear() / 2.).is_err());
    assert!(a.control_surface(&[], usize::MAX, usize::MAX).is_err());
    for _ in 0..32 {
        assert!(ca.evaluate(f64::INFINITY).is_err());
    }
    assert!(ca.evaluate(ca.domain()?[0])?.point.is_finite());
    Ok(())
}
proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]
    #[test]
    fn rational_surfaces_agree_in_random_frames(
        size in 0.1f64..20., angle in -3.0f64..3., offset in -10.0f64..10.,
        weight in 0.25f64..4., u in 0.0f64..1., v in 0.0f64..1.
    ) {
        let p=Precision::DEFAULT;
        let frame=Frame::from_axis(point(offset,-offset,offset),Vec3::new(angle.cos(),angle.sin(),1.),&p)?;
        let control=[Vec3::new(0.,0.,0.),Vec3::new(0.,size,0.),Vec3::new(size,0.,0.),Vec3::new(size,size,size/3.)];
        let source=BSurface::new(1,1,vec![0.,0.,1.,1.],vec![0.,0.,1.,1.],2,2,
            control.into_iter().map(|q|frame.to_world(q)).collect(),Some(vec![1.,weight,1.,weight]))?;
        let b=Backend::new(p)?;let s=b.import_surface(&source)?;let actual=s.evaluate([u,v])?;let expected=source.eval(u,v);
        prop_assert!(actual.point.distance(expected.point)<=p.linear());
        prop_assert!((actual.du-expected.du).length()<=p.linear());
        prop_assert!((actual.dv-expected.dv).length()<=p.linear());
        let exported=s.nurbs()?;prop_assert_eq!(exported.knots(),source.knots());
    }
}

#[test]
fn grids_sweep_revolution_and_helix() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let corners = [
        point(0., 0., 0.),
        point(1., 0., 0.),
        point(0., 1., 0.),
        point(1., 1., 0.),
    ];
    close(
        b.corner_surface(corners)?.evaluate([0.5, 0.5])?.point,
        point(0.5, 0.5, 0.),
        p,
    );
    let grid: Vec<_> = (0..4)
        .flat_map(|i| (0..4).map(move |j| point(f64::from(i) / 3., f64::from(j) / 3., 0.)))
        .collect();
    for s in [
        b.control_surface(&grid, 4, 4)?,
        b.interpolate_surface(&grid, 4, 4)?,
        b.approximate_surface(&grid)?,
    ] {
        let [u0, u1, v0, v1] = s.domain()?;
        let e = s.evaluate([(u0 + u1) / 2., (v0 + v1) / 2.])?;
        assert!(e.point.is_finite());
        assert!(e.point.z.abs() <= p.linear());
    }
    let profile = b.circle(point(0., 0., 0.), 0.25)?;
    let path = b.import_curve(&line(point(0., 0., 0.), point(0., 0., 2.))?)?;
    let swept = b.sweep(&profile, &path)?;
    let [u0, u1, v0, v1] = swept.domain()?;
    assert!(
        swept
            .evaluate([(u0 + u1) / 2., (v0 + v1) / 2.])?
            .point
            .is_finite()
    );
    let meridian = b.import_curve(&line(point(2., 0., 0.), point(2., 0., 1.))?)?;
    let revolved = b.revolve(
        &meridian,
        point(0., 0., 0.),
        Vec3::new(0., 0., 1.),
        std::f64::consts::TAU,
    )?;
    let [u0, u1, v0, v1] = revolved.domain()?;
    let pt = revolved.evaluate([(u0 + u1) / 2., (v0 + v1) / 2.])?.point;
    assert!((pt.x.hypot(pt.y) - 2.).abs() <= p.linear());
    let helix_backend = Backend::new(Precision::new(p.linear() * 1024., p.angular(), p.size())?)?;
    let helix = helix_backend.helix(point(0., 0., 0.), 1., [1., 1.], 1., true)?;
    let [lo, hi] = helix.domain()?;
    assert!(helix.evaluate((lo + hi) / 2.)?.point.is_finite());
    assert!(
        b.line_segment(point(0., 0., 0.), point(1., 0., 0.))?
            .evaluate(0.5)?
            .point
            .is_finite()
    );
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]
    #[test]
    fn near_degenerate_inputs_fail_without_poisoning_the_session(factor in 0.01f64..0.99, offset in -10.0f64..10.) {
        let p=Precision::DEFAULT;let b=Backend::new(p)?;
        prop_assert!(b.line_segment(point(offset,0.,0.),point(offset+p.linear()*factor,0.,0.)).is_err());
        let valid=b.import_surface(&plane(0.)?)?;
        prop_assert!(valid.evaluate([0.5,0.5])?.point.is_finite());
    }
}

#[test]
fn unequal_surface_dimensions_and_rational_curve_derivatives() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let points: Vec<_> = (0..3)
        .flat_map(|i| (0..4).map(move |j| point(f64::from(i), f64::from(j), f64::from(i * j) / 3.)))
        .collect();
    let source = BSurface::new(
        2,
        3,
        vec![0., 0., 0., 1., 1., 1.],
        vec![0., 0., 0., 0., 1., 1., 1., 1.],
        3,
        4,
        points,
        Some((0..12).map(|i| 1. + f64::from(i) / 12.).collect()),
    )?;
    let native = b.import_surface(&source)?;
    assert_eq!(native.nurbs()?.degrees(), source.degrees());
    assert_eq!(native.nurbs()?.counts(), source.counts());
    for uv in [[0., 0.], [0.2, 0.7], [1., 1.]] {
        let actual = native.evaluate(uv)?;
        let expected = source.eval(uv[0], uv[1]);
        close(actual.point, expected.point, p);
        assert!((actual.du - expected.du).length() <= p.linear());
        assert!((actual.dv - expected.dv).length() <= p.linear());
    }
    let rational = BCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![point(1., 0., 0.), point(1., 1., 0.), point(0., 1., 0.)],
        Some(vec![1., std::f64::consts::FRAC_1_SQRT_2, 1.]),
    )?;
    let curve = b.import_curve(&rational)?;
    let value = curve.evaluate(0.5)?;
    close(value.point, rational.eval(0.5).point, p);
    assert!(value.d2.is_finite());
    assert!(value.d2.length() > p.linear());
    let other = b.import_curve(&rational)?;
    assert!(
        curve
            .intersect_curve(&other)?
            .iter()
            .any(|s| s.native_kind == 2)
    );
    Ok(())
}

#[test]
fn lazy_analytic_geometry_can_be_edited_as_nurbs() -> Result<()> {
    let p = Precision::DEFAULT;
    let b = Backend::new(p)?;
    let native = b.line_segment(point(0., 0., 0.), point(2., 0., 0.))?;
    let raised = native.elevate_degree(3)?;
    let [lo, hi] = native.domain()?;
    close(
        raised.evaluate((lo + hi) / 2.)?.point,
        native.evaluate((lo + hi) / 2.)?.point,
        p,
    );
    close(
        native.reversed()?.evaluate(lo)?.point,
        native.evaluate(hi)?.point,
        p,
    );
    assert!(!native.lines_and_arcs()?.is_empty());
    Ok(())
}
