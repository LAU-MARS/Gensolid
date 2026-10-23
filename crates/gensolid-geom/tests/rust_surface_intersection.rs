//! Paired XYZ/UV contract for the pure Rust, algebraic surface intersection path.
#![allow(clippy::unwrap_used)]
use gensolid_core::{ErrorCode, Precision};
use gensolid_geom::{BSurface, rust_surface_intersection::intersect_extrusion_plane};
use gensolid_math::Point3;

fn target(linear_u: bool) -> BSurface {
    let mut points = vec![];
    if linear_u {
        for y in [0., 2.] {
            for (x, z) in [(0., 0.), (1.5, 0.4), (3., 0.)] {
                points.push(Point3::new(x, y, z));
            }
        }
        BSurface::new(
            1,
            2,
            vec![3., 3., 7., 7.],
            vec![2., 2., 2., 5., 5., 5.],
            2,
            3,
            points,
            None,
        )
        .unwrap()
    } else {
        for (x, z) in [(0., 0.), (1.5, 0.4), (3., 0.)] {
            for y in [0., 2.] {
                points.push(Point3::new(x, y, z));
            }
        }
        BSurface::new(
            2,
            1,
            vec![2., 2., 2., 5., 5., 5.],
            vec![3., 3., 7., 7.],
            3,
            2,
            points,
            None,
        )
        .unwrap()
    }
}
fn tool(offset: f64, slope: f64, bounded: bool) -> BSurface {
    let xs = if bounded { [-5., 8.] } else { [1., 2.] };
    let points = xs
        .into_iter()
        .flat_map(|x| [-5., 5.].map(|z| Point3::new(x, offset + slope * (x + z), z)))
        .collect();
    BSurface::new(
        1,
        1,
        vec![-2., -2., 4., 4.],
        vec![8., 8., 12., 12.],
        2,
        2,
        points,
        None,
    )
    .unwrap()
}

#[test]
fn paired_curves_match_both_domains_and_operand_orders() {
    let p = Precision::DEFAULT;
    for linear_u in [false, true] {
        for slope in [0., 0.2] {
            let a = target(linear_u);
            let b = tool(0.7, slope, true);
            for (a, b) in [(&a, &b), (&b, &a)] {
                let hits = intersect_extrusion_plane(a, b, &p).unwrap();
                assert_eq!(hits.len(), 1);
                let hit = &hits[0];
                let r = hit.curve.range();
                for i in 0..=64 {
                    let t = r.lo + (r.hi - r.lo) * (i as f64 / 64.);
                    let xyz = hit.curve.eval(t).point;
                    let auv = hit.first_uv.eval(t).point;
                    let buv = hit.second_uv.eval(t).point;
                    assert!(xyz.distance(a.eval(auv.x, auv.y).point) <= p.linear());
                    assert!(xyz.distance(b.eval(buv.x, buv.y).point) <= p.linear());
                    assert!(auv.z.abs() <= p.solver() && buv.z.abs() <= p.solver());
                }
            }
        }
    }
}

#[test]
fn disjoint_is_proved_but_partial_and_overlap_are_not_faked_as_empty() {
    let p = Precision::DEFAULT;
    let a = target(false);
    assert!(
        intersect_extrusion_plane(&a, &tool(4., 0., true), &p)
            .unwrap()
            .is_empty()
    );
    for b in [tool(0.7, 0., false), tool(-0.1, 0.2, true)] {
        assert_eq!(
            intersect_extrusion_plane(&a, &b, &p).unwrap_err().code(),
            ErrorCode::Unimplemented
        );
    }
    let b = tool(0.7, 0., true);
    assert_eq!(
        intersect_extrusion_plane(&b, &b, &p).unwrap_err().code(),
        ErrorCode::Unimplemented
    );
}

#[test]
fn unsupported_rational_and_non_affine_tool_are_explicit() {
    let p = Precision::DEFAULT;
    let a = target(false);
    let b = tool(0.7, 0., true);
    let rational = BSurface::new(
        2,
        1,
        a.knots().0.to_vec(),
        a.knots().1.to_vec(),
        3,
        2,
        a.points().to_vec(),
        Some(vec![1., 1., 0.5, 0.5, 1., 1.]),
    )
    .unwrap();
    assert_eq!(
        intersect_extrusion_plane(&rational, &b, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let mut points = b.points().to_vec();
    points[3].x += 0.5;
    let non_affine = BSurface::new(
        1,
        1,
        b.knots().0.to_vec(),
        b.knots().1.to_vec(),
        2,
        2,
        points,
        None,
    )
    .unwrap();
    assert_eq!(
        intersect_extrusion_plane(&a, &non_affine, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}

#[test]
fn unrepresentable_uv_in_large_translated_domains_is_refused() {
    let p = Precision::DEFAULT;
    let a = target(false);
    let b = tool(0.7, 0., true);
    let lo = 1_000_000_000_000_000.;
    for (u, v) in [
        (a.knots().0.to_vec(), vec![lo, lo, lo + 1., lo + 1.]),
        (
            vec![lo, lo, lo, lo + 1., lo + 1., lo + 1.],
            a.knots().1.to_vec(),
        ),
    ] {
        let shifted = BSurface::new(2, 1, u, v, 3, 2, a.points().to_vec(), None).unwrap();
        assert_eq!(
            intersect_extrusion_plane(&shifted, &b, &p)
                .unwrap_err()
                .code(),
            ErrorCode::Unimplemented
        );
    }
}

#[test]
fn exactly_representable_uv_controls_do_not_hide_interior_quantization() {
    let p = Precision::DEFAULT;
    let a = target(false);
    let lo = 1_000_000_000_000_000.;
    let shifted = BSurface::new(
        2,
        1,
        a.knots().0.to_vec(),
        vec![lo, lo, lo + 1., lo + 1.],
        3,
        2,
        a.points().to_vec(),
        None,
    )
    .unwrap();
    // Exact control fractions 1/4, 1/2, 3/4 still have unrepresentable interior values.
    let plane = BSurface::new(
        1,
        1,
        vec![0., 0., 1., 1.],
        vec![0., 0., 1., 1.],
        2,
        2,
        vec![
            Point3::new(-3., -0.5, -5.),
            Point3::new(-3., -0.5, 5.),
            Point3::new(6., 2.5, -5.),
            Point3::new(6., 2.5, 5.),
        ],
        None,
    )
    .unwrap();
    assert_eq!(
        intersect_extrusion_plane(&shifted, &plane, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
    let shifted_tool = BSurface::new(
        1,
        1,
        vec![lo, lo, lo + 1., lo + 1.],
        vec![lo, lo, lo + 1., lo + 1.],
        2,
        2,
        plane.points().to_vec(),
        None,
    )
    .unwrap();
    assert_eq!(
        intersect_extrusion_plane(&a, &shifted_tool, &p)
            .unwrap_err()
            .code(),
        ErrorCode::Unimplemented
    );
}
