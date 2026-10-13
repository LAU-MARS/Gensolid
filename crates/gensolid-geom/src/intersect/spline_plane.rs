//! Complete degree-one through degree-three rational spline sections by a plane.
//!
//! Each nonempty knot span is converted algebraically to a homogeneous Bernstein
//! polynomial. Positive weights keep its denominator positive, so all plane roots
//! are exactly the roots of its polynomial numerator. No search seeds or
//! geometric sampling determine the number of branches.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Point3, Sign};

use super::certified::{Scalar, Vector};
use super::{IntersectionSet, squared_length_tolerance};
use crate::{BCurve, Curve, Plane};

mod cubic;

struct PointRoot {
    point: Point3,
    /// Shared knot identity establishes continuity even if the two independently
    /// evaluated one-sided points differ slightly through roundoff.
    endpoint: Option<usize>,
}

#[derive(Clone, Copy)]
enum SpanRoot {
    Start,
    End,
    Interior(Scalar),
}

/// Map a knot endpoint to its shared identity only when its multiplicity guarantees
/// continuity. Full multiplicity can represent two distinct one-sided points.
fn endpoint_id(curve: &BCurve, index: usize) -> Option<usize> {
    let knots = curve.knots();
    let mut first = index;
    let mut last = index;
    while first != 0 && knots[first - 1] >= knots[index] {
        first -= 1;
    }
    while last + 1 < knots.len() && knots[last + 1] <= knots[index] {
        last += 1;
    }
    (last - first < curve.degree()).then_some(first)
}

fn add_point(
    points: &mut Vec<PointRoot>,
    candidate: PointRoot,
    precision: &Precision,
) -> Result<()> {
    if !candidate.point.is_finite() {
        return Err(Error::unsupported(
            "spline-plane intersection point is not finite",
        ));
    }
    for existing in points.iter() {
        if candidate.endpoint.is_some() && candidate.endpoint == existing.endpoint {
            return Ok(());
        }
        let delta = Vector::between(candidate.point, existing.point);
        if delta.zero(precision.linear())? {
            return Ok(());
        }
        let separation = delta.norm2() - Scalar::new(precision.linear()).sqr();
        if separation.sign(precision.linear() * precision.linear())? != Sign::Positive {
            return Err(Error::unsupported(
                "spline-plane point branches are below session resolution",
            ));
        }
    }
    points.push(candidate);
    Ok(())
}

fn interior(root: Scalar) -> Result<Option<SpanRoot>> {
    match root.strict_sign()? {
        Sign::Negative => return Ok(None),
        Sign::Zero => {
            return Err(Error::unsupported(
                "unexpected unresolved spline start root",
            ));
        }
        Sign::Positive => {}
    }
    match (Scalar::new(1.0) - root).strict_sign()? {
        Sign::Negative => Ok(None),
        Sign::Zero => Err(Error::unsupported("unexpected unresolved spline end root")),
        Sign::Positive => Ok(Some(SpanRoot::Interior(root))),
    }
}

fn quadratic_roots(bernstein: [Scalar; 3], precision: &Precision) -> Result<Option<Vec<SpanRoot>>> {
    let [c0, c1, c2] = bernstein;
    let signs = [
        c0.sign(precision.linear())?,
        c1.sign(precision.linear())?,
        c2.sign(precision.linear())?,
    ];
    if signs.iter().all(|sign| *sign == Sign::Zero) {
        return Ok(None);
    }
    // Bernstein basis functions are nonnegative and sum to one. This certifies
    // empty spans without relying on cancellation in power-basis coefficients.
    if signs.iter().all(|sign| *sign == Sign::Positive)
        || signs.iter().all(|sign| *sign == Sign::Negative)
    {
        return Ok(Some(Vec::new()));
    }
    let a = c0 - Scalar::new(2.0) * c1 + c2;
    let b = Scalar::new(2.0) * (c1 - c0);
    let mut found = Vec::new();
    if a.sign(precision.linear())? == Sign::Zero {
        if b.sign(precision.linear())? == Sign::Zero {
            return Err(Error::unsupported(
                "spline-plane polynomial degree is numerically unresolved",
            ));
        }
        if signs[0] == Sign::Zero {
            found.push(SpanRoot::Start);
        } else if signs[2] == Sign::Zero {
            found.push(SpanRoot::End);
        } else if let Some(root) = interior(-c0 / b)? {
            found.push(root);
        }
        return Ok(Some(found));
    }
    // Endpoint factorization avoids cancellation at shared knots. The endpoint
    // residual uses Scalar's explicit linear fallback only when its interval
    // cannot decide; inclusion is never inferred from a rounded parameter estimate.
    if signs[0] == Sign::Zero {
        found.push(SpanRoot::Start);
        if signs[2] == Sign::Zero {
            found.push(SpanRoot::End);
        } else if b.sign(precision.linear())? != Sign::Zero
            && let Some(root) = interior(-b / a)?
        {
            found.push(root);
        }
        return Ok(Some(found));
    }
    if signs[2] == Sign::Zero {
        found.push(SpanRoot::End);
        if (a - c0).sign(precision.linear())? != Sign::Zero
            && let Some(root) = interior(c0 / a)?
        {
            found.push(root);
        }
        return Ok(Some(found));
    }
    let discriminant = b.sqr() - Scalar::new(4.0) * a * c0;
    let scale = b.value.abs() + (a.value * c0.value).abs().sqrt();
    match discriminant.sign(squared_length_tolerance(scale, precision))? {
        Sign::Negative => {}
        Sign::Zero => {
            if let Some(root) = interior(-b / (Scalar::new(2.0) * a))? {
                found.push(root);
            }
        }
        Sign::Positive => {
            let sqrt = discriminant.sqrt();
            let signed_sqrt = if b.sign(precision.linear())? == Sign::Negative {
                -sqrt
            } else {
                sqrt
            };
            let q = -Scalar::new(0.5) * (b + signed_sqrt);
            for root in [q / a, c0 / q] {
                if let Some(root) = interior(root)? {
                    found.push(root);
                }
            }
        }
    }
    Ok(Some(found))
}

fn linear_span(
    curve: &BCurve,
    span: usize,
    plane: &Plane,
    precision: &Precision,
) -> Result<Option<Vec<PointRoot>>> {
    let (a, b) = (curve.points()[span - 1], curve.points()[span]);
    let normal = Vector::new(plane.frame().z());
    let da = Vector::between(a, plane.frame().origin()).dot(normal);
    let db = Vector::between(b, plane.frame().origin()).dot(normal);
    let (point, endpoint) = match (da.sign(precision.linear())?, db.sign(precision.linear())?) {
        (Sign::Zero, Sign::Zero) => return Ok(None),
        (Sign::Zero, _) => (a, endpoint_id(curve, span)),
        (_, Sign::Zero) => (b, endpoint_id(curve, span + 1)),
        (Sign::Negative, Sign::Positive) | (Sign::Positive, Sign::Negative) => {
            (a + (b - a) * (da / (da - db)).value, None)
        }
        _ => return Ok(Some(Vec::new())),
    };
    Ok(Some(vec![PointRoot { point, endpoint }]))
}

fn quadratic_span(
    curve: &BCurve,
    span: usize,
    plane: &Plane,
    precision: &Precision,
) -> Result<Option<Vec<PointRoot>>> {
    let knots = curve.knots();
    let width = Scalar::new(knots[span + 1]) - Scalar::new(knots[span]);
    let alpha = width / (Scalar::new(knots[span + 1]) - Scalar::new(knots[span - 1]));
    let beta = width / (Scalar::new(knots[span + 2]) - Scalar::new(knots[span]));
    let one = Scalar::new(1.0);
    let active_weights = curve.weights().map(|weights| &weights[span - 2..=span]);
    let weight_scale =
        active_weights.map_or(1.0, |weights| weights.iter().copied().fold(0.0, f64::max));
    let weights: [Scalar; 3] = std::array::from_fn(|i| {
        Scalar::new(active_weights.map_or(1.0, |weights| weights[i])) / Scalar::new(weight_scale)
    });
    let normal = Vector::new(plane.frame().z());
    let heights: [Scalar; 3] = std::array::from_fn(|i| {
        Vector::between(curve.points()[span - 2 + i], plane.frame().origin()).dot(normal)
            * weights[i]
    });
    let bernstein = [
        alpha * heights[0] + (one - alpha) * heights[1],
        heights[1],
        (one - beta) * heights[1] + beta * heights[2],
    ];
    let Some(roots) = quadratic_roots(bernstein, precision)? else {
        return Ok(None);
    };
    let w0 = (alpha * weights[0] + (one - alpha) * weights[1]).value;
    let w1 = weights[1].value;
    let w2 = ((one - beta) * weights[1] + beta * weights[2]).value;
    let reference = curve.points()[span - 1];
    let q0 = (curve.points()[span - 2] - reference) * (alpha * weights[0]).value;
    let q2 = (curve.points()[span] - reference) * (beta * weights[2]).value;
    let mut points = Vec::new();
    for root in roots {
        let (t, endpoint) = match root {
            SpanRoot::Start => (0.0, endpoint_id(curve, span)),
            SpanRoot::End => (1.0, endpoint_id(curve, span + 1)),
            SpanRoot::Interior(t) => (t.value, None),
        };
        let s = 1.0 - t;
        let denominator = w0 * s * s + 2.0 * w1 * s * t + w2 * t * t;
        let point = reference + (q0 * (s * s) + q2 * (t * t)) / denominator;
        let residual = (point - plane.frame().origin())
            .dot(plane.frame().z())
            .abs();
        if !point.is_finite() || residual > precision.linear() {
            return Err(Error::unsupported(
                "spline-plane root evaluation exceeds session precision",
            ));
        }
        points.push(PointRoot { point, endpoint });
    }
    Ok(Some(points))
}

pub(super) fn intersect(
    curve: &BCurve,
    plane: &Plane,
    precision: &Precision,
) -> Result<IntersectionSet> {
    if curve.degree() > 3 {
        return Err(Error::unimplemented(
            "spline-plane intersection requires degree one, two, or three",
        ));
    }
    let mut points = Vec::new();
    let mut span_count = 0;
    let mut overlaps = 0;
    for span in curve.degree()..curve.points().len() {
        if curve.knots()[span] >= curve.knots()[span + 1] {
            continue;
        }
        span_count += 1;
        let roots = if curve.degree() == 1 {
            linear_span(curve, span, plane, precision)?
        } else if curve.degree() == 2 {
            quadratic_span(curve, span, plane, precision)?
        } else {
            cubic::span(curve, span, plane, precision)?
        };
        if let Some(roots) = roots {
            for root in roots {
                add_point(&mut points, root, precision)?;
            }
        } else {
            overlaps += 1;
        }
    }
    if overlaps != 0 {
        let range = curve.range();
        let disconnected = curve.knots().windows(curve.degree() + 1).any(|window| {
            window[0] >= window[curve.degree()] && window[0] > range.lo && window[0] < range.hi
        });
        if overlaps == span_count && !disconnected {
            return Ok(IntersectionSet::Curves(vec![Curve::BCurve(curve.clone())]));
        }
        return Err(Error::unimplemented(
            "spline-plane intersection with partial overlap or disconnected curves is not implemented",
        ));
    }
    Ok(if points.is_empty() {
        IntersectionSet::Empty
    } else {
        IntersectionSet::Points(points.into_iter().map(|root| root.point).collect())
    })
}

#[cfg(test)]
mod tests {
    use gensolid_math::{Frame, Vec3};
    use proptest::prelude::*;

    use super::*;
    use crate::GeomStore;

    pub(super) fn frame_strategy() -> impl Strategy<Value = Frame> {
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

    pub(super) fn check_points(curve: &BCurve, frame: Frame, expected_parameters: &[f64]) {
        let p = Precision::DEFAULT;
        let hit = intersect(curve, &Plane::new(frame), &p).unwrap();
        let IntersectionSet::Points(points) = hit else {
            panic!("expected isolated points")
        };
        assert_eq!(points.len(), expected_parameters.len());
        let mut store = GeomStore::new();
        let id = store.add_curve(Curve::BCurve(curve.clone())).unwrap();
        for point in &points {
            assert!((*point - frame.origin()).dot(frame.z()).abs() <= p.linear());
            assert!(store.project_to_curve(id, *point, &p).unwrap().distance <= p.linear());
        }
        for parameter in expected_parameters {
            let expected = curve.eval(*parameter).point;
            assert!(
                points
                    .iter()
                    .any(|point| point.distance(expected) <= p.linear()),
                "parameter={parameter}, points={points:?}, expected={expected:?}"
            );
        }
    }

    #[test]
    fn quadratic_tangency_and_linear_degree_reduction_are_distinct() {
        let p = Precision::DEFAULT;
        let make = |heights: [f64; 3]| {
            BCurve::new(
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                heights
                    .into_iter()
                    .enumerate()
                    .map(|(i, h)| Point3::new(i as f64, 0.0, h))
                    .collect(),
                None,
            )
            .unwrap()
        };
        check_points(&make([0.25, -0.25, 0.25]), Frame::WORLD, &[0.5]);
        check_points(&make([-1.0, 0.0, 1.0]), Frame::WORLD, &[0.5]);
        check_points(&make([0.0, -1.0, 0.0]), Frame::WORLD, &[0.0, 1.0]);
        check_points(&make([0.0, 0.0, 1.0]), Frame::WORLD, &[0.0]);
        check_points(&make([1.0, 0.0, 0.0]), Frame::WORLD, &[1.0]);
        assert_eq!(
            intersect(&make([1.0, 1.0, 1.0]), &Plane::new(Frame::WORLD), &p).unwrap(),
            IntersectionSet::Empty
        );
        assert!(matches!(
            intersect(&make([0.0, 0.0, 0.0]), &Plane::new(Frame::WORLD), &p).unwrap(),
            IntersectionSet::Curves(_)
        ));
    }

    #[test]
    fn shared_knots_merge_only_continuous_endpoint_roots() {
        let p = Precision::DEFAULT;
        let plane = Plane::new(Frame::WORLD);
        let continuous = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
            vec![
                Point3::new(-2.0, 0.0, -1.0),
                Point3::new(-1.0, 0.0, -0.5),
                Point3::ORIGIN,
                Point3::new(1.0, 0.0, 0.5),
                Point3::new(2.0, 0.0, 1.0),
            ],
            Some(vec![1.0, 2.0, 3.0, 4.0, 5.0]),
        )
        .unwrap();
        check_points(&continuous, Frame::WORLD, &[0.5]);
        let discontinuous = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0],
            vec![
                Point3::new(-2.0, 0.0, -1.0),
                Point3::new(-1.0, 0.0, -0.5),
                Point3::ORIGIN,
                Point3::new(3.0, 0.0, 0.0),
                Point3::new(4.0, 0.0, 0.5),
                Point3::new(5.0, 0.0, 1.0),
            ],
            None,
        )
        .unwrap();
        assert!(
            matches!(intersect(&discontinuous, &plane, &p).unwrap(), IntersectionSet::Points(points) if points.len() == 2)
        );
        let coplanar = BCurve::new(
            2,
            continuous.knots().to_vec(),
            continuous
                .points()
                .iter()
                .map(|point| Point3::new(point.x, point.y, 0.0))
                .collect(),
            None,
        )
        .unwrap();
        assert!(matches!(
            intersect(&coplanar, &plane, &p).unwrap(),
            IntersectionSet::Curves(_)
        ));
        let partial = BCurve::new(
            2,
            continuous.knots().to_vec(),
            vec![
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::ORIGIN,
                Point3::new(1.0, 0.0, -1.0),
                Point3::new(2.0, 0.0, 1.0),
            ],
            None,
        )
        .unwrap();
        assert_eq!(
            intersect(&partial, &plane, &p).unwrap_err().code(),
            gensolid_core::ErrorCode::Unimplemented
        );
    }

    #[test]
    fn simple_knot_endpoint_roots_are_not_reported_twice() {
        let knots = vec![0.0, 0.0, 0.0, 0.25, 0.75, 1.0, 1.0, 1.0];
        let points = (0..5)
            .map(|i| {
                let (s, t) = (knots[i + 1], knots[i + 2]);
                Point3::new((s + t) * 0.5, 0.0, s * t - (s + t) * 0.5 + 0.25 * 0.75)
            })
            .collect();
        let curve = BCurve::new(2, knots, points, None).unwrap();
        check_points(&curve, Frame::WORLD, &[0.25, 0.75]);
    }

    #[test]
    fn distinct_point_branches_below_resolution_are_rejected() {
        let p = Precision::DEFAULT;
        let curve = BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(p.linear() * 0.125, 0.0, -3.0),
                Point3::new(p.linear() * 0.25, 0.0, 1.0),
            ],
            None,
        )
        .unwrap();
        assert_eq!(
            intersect(&curve, &Plane::new(Frame::WORLD), &p)
                .unwrap_err()
                .code(),
            gensolid_core::ErrorCode::Unsupported
        );
    }

    #[test]
    fn seeded_development_quadratic_catalog_finds_both_roots() {
        let p = Precision::DEFAULT;
        let mut seed = 20261013_u64;
        let mut random = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1_u64 << 53) as f64
        };
        for _ in 0..64 {
            let origin = Point3::new(
                (2.0 * random() - 1.0) * 25.0,
                (2.0 * random() - 1.0) * 25.0,
                (2.0 * random() - 1.0) * 25.0,
            );
            let axis = Vec3::new(
                2.0 * random() - 1.0,
                2.0 * random() - 1.0,
                2.0 * random() - 1.0,
            );
            let frame = Frame::from_axis(origin, axis, &p).unwrap();
            let radius = 0.1 + 10.0 * random();
            let curve = BCurve::new(
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                vec![
                    origin - frame.x() * radius + frame.z() * radius,
                    origin - frame.z() * (radius * 3.0),
                    origin + frame.x() * radius + frame.z() * radius,
                ],
                Some(vec![1.0, 1.0, 1.0]),
            )
            .unwrap();
            let left = 0.5 - std::f64::consts::FRAC_1_SQRT_2 * 0.5;
            check_points(&curve, frame, &[left, 1.0 - left]);
        }
    }

    proptest! {
        #[test]
        fn weighted_bezier_finds_both_known_roots(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            weights in prop::array::uniform3(0.2f64..5.0),
            log_common_weight in -200.0f64..200.0,
            r1 in 0.1f64..0.4,
            r2 in 0.6f64..0.9,
        ) {
            let scale = 10.0f64.powf(log_scale);
            let heights = [r1 * r2, r1 * r2 - (r1 + r2) * 0.5, (1.0 - r1) * (1.0 - r2)];
            let points = (0..3).map(|i| frame.origin() + frame.x() * (scale * i as f64) + frame.z() * (scale * heights[i] / weights[i])).collect();
            let curve = BCurve::new(2, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], points,
                Some(weights.iter().map(|weight| weight * 10.0f64.powf(log_common_weight)).collect())).unwrap();
            check_points(&curve, frame, &[r1, r2]);
        }

        #[test]
        fn nonuniform_simple_knots_preserve_global_quadratic_roots(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            k1 in 0.1f64..0.4,
            k2 in 0.6f64..0.9,
            weights in prop::array::uniform5(0.2f64..5.0),
        ) {
            let scale = 10.0f64.powf(log_scale);
            let knots = vec![0.0, 0.0, 0.0, k1, k2, 1.0, 1.0, 1.0];
            let (r1, r2) = (k1 * 0.5, (k2 + 1.0) * 0.5);
            // The polar form of (u-r1)(u-r2) is st-(r1+r2)(s+t)/2+r1r2.
            // Its values at each adjacent knot pair are the B-spline coefficients.
            let points = (0..5).map(|i| {
                let (s, t) = (knots[i + 1], knots[i + 2]);
                let h = s * t - (r1 + r2) * (s + t) * 0.5 + r1 * r2;
                frame.origin() + frame.x() * (scale * (s + t) * 0.5) + frame.z() * (scale * h / weights[i])
            }).collect();
            let curve = BCurve::new(2, knots, points, Some(weights.to_vec())).unwrap();
            check_points(&curve, frame, &[r1, r2]);
        }

        #[test]
        fn repeated_nonuniform_knots_find_four_independent_span_roots(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            log_knot_scale in -6.0f64..3.0,
            fraction in 0.05f64..0.95,
            weights in prop::array::uniform5(0.2f64..5.0),
        ) {
            let scale = 10.0f64.powf(log_scale);
            let end = 10.0f64.powf(log_knot_scale);
            let middle = end * fraction;
            let points = (0..5).map(|i| frame.origin() + frame.x() * (scale * i as f64)
                + frame.z() * (scale * (if i % 2 == 0 { 1.0 } else { -3.0 }) / weights[i])).collect();
            let curve = BCurve::new(2, vec![0.0, 0.0, 0.0, middle, middle, end, end, end], points, Some(weights.to_vec())).unwrap();
            let left = 0.5 - std::f64::consts::FRAC_1_SQRT_2 * 0.5;
            let right = 1.0 - left;
            check_points(&curve, frame, &[middle * left, middle * right, middle + (end - middle) * left, middle + (end - middle) * right]);
        }

        #[test]
        fn tangency_both_sides_keep_distinct_root_counts(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            gap_factor in 4.0f64..64.0,
            weights in prop::array::uniform3(0.2f64..5.0),
        ) {
            let p = Precision::DEFAULT;
            let scale = 10.0f64.powf(log_scale);
            for side in [-1.0, 1.0] {
                let gap = side * gap_factor * p.linear();
                let heights = [0.25 * scale + gap, -0.25 * scale + gap, 0.25 * scale + gap];
                let points = (0..3).map(|i| frame.origin() + frame.x() * (scale * i as f64) + frame.z() * (heights[i] / weights[i])).collect();
                let curve = BCurve::new(2, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], points, Some(weights.to_vec())).unwrap();
                let hit = intersect(&curve, &Plane::new(frame), &p).unwrap();
                prop_assert_eq!(hit.branch_count(), if side < 0.0 { 2 } else { 0 });
                if side < 0.0 {
                    let delta = (-gap / scale).sqrt();
                    check_points(&curve, frame, &[0.5 - delta, 0.5 + delta]);
                }
            }
        }
    }
}
