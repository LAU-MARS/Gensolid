//! Algebraic cubic span extraction and certified monotonic root classification.
//!
//! The derivative is quadratic: every one of its real roots partitions the cubic
//! into monotonic intervals. Certified signs at that complete list of extrema and
//! endpoints determine all distinct roots, including tangencies. Position refinement
//! occurs only after a monotonic interval is certified to contain one root.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Point3, Sign, Vec3};

use super::super::certified::{Scalar, Vector};
use super::{PointRoot, SpanRoot, endpoint_id, quadratic_roots};
use crate::{BCurve, Plane};

const DEGREE: usize = 3;
const MAX_ROOT_ITERATIONS: usize = 128;

/// Polynomial de Boor recursion, carried out in Bernstein form on this span.
/// Multiplication by a linear blend raises degree using the exact Bernstein
/// product weights (r-i)/r and i/r, retaining interval enclosures throughout.
fn extract(curve: &BCurve, span: usize, values: [Scalar; 4]) -> Vec<Scalar> {
    let knots = curve.knots();
    let lo = Scalar::new(knots[span]);
    let hi = Scalar::new(knots[span + 1]);
    let one = Scalar::new(1.0);
    let mut work: Vec<Vec<Scalar>> = values.into_iter().map(|value| vec![value]).collect();
    for r in 1..=DEGREE {
        for j in (r..=DEGREE).rev() {
            let left_knot = Scalar::new(knots[span - DEGREE + j]);
            let denominator = Scalar::new(knots[span + j - r + 1]) - left_knot;
            let alpha0 = (lo - left_knot) / denominator;
            let alpha1 = (hi - left_knot) / denominator;
            let mut next = Vec::with_capacity(r + 1);
            for i in 0..=r {
                let mut value = Scalar::new(0.0);
                if i < r {
                    value = value
                        + (Scalar::new((r - i) as f64) / Scalar::new(r as f64))
                            * ((one - alpha0) * work[j - 1][i] + alpha0 * work[j][i]);
                }
                if i != 0 {
                    value = value
                        + (Scalar::new(i as f64) / Scalar::new(r as f64))
                            * ((one - alpha1) * work[j - 1][i - 1] + alpha1 * work[j][i - 1]);
                }
                next.push(value);
            }
            work[j] = next;
        }
    }
    work.remove(DEGREE)
}

fn evaluate(coefficients: &[Scalar], t: Scalar) -> Scalar {
    let mut work = coefficients.to_vec();
    let one_minus_t = Scalar::new(1.0) - t;
    for remaining in (1..work.len()).rev() {
        for i in 0..remaining {
            work[i] = one_minus_t * work[i] + t * work[i + 1];
        }
    }
    work[0]
}

fn numerical_value(coefficients: &[Scalar], t: f64) -> f64 {
    let mut work: Vec<f64> = coefficients.iter().map(|value| value.value).collect();
    for remaining in (1..work.len()).rev() {
        for i in 0..remaining {
            work[i] = (1.0 - t) * work[i] + t * work[i + 1];
        }
    }
    work[0]
}

struct Geometry {
    reference: Point3,
    coordinates: Vec<Vec3>,
    weights: Vec<f64>,
    speed_bound: f64,
}

impl Geometry {
    fn evaluate(&self, t: f64) -> Point3 {
        let mut coordinates = self.coordinates.clone();
        let mut weights = self.weights.clone();
        for remaining in (1..coordinates.len()).rev() {
            for i in 0..remaining {
                coordinates[i] = coordinates[i] * (1.0 - t) + coordinates[i + 1] * t;
                weights[i] = weights[i] * (1.0 - t) + weights[i + 1] * t;
            }
        }
        self.reference + coordinates[0] / weights[0]
    }
}

fn geometry(curve: &BCurve, span: usize, weights: [Scalar; 4]) -> Result<Geometry> {
    let reference = curve.points()[span - DEGREE];
    let offsets: [Vector; 4] = std::array::from_fn(|i| {
        Vector::between(curve.points()[span - DEGREE + i], reference) * weights[i]
    });
    let xs = extract(curve, span, std::array::from_fn(|i| offsets[i].x));
    let ys = extract(curve, span, std::array::from_fn(|i| offsets[i].y));
    let zs = extract(curve, span, std::array::from_fn(|i| offsets[i].z));
    let extracted_weights = extract(curve, span, weights);
    let weights: Vec<f64> = extracted_weights
        .iter()
        .map(|weight| weight.value)
        .collect();
    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight <= 0.0)
    {
        return Err(Error::unsupported(
            "cubic spline homogeneous weights cannot be evaluated",
        ));
    }
    let coordinates: Vec<Vec3> = (0..=DEGREE)
        .map(|i| Vec3::new(xs[i].value, ys[i].value, zs[i].value))
        .collect();
    let control_points: Vec<Vec3> = coordinates
        .iter()
        .zip(&weights)
        .map(|(coordinate, weight)| *coordinate / *weight)
        .collect();
    let mut diameter: f64 = 0.0;
    for a in &control_points {
        for b in &control_points {
            diameter = diameter.max((*a - *b).length());
        }
    }
    let min_weight = weights.iter().copied().fold(f64::INFINITY, f64::min);
    let max_weight = weights.iter().copied().fold(0.0, f64::max);
    // The rational derivative's pairwise control-point expression is bounded by
    // 2*p*diameter*(max_weight/min_weight)^2 throughout [0,1].
    let ratio = max_weight / min_weight;
    let speed_bound = 2.0 * DEGREE as f64 * diameter * ratio * ratio;
    if !speed_bound.is_finite() {
        return Err(Error::unsupported(
            "cubic spline spatial refinement bound is not finite",
        ));
    }
    Ok(Geometry {
        reference,
        coordinates,
        weights,
        speed_bound,
    })
}

/// The branch count has already been certified from derivative roots and interval
/// signs. These floating comparisons refine its position only; they never choose
/// whether a root exists, split a branch, or merge two branches.
fn refine(
    coefficients: &[Scalar],
    lo: Scalar,
    hi: Scalar,
    geometry: &Geometry,
    precision: &Precision,
) -> Result<Scalar> {
    let (mut left, mut right) = (lo.value, hi.value);
    let mut left_value = numerical_value(coefficients, left);
    for _ in 0..MAX_ROOT_ITERATIONS {
        let midpoint = left + (right - left) * 0.5;
        let value = numerical_value(coefficients, midpoint);
        if !value.is_finite() {
            return Err(Error::unsupported(
                "cubic spline root refinement is not finite",
            ));
        }
        if (right - left) * geometry.speed_bound <= precision.solver()
            || midpoint <= left
            || midpoint >= right
        {
            return Ok(Scalar::new(midpoint));
        }
        if value.is_sign_negative() == left_value.is_sign_negative() {
            left = midpoint;
            left_value = value;
        } else {
            right = midpoint;
        }
    }
    Err(Error::unsupported(
        "cubic spline root position did not converge",
    ))
}

fn roots(
    coefficients: &[Scalar],
    geometry: &Geometry,
    precision: &Precision,
) -> Result<Option<Vec<SpanRoot>>> {
    let tolerance = precision.linear()
        * geometry
            .weights
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
    let signs: Vec<Sign> = coefficients
        .iter()
        .map(|value| value.sign(tolerance))
        .collect::<Result<_>>()?;
    if signs.iter().all(|sign| *sign == Sign::Zero) {
        return Ok(None);
    }
    if signs.iter().all(|sign| *sign == Sign::Positive)
        || signs.iter().all(|sign| *sign == Sign::Negative)
    {
        return Ok(Some(Vec::new()));
    }
    let derivative: [Scalar; 3] = std::array::from_fn(|i| {
        Scalar::new(DEGREE as f64) * (coefficients[i + 1] - coefficients[i])
    });
    let Some(critical_roots) = quadratic_roots(derivative, precision)? else {
        return Err(Error::unsupported(
            "cubic spline derivative degree is numerically unresolved",
        ));
    };
    let mut critical: Vec<Scalar> = critical_roots
        .into_iter()
        .filter_map(|root| match root {
            SpanRoot::Interior(t) => Some(t),
            SpanRoot::Start | SpanRoot::End => None,
        })
        .collect();
    if critical.len() == 2 {
        match (critical[1] - critical[0]).strict_sign()? {
            Sign::Negative => critical.swap(0, 1),
            Sign::Positive => {}
            Sign::Zero => return Err(Error::unsupported("cubic spline extrema are not separated")),
        }
    }
    let mut nodes = vec![(Scalar::new(0.0), SpanRoot::Start, signs[0])];
    for t in critical {
        let sign = evaluate(coefficients, t).sign(tolerance)?;
        nodes.push((t, SpanRoot::Interior(t), sign));
    }
    nodes.push((Scalar::new(1.0), SpanRoot::End, signs[DEGREE]));
    // A stationary zero has multiplicity at least two. If tolerance fallbacks
    // would require more multiplicity than this cubic possesses, the cluster
    // cannot be classified: do not turn a near-triple root into two tangencies.
    let minimum_multiplicity: usize = nodes
        .iter()
        .filter(|node| node.2 == Sign::Zero)
        .map(|node| {
            if matches!(node.1, SpanRoot::Interior(_)) {
                2
            } else {
                1
            }
        })
        .sum();
    if minimum_multiplicity > DEGREE {
        return Err(Error::unsupported(
            "cubic spline multiple-root cluster is numerically unresolved",
        ));
    }
    let mut roots = Vec::new();
    for (_, root, sign) in &nodes {
        if *sign == Sign::Zero {
            roots.push(*root);
        }
    }
    for pair in nodes.windows(2) {
        if matches!(
            (pair[0].2, pair[1].2),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        ) {
            let t = refine(coefficients, pair[0].0, pair[1].0, geometry, precision)?;
            roots.push(SpanRoot::Interior(t));
        }
    }
    Ok(Some(roots))
}

pub(super) fn span(
    curve: &BCurve,
    span: usize,
    plane: &Plane,
    precision: &Precision,
) -> Result<Option<Vec<PointRoot>>> {
    let active_weights = curve
        .weights()
        .map(|weights| &weights[span - DEGREE..=span]);
    let scale = active_weights.map_or(1.0, |weights| weights.iter().copied().fold(0.0, f64::max));
    let weights: [Scalar; 4] = std::array::from_fn(|i| {
        Scalar::new(active_weights.map_or(1.0, |weights| weights[i])) / Scalar::new(scale)
    });
    let normal = Vector::new(plane.frame().z());
    let heights = std::array::from_fn(|i| {
        Vector::between(curve.points()[span - DEGREE + i], plane.frame().origin()).dot(normal)
            * weights[i]
    });
    let coefficients = extract(curve, span, heights);
    let geometry = geometry(curve, span, weights)?;
    let Some(roots) = roots(&coefficients, &geometry, precision)? else {
        return Ok(None);
    };
    let mut points = Vec::new();
    for root in roots {
        let (t, endpoint) = match root {
            SpanRoot::Start => (0.0, endpoint_id(curve, span)),
            SpanRoot::End => (1.0, endpoint_id(curve, span + 1)),
            SpanRoot::Interior(t) => (t.value, None),
        };
        let point = geometry.evaluate(t);
        if !point.is_finite()
            || (point - plane.frame().origin())
                .dot(plane.frame().z())
                .abs()
                > precision.linear()
        {
            return Err(Error::unsupported(
                "cubic spline root evaluation exceeds session precision",
            ));
        }
        points.push(PointRoot { point, endpoint });
    }
    Ok(Some(points))
}

#[cfg(test)]
mod tests {
    use gensolid_math::Frame;
    use proptest::prelude::*;

    use super::super::{
        intersect,
        tests::{check_points, frame_strategy},
    };
    use super::*;
    use crate::{Curve, IntersectionCase, IntersectionSet, Surface, intersection_score};

    fn bernstein([constant, linear, quadratic, cubic]: [f64; 4]) -> [f64; 4] {
        [
            constant,
            constant + linear / 3.0,
            constant + 2.0 * linear / 3.0 + quadratic / 3.0,
            constant + linear + quadratic + cubic,
        ]
    }

    fn factored(roots: [f64; 3]) -> [f64; 4] {
        let [a, b, c] = roots;
        bernstein([-a * b * c, a * b + a * c + b * c, -a - b - c, 1.0])
    }

    fn bezier(frame: Frame, scale: f64, heights: [f64; 4], weights: [f64; 4]) -> BCurve {
        let points = (0..4)
            .map(|i| {
                frame.origin()
                    + frame.x() * (scale * i as f64 / 3.0)
                    + frame.z() * (scale * heights[i] / weights[i])
            })
            .collect();
        BCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            points,
            Some(weights.to_vec()),
        )
        .unwrap()
    }

    #[test]
    fn cubic_double_triple_and_endpoint_roots_are_distinct() {
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([0.25, 0.25, 0.75]), [1.0; 4]),
            Frame::WORLD,
            &[0.25, 0.75],
        );
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([0.5, 0.5, 0.5]), [1.0; 4]),
            Frame::WORLD,
            &[0.5],
        );
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([0.0, 0.5, 1.0]), [1.0; 4]),
            Frame::WORLD,
            &[0.0, 0.5, 1.0],
        );
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([0.0, 0.0, 0.0]), [1.0; 4]),
            Frame::WORLD,
            &[0.0],
        );
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([1.0, 1.0, 1.0]), [1.0; 4]),
            Frame::WORLD,
            &[1.0],
        );
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([0.0, 0.0, 0.75]), [1.0; 4]),
            Frame::WORLD,
            &[0.0, 0.75],
        );
        check_points(
            &bezier(Frame::WORLD, 1.0, factored([1.0, 1.0, 0.25]), [1.0; 4]),
            Frame::WORLD,
            &[0.25, 1.0],
        );
        check_points(
            &bezier(
                Frame::WORLD,
                1.0,
                bernstein([-0.125, 0.25, -0.5, 1.0]),
                [1.0; 4],
            ),
            Frame::WORLD,
            &[0.5],
        );
        let empty = bezier(Frame::WORLD, 1.0, [1.0; 4], [1.0; 4]);
        assert_eq!(
            intersect(&empty, &Plane::new(Frame::WORLD), &Precision::DEFAULT).unwrap(),
            IntersectionSet::Empty
        );
        let coplanar = bezier(Frame::WORLD, 1.0, [0.0; 4], [1.0, 2.0, 3.0, 4.0]);
        assert!(matches!(
            intersect(&coplanar, &Plane::new(Frame::WORLD), &Precision::DEFAULT).unwrap(),
            IntersectionSet::Curves(_)
        ));
        let positive_cubic_with_extremum = bezier(
            Frame::WORLD,
            1.0,
            bernstein([0.3, -0.7, 0.0, 1.0]),
            [1.0; 4],
        );
        assert_eq!(
            intersect(
                &positive_cubic_with_extremum,
                &Plane::new(Frame::WORLD),
                &Precision::DEFAULT
            )
            .unwrap(),
            IntersectionSet::Empty
        );
    }

    #[test]
    fn unclamped_span_endpoints_and_interior_roots_are_complete() {
        let knots = vec![-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0];
        let points = (0..4)
            .map(|i| {
                let (s, t, u) = (knots[i + 1], knots[i + 2], knots[i + 3]);
                let height = s * t * u - (s * t + s * u + t * u) * 0.5 + (s + t + u) / 6.0;
                Point3::new((s + t + u) / 3.0, 0.0, height)
            })
            .collect();
        let curve = BCurve::new(3, knots, points, None).unwrap();
        check_points(&curve, Frame::WORLD, &[0.0, 0.5, 1.0]);
    }

    #[test]
    fn cubic_continuous_and_full_multiplicity_endpoints_remain_distinct() {
        let plane = Plane::new(Frame::WORLD);
        let p = Precision::DEFAULT;
        let continuous = BCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0],
            (0..7)
                .map(|i| Point3::new(i as f64, 0.0, i as f64 - 3.0))
                .collect(),
            Some(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
        )
        .unwrap();
        check_points(&continuous, Frame::WORLD, &[0.5]);
        let discontinuous = BCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0],
            (0..8)
                .map(|i| {
                    Point3::new(
                        i as f64,
                        0.0,
                        if i < 4 {
                            i as f64 - 3.0
                        } else {
                            i as f64 - 4.0
                        },
                    )
                })
                .collect(),
            None,
        )
        .unwrap();
        assert!(
            matches!(intersect(&discontinuous, &plane, &p).unwrap(), IntersectionSet::Points(points) if points.len() == 2)
        );
        let partial = BCurve::new(
            3,
            continuous.knots().to_vec(),
            (0..7)
                .map(|i| Point3::new(i as f64, 0.0, if i <= 3 { 0.0 } else { i as f64 - 3.0 }))
                .collect(),
            None,
        )
        .unwrap();
        assert_eq!(
            intersect(&partial, &plane, &p).unwrap_err().code(),
            gensolid_core::ErrorCode::Unimplemented
        );
    }

    #[test]
    fn original_sixty_four_cubic_catalog_cases_find_both_branches() {
        let p = Precision::DEFAULT;
        let mut seed = 20261013_u64;
        let mut random = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1_u64 << 53) as f64
        };
        let mut cases = Vec::new();
        for _ in 0..64 {
            let origin = Point3::new(
                (2.0 * random() - 1.0) * 25.0,
                (2.0 * random() - 1.0) * 25.0,
                (2.0 * random() - 1.0) * 25.0,
            );
            let frame = Frame::from_axis(
                origin,
                Vec3::new(
                    2.0 * random() - 1.0,
                    2.0 * random() - 1.0,
                    2.0 * random() - 1.0,
                ),
                &p,
            )
            .unwrap();
            let r = 0.1 + 10.0 * random();
            let curve = BCurve::new(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                vec![
                    origin - frame.x() * r + frame.z() * r,
                    origin - frame.x() * (r / 3.0) - frame.z() * (5.0 * r / 3.0),
                    origin + frame.x() * (r / 3.0) - frame.z() * (5.0 * r / 3.0),
                    origin + frame.x() * r + frame.z() * r,
                ],
                None,
            )
            .unwrap();
            cases.push(IntersectionCase::CurveSurface {
                curve: Curve::BCurve(curve),
                surface: Surface::Plane(Plane::new(frame)),
                expected_branches: 2,
            });
        }
        let score = intersection_score(&cases, &p);
        assert_eq!(score.matched, 64, "{score:?}");
        assert_eq!(score.failed, 0);
        assert_eq!(score.unimplemented, 0);
    }

    proptest! {
        #[test]
        fn common_weight_scale_and_parameter_frame_preserve_cubic_geometry(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            weights in prop::array::uniform4(0.2f64..5.0),
            common_log_weight in -200.0f64..200.0,
            log_knot_scale in -6.0f64..3.0,
            offset_factor in -10.0f64..10.0,
        ) {
            let scale = 10.0f64.powf(log_scale);
            let base = bezier(frame, scale, factored([0.2, 0.5, 0.8]), weights);
            let knot_scale = 10.0f64.powf(log_knot_scale);
            let lo = offset_factor * knot_scale;
            let curve = BCurve::new(3, vec![lo, lo, lo, lo, lo + knot_scale, lo + knot_scale, lo + knot_scale, lo + knot_scale],
                base.points().to_vec(), Some(weights.iter().map(|weight| weight * 10.0f64.powf(common_log_weight)).collect())).unwrap();
            check_points(&curve, frame, &[lo + knot_scale * 0.2, lo + knot_scale * 0.5, lo + knot_scale * 0.8]);
        }

        #[test]
        fn arbitrary_positive_weights_find_three_independent_known_roots(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            weights in prop::array::uniform4(0.2f64..5.0),
            r1 in 0.05f64..0.25, r2 in 0.35f64..0.65, r3 in 0.75f64..0.95,
        ) {
            let roots = [r1, r2, r3];
            let curve = bezier(frame, 10.0f64.powf(log_scale), factored(roots), weights);
            check_points(&curve, frame, &roots);
        }

        #[test]
        fn nonuniform_simple_knots_preserve_a_global_genuine_cubic(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            weights in prop::array::uniform7(0.2f64..5.0),
            k1 in 0.1f64..0.3, k2 in 0.4f64..0.6, k3 in 0.7f64..0.9,
        ) {
            let scale = 10.0f64.powf(log_scale);
            let roots = [k1 * 0.5, (k1 + k2) * 0.5, (k3 + 1.0) * 0.5];
            let sum = roots.iter().sum::<f64>();
            let pairs = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2];
            let product = roots.iter().product::<f64>();
            let knots = vec![0.0, 0.0, 0.0, 0.0, k1, k2, k3, 1.0, 1.0, 1.0, 1.0];
            let points = (0..7).map(|i| {
                let (s, t, u) = (knots[i + 1], knots[i + 2], knots[i + 3]);
                // The symmetric polar form of the factored cubic reproduces that
                // polynomial on arbitrary knot spans, independently of extraction.
                let height = s * t * u - sum * (s * t + s * u + t * u) / 3.0 + pairs * (s + t + u) / 3.0 - product;
                frame.origin() + frame.x() * (scale * (s + t + u) / 3.0) + frame.z() * (scale * height / weights[i])
            }).collect();
            let curve = BCurve::new(3, knots, points, Some(weights.to_vec())).unwrap();
            check_points(&curve, frame, &roots);
        }

        #[test]
        fn repeated_nonuniform_knots_find_all_six_genuine_cubic_roots(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            log_knot_scale in -6.0f64..3.0,
            weights in prop::array::uniform7(0.2f64..5.0),
            fraction in 0.05f64..0.95,
        ) {
            let scale = 10.0f64.powf(log_scale);
            let end = 10.0f64.powf(log_knot_scale);
            let middle = end * fraction;
            let heights = factored([0.2, 0.5, 0.8]);
            let points = (0..7).map(|i| {
                let height = if i <= 3 { heights[i] } else { -heights[i - 3] };
                frame.origin() + frame.x() * (scale * i as f64 / 3.0) + frame.z() * (scale * height / weights[i])
            }).collect();
            let curve = BCurve::new(3, vec![0.0, 0.0, 0.0, 0.0, middle, middle, middle, end, end, end, end], points, Some(weights.to_vec())).unwrap();
            check_points(&curve, frame, &[middle * 0.2, middle * 0.5, middle * 0.8,
                middle + (end - middle) * 0.2, middle + (end - middle) * 0.5, middle + (end - middle) * 0.8]);
        }

        #[test]
        fn near_span_boundaries_are_not_clamped_into_the_domain(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            weights in prop::array::uniform4(0.2f64..5.0),
            factor in 4.0f64..64.0,
        ) {
            let p = Precision::DEFAULT;
            let scale = 10.0f64.powf(log_scale);
            let gap = factor * p.linear() / scale;
            for sign in [-1.0, 1.0] {
                let root = sign * gap;
                let curve = bezier(frame, scale, factored([root, -1.0, -2.0]), weights);
                let hit = intersect(&curve, &Plane::new(frame), &p).unwrap();
                prop_assert_eq!(hit.branch_count(), if sign < 0.0 { 0 } else { 1 });
                if sign > 0.0 { check_points(&curve, frame, &[root]); }
            }
        }

        #[test]
        fn near_double_root_keeps_both_sides_and_the_other_simple_root(
            frame in frame_strategy(),
            log_scale in -5.0f64..2.0,
            weights in prop::array::uniform4(0.2f64..5.0),
            factor in 4.0f64..64.0,
        ) {
            let p = Precision::DEFAULT;
            let scale = 10.0f64.powf(log_scale);
            let gap = factor * p.linear() / scale;
            for side in [-1.0, 1.0] {
                let constant = 0.25 + side * gap;
                let heights = bernstein([-0.8 * constant, 0.8 + constant, -1.8, 1.0]);
                let curve = bezier(frame, scale, heights, weights);
                let hit = intersect(&curve, &Plane::new(frame), &p).unwrap();
                prop_assert_eq!(hit.branch_count(), if side < 0.0 { 3 } else { 1 });
                let roots = if side < 0.0 { vec![0.5 - gap.sqrt(), 0.5 + gap.sqrt(), 0.8] } else { vec![0.8] };
                check_points(&curve, frame, &roots);
            }
        }
    }
}
