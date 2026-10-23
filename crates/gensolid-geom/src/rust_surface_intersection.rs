//! Algebraic candidates for a polynomial Bezier extrusion and an affine plane patch.
//!
//! Original Rust implementation of S(u,v)=P(u)+v D and its affine plane projection.
//! No SMLib or OCCT source is translated. These are paired XYZ/UV candidates;
//! the modelling layer retains responsibility for whole-branch topology proofs.
use crate::{BCurve, BSurface, backend::SurfaceIntersectionCandidate};
use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Interval, Point3, Sign, Vec3};

mod spline;
pub use spline::intersect_spline_extrusion_plane;

fn uniform_weights(s: &BSurface) -> bool {
    s.weights()
        .is_none_or(|w| w.windows(2).all(|pair| pair[0].total_cmp(&pair[1]).is_eq()))
}

fn dot(a: Vec3, b: Vec3) -> Interval {
    Interval::point(a.x) * Interval::point(b.x)
        + Interval::point(a.y) * Interval::point(b.y)
        + Interval::point(a.z) * Interval::point(b.z)
}

fn dot_delta(a: Point3, b: Point3, n: Vec3) -> Interval {
    (Interval::point(a.x) - Interval::point(b.x)) * Interval::point(n.x)
        + (Interval::point(a.y) - Interval::point(b.y)) * Interval::point(n.y)
        + (Interval::point(a.z) - Interval::point(b.z)) * Interval::point(n.z)
}

fn inverse(value: Interval) -> Result<Interval> {
    if !matches!(value.sign(), Some(Sign::Positive | Sign::Negative)) {
        return Err(Error::unimplemented("plane section is not transverse"));
    }
    Interval::new((1.0 / value.hi()).next_down(), (1.0 / value.lo()).next_up())
        .filter(|v| v.lo().is_finite() && v.hi().is_finite())
        .ok_or_else(|| Error::unimplemented("plane section denominator overflow"))
}

fn midpoint(value: Interval) -> f64 {
    value.lo() + (value.hi() - value.lo()) * 0.5
}

fn magnitude(value: Interval) -> f64 {
    value.lo().abs().max(value.hi().abs())
}

fn parameter_spacing(lo: f64, hi: f64) -> f64 {
    let scale = lo.abs().max(hi.abs());
    scale.next_up() - scale
}

struct AffinePatch {
    origin: Point3,
    u: [Interval; 3],
    v: [Interval; 3],
    normal: Vec3,
    inverse_area: Interval,
}

fn interval_delta(a: Point3, b: Point3) -> [Interval; 3] {
    [
        Interval::point(a.x) - Interval::point(b.x),
        Interval::point(a.y) - Interval::point(b.y),
        Interval::point(a.z) - Interval::point(b.z),
    ]
}

fn cross_dot(a: [Interval; 3], b: [Interval; 3], n: Vec3) -> Interval {
    (a[1] * b[2] - a[2] * b[1]) * Interval::point(n.x)
        + (a[2] * b[0] - a[0] * b[2]) * Interval::point(n.y)
        + (a[0] * b[1] - a[1] * b[0]) * Interval::point(n.z)
}

fn affine_patch(s: &BSurface, p: &Precision) -> Result<AffinePatch> {
    if s.degrees() != (1, 1) || s.counts() != (2, 2) || !uniform_weights(s) {
        return Err(Error::unimplemented(
            "intersection tool needs an affine plane patch",
        ));
    }
    clamped(s.knots().0, 1, s.ranges().0.lo, s.ranges().0.hi)?;
    clamped(s.knots().1, 1, s.ranges().1.lo, s.ranges().1.hi)?;
    let q = s.points();
    let origin = q[0];
    let u = q[2] - origin;
    let v = q[1] - origin;
    // This is a representation check, not a decision that a body is planar.
    // A failed exact identity is accepted only within the explicit session solver bound.
    if (q[3] - q[2] - v).length() > p.solver() {
        return Err(Error::unimplemented(
            "planar tool has non-affine parameterization",
        ));
    }
    let normal = u.cross(v);
    let area = dot(normal, normal);
    if (area - Interval::point(p.degenerate_length()).sqr().sqr()).sign() != Some(Sign::Positive) {
        return Err(Error::degenerate("affine plane patch has degenerate axes"));
    }
    let u = interval_delta(q[2], origin);
    let v = interval_delta(q[1], origin);
    Ok(AffinePatch {
        origin,
        u,
        v,
        normal,
        inverse_area: inverse(cross_dot(u, v, normal))?,
    })
}

fn clamped(knots: &[f64], degree: usize, lo: f64, hi: f64) -> Result<()> {
    if knots.len() != 2 * (degree + 1)
        || knots[..=degree].iter().any(|x| !x.total_cmp(&lo).is_eq())
        || knots[degree + 1..]
            .iter()
            .any(|x| !x.total_cmp(&hi).is_eq())
    {
        return Err(Error::unimplemented("section needs clamped Bezier domains"));
    }
    Ok(())
}

fn section(
    first: &BSurface,
    second: &BSurface,
    linear_u: bool,
    p: &Precision,
) -> Result<Vec<SurfaceIntersectionCandidate>> {
    let plane = affine_patch(second, p)?;
    let (nu, nv) = first.counts();
    let (pu, pv) = first.degrees();
    let (uk, vk) = first.knots();
    let (ur, vr) = first.ranges();
    let (degree, count, knots, domain, linear, linear_knots) = if linear_u {
        if pu != 1 || nu != 2 {
            return Err(Error::unimplemented("surface needs one linear direction"));
        }
        (pv, nv, vk, vr, ur, uk)
    } else {
        if pv != 1 || nv != 2 {
            return Err(Error::unimplemented("surface needs one linear direction"));
        }
        (pu, nu, uk, ur, vr, vk)
    };
    if count != degree + 1 || degree > 16 || !uniform_weights(first) {
        return Err(Error::unimplemented(
            "section needs a polynomial Bezier extrusion",
        ));
    }
    clamped(knots, degree, domain.lo, domain.hi)?;
    clamped(linear_knots, 1, linear.lo, linear.hi)?;
    let pair = |i| {
        if linear_u {
            (first.points()[i], first.points()[nv + i])
        } else {
            (first.points()[i * nv], first.points()[i * nv + 1])
        }
    };
    let (a, b) = pair(0);
    let direction = b - a;
    let mut extrusion_error = 0.0_f64;
    for i in 0..count {
        let (a, b) = pair(i);
        extrusion_error = extrusion_error.max((b - a - direction).length());
    }
    if extrusion_error > p.solver() {
        return Err(Error::unimplemented(
            "ruled surface is not a translation extrusion",
        ));
    }
    let slab = second
        .points()
        .iter()
        .map(|q| dot_delta(*q, plane.origin, plane.normal))
        .fold(Interval::point(0.0), Interval::hull);
    let distances = first
        .points()
        .iter()
        .map(|q| dot_delta(*q, plane.origin, plane.normal))
        .collect::<Vec<_>>();
    // Positive Bernstein weights put both complete patches in their control hulls.
    // Compare with the tool's full slab, including its representation roundoff.
    if distances
        .iter()
        .all(|d| (*d - Interval::point(slab.hi())).sign() == Some(Sign::Positive))
        || distances
            .iter()
            .all(|d| (*d - Interval::point(slab.lo())).sign() == Some(Sign::Negative))
    {
        return Ok(Vec::new());
    }
    let reciprocal = inverse(dot_delta(b, a, plane.normal))?;
    let (tool_u, tool_v) = second.ranges();
    let mut xyz = Vec::with_capacity(count);
    let mut first_uv = Vec::with_capacity(count);
    let mut second_uv = Vec::with_capacity(count);
    let mut deviation = extrusion_error;
    let inverse_profile_span = inverse(Interval::point(domain.hi) - Interval::point(domain.lo))?;
    let inverse_linear_span = inverse(Interval::point(linear.hi) - Interval::point(linear.lo))?;
    let inverse_degree = inverse(Interval::point(degree as f64))?;
    let mut derivative_bound = 0.0_f64;
    for i in 1..count {
        let (a0, a1) = pair(i - 1);
        let (b0, b1) = pair(i);
        derivative_bound = derivative_bound
            .max((b0 - a0).length() * degree as f64)
            .max((b1 - a1).length() * degree as f64);
    }
    let mut first_uv_error = 0.0_f64;
    // Interior evaluations also round when adding the parameter origin back.
    // Exact control-point fractions do not suffice for a highly translated domain.
    let first_quantization = parameter_spacing(domain.lo, domain.hi)
        * magnitude(inverse_profile_span)
        * derivative_bound
        + parameter_spacing(linear.lo, linear.hi)
            * magnitude(inverse_linear_span)
            * direction.length();
    let second_quantization = parameter_spacing(tool_u.lo, tool_u.hi)
        * magnitude(inverse(
            Interval::point(tool_u.hi) - Interval::point(tool_u.lo),
        )?)
        * (second.points()[2] - second.points()[0]).length()
        + parameter_spacing(tool_v.lo, tool_v.hi)
            * magnitude(inverse(
                Interval::point(tool_v.hi) - Interval::point(tool_v.lo),
            )?)
            * (second.points()[1] - second.points()[0]).length();
    for i in 0..count {
        let (a, _) = pair(i);
        let fraction = -dot_delta(a, plane.origin, plane.normal) * reciprocal;
        if fraction.sign() != Some(Sign::Positive)
            || (Interval::point(1.0) - fraction).sign() != Some(Sign::Positive)
        {
            return Err(Error::unimplemented(
                "section needs a complete interior branch; clipping is not implemented",
            ));
        }
        let t = midpoint(fraction);
        let point = a + direction * t;
        let profile_t = domain.lo + (domain.hi - domain.lo) * (i as f64 / degree as f64);
        let linear_t = linear.lo + (linear.hi - linear.lo) * t;
        let uv = if linear_u {
            Point3::new(linear_t, profile_t, 0.0)
        } else {
            Point3::new(profile_t, linear_t, 0.0)
        };
        if !(linear_t > linear.lo && linear_t < linear.hi) {
            return Err(Error::unimplemented(
                "rounded section UV leaves the target domain",
            ));
        }
        // The parameter identity must survive f64 storage. Large translated domains
        // can round an interior parameter to a different fraction even for small geometry.
        // Bernstein convexity bounds the error over the whole parameter curve.
        let profile_fraction =
            (Interval::point(profile_t) - Interval::point(domain.lo)) * inverse_profile_span;
        let ideal_profile = Interval::point(i as f64) * inverse_degree;
        let linear_fraction =
            (Interval::point(linear_t) - Interval::point(linear.lo)) * inverse_linear_span;
        let lift = interval_delta(point, a);
        let d = interval_delta(b, pair(0).0);
        let residual = [
            lift[0] - d[0] * linear_fraction,
            lift[1] - d[1] * linear_fraction,
            lift[2] - d[2] * linear_fraction,
        ];
        let lift_error = magnitude(residual[0])
            .hypot(magnitude(residual[1]))
            .hypot(magnitude(residual[2]));
        first_uv_error = first_uv_error.max(
            lift_error
                + magnitude(profile_fraction - ideal_profile) * derivative_bound
                + extrusion_error,
        );
        let delta = interval_delta(point, plane.origin);
        let u = cross_dot(delta, plane.v, plane.normal) * plane.inverse_area;
        let v = cross_dot(plane.u, delta, plane.normal) * plane.inverse_area;
        for coordinate in [u, v] {
            // Conservative control-hull containment: partial branches are refused.
            if coordinate.sign() != Some(Sign::Positive)
                || (Interval::point(1.0) - coordinate).sign() != Some(Sign::Positive)
            {
                return Err(Error::unimplemented(
                    "section is not wholly inside the bounded affine tool",
                ));
            }
        }
        let tool_uv = Point3::new(
            tool_u.lo + (tool_u.hi - tool_u.lo) * midpoint(u),
            tool_v.lo + (tool_v.hi - tool_v.lo) * midpoint(v),
            0.0,
        );
        if !(tool_uv.x > tool_u.lo
            && tool_uv.x < tool_u.hi
            && tool_uv.y > tool_v.lo
            && tool_uv.y < tool_v.hi)
        {
            return Err(Error::unimplemented(
                "rounded section UV leaves the tool domain",
            ));
        }
        deviation = deviation.max(second.eval(tool_uv.x, tool_uv.y).point.distance(point));
        xyz.push(point);
        first_uv.push(uv);
        second_uv.push(tool_uv);
    }
    deviation = (deviation + second_quantization).max(first_uv_error + first_quantization);
    if !deviation.is_finite() || deviation > p.linear() {
        return Err(Error::unimplemented(
            "section representation exceeds session precision",
        ));
    }
    Ok(vec![SurfaceIntersectionCandidate {
        curve: BCurve::new(degree, knots.to_vec(), xyz, None)?,
        first_uv: BCurve::new(degree, knots.to_vec(), first_uv, None)?,
        second_uv: BCurve::new(degree, knots.to_vec(), second_uv, None)?,
        reported_deviation: deviation,
    }])
}

/// Paired algebraic intersection candidates for a supported extrusion/affine-plane pair.
/// Both surface orders and either extrusion parameter direction are supported.
/// Unsupported, partial, tangent and overlapping branches return `Unimplemented`.
pub fn intersect_extrusion_plane(
    first: &BSurface,
    second: &BSurface,
    p: &Precision,
) -> Result<Vec<SurfaceIntersectionCandidate>> {
    for linear_u in [false, true] {
        if let Ok(result) = section(first, second, linear_u, p) {
            return Ok(result);
        }
    }
    for linear_u in [false, true] {
        if let Ok(mut result) = section(second, first, linear_u, p) {
            for branch in &mut result {
                std::mem::swap(&mut branch.first_uv, &mut branch.second_uv);
            }
            return Ok(result);
        }
    }
    Err(Error::unimplemented(
        "Rust intersection supports complete polynomial Bezier extrusion / affine plane branches only",
    ))
}
