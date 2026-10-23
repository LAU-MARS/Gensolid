//! Independent homogeneous certificate for complete rational extrusion cuts.
//! Candidate construction is reused; no provider diagnostic is accepted as proof.
use super::*;
use crate::boundary::BoundaryEdge;
use gensolid_geom::{BCurve, rust_surface_intersection::intersect_spline_extrusion_plane};

fn absolute(x: Interval) -> f64 {
    x.lo().abs().max(x.hi().abs())
}
fn reciprocal(x: Interval) -> Result<Interval> {
    if !matches!(x.sign(), Some(Sign::Positive | Sign::Negative)) {
        return Err(Error::unimplemented(
            "rational cut denominator is unresolved",
        ));
    }
    Interval::new((1.0 / x.hi()).next_down(), (1.0 / x.lo()).next_up())
        .filter(|q| q.lo().is_finite() && q.hi().is_finite())
        .ok_or_else(|| Error::unimplemented("rational cut denominator overflow"))
}
fn coordinates(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn vector(p: Point3) -> [Interval; 3] {
    coordinates(p).map(Interval::point)
}
fn weight(c: &BCurve, i: usize) -> f64 {
    c.weights().map_or(1.0, |w| w[i])
}
fn same(a: f64, b: f64) -> bool {
    a.total_cmp(&b).is_eq()
}
fn clamped(knots: &[f64], degree: usize, r: gensolid_geom::ParamRange) -> bool {
    knots.len() == 2 * (degree + 1)
        && knots[..=degree].iter().all(|&k| same(k, r.lo))
        && knots[degree + 1..].iter().all(|&k| same(k, r.hi))
}
fn binomial(n: usize, k: usize) -> f64 {
    // The carrier degree cap (16) makes the largest product degree 33.
    // All binomial values and the integer recurrence intermediates fit u64,
    // and each final value is exactly representable in f64 at this degree.
    let mut result = 1_u64;
    for i in 1..=k.min(n - k) {
        result = result * (n + 1 - i) as u64 / i as u64;
    }
    result as f64
}
fn product(a: &[Interval], b: &[Interval]) -> Result<Vec<Interval>> {
    let p = a.len() - 1;
    let q = b.len() - 1;
    let mut out = vec![Interval::point(0.0); p + q + 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            let factor = Interval::point(binomial(p, i))
                * Interval::point(binomial(q, j))
                * reciprocal(Interval::point(binomial(p + q, i + j)))?;
            out[i + j] = out[i + j] + x * y * factor;
        }
    }
    Ok(out)
}
fn max_abs(values: &[Interval]) -> f64 {
    values.iter().copied().map(absolute).fold(0.0, f64::max)
}
fn spacing(lo: f64, hi: f64) -> f64 {
    let scale = lo.abs().max(hi.abs());
    scale.next_up() - scale
}

pub(super) fn section(
    s: &BSurface,
    tool: &BSurface,
    linear_u: bool,
    plane: Plane,
    p: &Precision,
) -> Result<(BoundaryEdge, [f64; 2])> {
    let (nu, nv) = s.counts();
    let (pu, pv) = s.degrees();
    let (uk, vk) = s.knots();
    let (u, v) = s.ranges();
    let (degree, count, knots, domain, linear, linear_knots) = if linear_u {
        if pu != 1 || nu != 2 {
            return Err(Error::unimplemented("rational cut needs linear rulings"));
        }
        (pv, nv, vk, v, u, uk)
    } else {
        if pv != 1 || nv != 2 {
            return Err(Error::unimplemented("rational cut needs linear rulings"));
        }
        (pu, nu, uk, u, v, vk)
    };
    if degree > 16
        || count != degree + 1
        || !clamped(knots, degree, domain)
        || !clamped(linear_knots, 1, linear)
    {
        return Err(Error::unimplemented(
            "rational oblique cut needs a single clamped Bezier span",
        ));
    }
    crate::patch::certify_graph(s, p)?;
    let indices = |i| {
        if linear_u {
            (i, nv + i)
        } else {
            (i * nv, i * nv + 1)
        }
    };
    let normal = plane.frame().z();
    let origin = plane.frame().origin();
    let mut source_weights = Vec::with_capacity(count);
    let mut a_points = Vec::with_capacity(count);
    let mut b_points = Vec::with_capacity(count);
    let mut slope = f64::INFINITY;
    let mut orientation = None;
    let mut ruling_bound = [0.0_f64; 3];
    let mut world_spacing = [0.0_f64; 3];
    for i in 0..count {
        let (ia, ib) = indices(i);
        let wa = s.weights().map_or(1.0, |w| w[ia]);
        let wb = s.weights().map_or(1.0, |w| w[ib]);
        if !same(wa, wb) {
            return Err(Error::unimplemented("rational ruling weights must match"));
        }
        let a = s.points()[ia];
        let b = s.points()[ib];
        let da = dot_delta(a, origin, normal);
        let db = dot_delta(b, origin, normal);
        let margin = Interval::point(p.linear());
        let sign = if (da + margin).sign() == Some(Sign::Negative)
            && (db - margin).sign() == Some(Sign::Positive)
        {
            Sign::Positive
        } else if (da - margin).sign() == Some(Sign::Positive)
            && (db + margin).sign() == Some(Sign::Negative)
        {
            Sign::Negative
        } else {
            return Err(Error::unimplemented(
                "rational cut is not globally transverse",
            ));
        };
        if orientation.is_some_and(|old| old != sign) {
            return Err(Error::unimplemented("rational cut changes orientation"));
        }
        orientation = Some(sign);
        let delta = db - da;
        slope = slope.min(if sign == Sign::Positive {
            delta.lo()
        } else {
            -delta.hi()
        });
        for k in 0..3 {
            ruling_bound[k] = ruling_bound[k].max(absolute(vector(b)[k] - vector(a)[k]));
            world_spacing[k] = world_spacing[k].max(spacing(coordinates(a)[k], coordinates(b)[k]));
        }
        source_weights.push(Interval::point(wa));
        a_points.push(a);
        b_points.push(b);
    }
    // The optional provider's curves are deliberately not used here. Generate
    // owned candidates from the original carriers, then certify their equations.
    let candidates = intersect_spline_extrusion_plane(s, tool, p)?;
    let [candidate] = candidates.as_slice() else {
        return Err(Error::unimplemented(
            "rational cut needs exactly one complete span",
        ));
    };
    let curve = &candidate.curve;
    let uv = &candidate.first_uv;
    if curve.degree() != degree
        || curve.knots() != knots
        || curve.points().len() != count
        || uv.degree() != degree + 1
        || !clamped(uv.knots(), degree + 1, domain)
        || (0..count).any(|i| !same(weight(curve, i), source_weights[i].lo()))
    {
        return Err(Error::unimplemented(
            "rational candidate representation does not match its carrier",
        ));
    }
    let profile_coordinate = |q: Point3| if linear_u { q.y } else { q.x };
    let linear_coordinate = |q: Point3| if linear_u { q.x } else { q.y };
    if !same(profile_coordinate(uv.points()[0]), domain.lo)
        || !same(profile_coordinate(uv.points()[degree + 1]), domain.hi)
    {
        return Err(Error::unimplemented(
            "rational UV does not cover both profile boundaries",
        ));
    }
    let inv_profile = reciprocal(Interval::point(domain.hi) - Interval::point(domain.lo))?;
    let inv_linear = reciprocal(Interval::point(linear.hi) - Interval::point(linear.lo))?;
    let uv_weights: Vec<_> = (0..uv.points().len())
        .map(|i| Interval::point(weight(uv, i)))
        .collect();
    let profile_numerator: Vec<_> = uv
        .points()
        .iter()
        .enumerate()
        .map(|(i, &q)| {
            (Interval::point(profile_coordinate(q)) - Interval::point(domain.lo))
                * inv_profile
                * uv_weights[i]
        })
        .collect();
    // Control ordinates need not be ordered for the rational identity u(t)=t.
    // Certify the derivative numerator U'W-UW' over its full Bernstein hull.
    let derivative = |coefficients: &[Interval]| {
        coefficients
            .windows(2)
            .map(|c| (c[1] - c[0]) * Interval::point((degree + 1) as f64))
            .collect::<Vec<_>>()
    };
    let positive = product(&derivative(&profile_numerator), &uv_weights)?;
    let negative = product(&profile_numerator, &derivative(&uv_weights))?;
    if positive
        .iter()
        .zip(negative)
        .any(|(&a, b)| (a - b).sign() != Some(Sign::Positive))
    {
        return Err(Error::unimplemented(
            "rational profile UV derivative certificate is unresolved",
        ));
    }
    let mut linear_numerator = Vec::with_capacity(uv.points().len());
    for (i, &q) in uv.points().iter().enumerate() {
        let t = (Interval::point(linear_coordinate(q)) - Interval::point(linear.lo)) * inv_linear;
        if t.sign() != Some(Sign::Positive)
            || (Interval::point(1.0) - t).sign() != Some(Sign::Positive)
        {
            return Err(Error::unimplemented(
                "rational UV leaves the ruling interval",
            ));
        }
        linear_numerator.push(t * uv_weights[i]);
    }
    let min_source = source_weights
        .iter()
        .map(|w| w.lo())
        .fold(f64::INFINITY, f64::min);
    let max_source = source_weights.iter().map(|w| w.hi()).fold(0.0, f64::max);
    let min_uv = uv_weights
        .iter()
        .map(|w| w.lo())
        .fold(f64::INFINITY, f64::min);
    let inverse_uv = reciprocal(Interval::point(min_uv))?;
    let inverse_product = reciprocal(Interval::point(min_source) * Interval::point(min_uv))?;
    let ratio = Interval::point(max_source) * reciprocal(Interval::point(min_source))?;
    // Compare U/W with t directly: (U elevated) - t W. This uses the
    // candidate's stored weights, independently of its construction algorithm.
    let elevated = product(&profile_numerator, &[Interval::point(1.0); 2])?;
    let identity = product(&uv_weights, &[Interval::point(0.0), Interval::point(1.0)])?;
    let residual: Vec<_> = elevated.iter().zip(identity).map(|(&a, b)| a - b).collect();
    let parameter_error = Interval::point(max_abs(&residual)) * inverse_uv
        + Interval::point(spacing(domain.lo, domain.hi)) * inv_profile;
    let linear_quantization = Interval::point(spacing(linear.lo, linear.hi)) * inv_linear;
    let mut lift = [Interval::point(0.0); 3];
    let mut plane_error = 0.0_f64;
    for point in curve.points() {
        plane_error = plane_error.max(absolute(dot_delta(*point, origin, normal)));
    }
    for (k, lift_axis) in lift.iter_mut().enumerate() {
        let displacement: Vec<_> = (0..count)
            .map(|i| (vector(curve.points()[i])[k] - vector(a_points[i])[k]) * source_weights[i])
            .collect();
        let ruling: Vec<_> = (0..count)
            .map(|i| (vector(b_points[i])[k] - vector(a_points[i])[k]) * source_weights[i])
            .collect();
        let left = product(&displacement, &uv_weights)?;
        let right = product(&ruling, &linear_numerator)?;
        let residual: Vec<_> = left.iter().zip(right).map(|(&a, b)| a - b).collect();
        let lo = a_points
            .iter()
            .chain(&b_points)
            .map(|q| coordinates(*q)[k])
            .fold(f64::INFINITY, f64::min);
        let hi = a_points
            .iter()
            .chain(&b_points)
            .map(|q| coordinates(*q)[k])
            .fold(f64::NEG_INFINITY, f64::max);
        // Quotient differentiation and positive Bernstein bases bound the
        // normalized-profile derivative by 2 p (w_max/w_min) hull_diameter.
        let derivative = Interval::point((2 * degree) as f64)
            * ratio
            * (Interval::point(hi) - Interval::point(lo));
        for q in curve.points() {
            world_spacing[k] =
                world_spacing[k].max(spacing(coordinates(*q)[k], coordinates(*q)[k]));
        }
        *lift_axis = Interval::point(max_abs(&residual)) * inverse_product
            + parameter_error * derivative
            + linear_quantization * Interval::point(ruling_bound[k])
            + Interval::point(world_spacing[k]) * Interval::point(2.0);
    }
    let mut source_plane_error = Interval::point(plane_error) + Interval::point(p.solver());
    for k in 0..3 {
        source_plane_error =
            source_plane_error + lift[k] * Interval::point([normal.x, normal.y, normal.z][k].abs());
    }
    let root_error = source_plane_error * reciprocal(Interval::point(slope))?;
    let mut squared = Interval::point(0.0);
    for k in 0..3 {
        squared = squared + (lift[k] + root_error * Interval::point(ruling_bound[k])).sqr();
    }
    if (squared - Interval::point(p.default_local()).sqr()).sign() != Some(Sign::Negative) {
        return Err(Error::unimplemented(
            "rational section lacks a whole-branch spatial error certificate",
        ));
    }
    let fraction = |q: Point3| (linear_coordinate(q) - linear.lo) / (linear.hi - linear.lo);
    Ok((
        BoundaryEdge {
            curve: curve.clone(),
            uv: Some(uv.clone()),
        },
        [fraction(uv.points()[0]), fraction(uv.points()[degree + 1])],
    ))
}
