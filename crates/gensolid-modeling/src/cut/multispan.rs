//! Complete multi-span rational sections, certified against the raw stored basis.
//!
//! Floating knot insertion constructs candidates only. Independent outward-rounded
//! homogeneous insertion below encloses the original input and the final stored
//! curves. Each output span has a unit parameter interval; its UV retains the
//! original source knot endpoints. Integer C0 breaks make both topology reversals
//! representation-stable without changing the carrier or relaxing sewing.
use super::*;
use crate::{boundary::BoundaryEdge, spline_bounds::homogeneous_spans as raw_spans};
use gensolid_geom::{BCurve, ParamRange};

fn refused(message: &str) -> Error {
    Error::unimplemented(message)
}
fn absolute(x: Interval) -> f64 {
    x.lo().abs().max(x.hi().abs())
}
fn reciprocal(x: Interval) -> Result<Interval> {
    if !matches!(x.sign(), Some(Sign::Positive | Sign::Negative)) {
        return Err(refused("multi-span denominator is unresolved"));
    }
    Interval::new((1.0 / x.hi()).next_down(), (1.0 / x.lo()).next_up())
        .filter(|q| q.lo().is_finite() && q.hi().is_finite())
        .ok_or_else(|| refused("multi-span denominator overflow"))
}
fn coordinates(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn weight(c: &BCurve, i: usize) -> f64 {
    c.weights().map_or(1.0, |w| w[i])
}
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}
fn spacing(lo: f64, hi: f64) -> f64 {
    let scale = lo.abs().max(hi.abs());
    scale.next_up() - scale
}
fn max_abs(values: &[Interval]) -> f64 {
    values.iter().copied().map(absolute).fold(0.0, f64::max)
}
fn binomial(n: usize, k: usize) -> f64 {
    // degree <= 16 implies product degree <= 49. The largest binomial is
    // 63_205_303_218_876; all recurrence intermediates are below 2^53 and u64.
    // Casting the final integer is exact. Products of binomials use Interval.
    let mut result = 1_u64;
    for i in 1..=k.min(n - k) {
        result = result * (n + 1 - i) as u64 / i as u64;
    }
    result as f64
}
fn product(a: &[Interval], b: &[Interval]) -> Result<Vec<Interval>> {
    let (p, q) = (a.len() - 1, b.len() - 1);
    if p + q > 49 {
        return Err(refused("multi-span Bernstein product degree exceeds 49"));
    }
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
fn difference(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    a.iter().zip(b).map(|(&a, &b)| a - b).collect()
}
fn derivative(a: &[Interval]) -> Vec<Interval> {
    a.windows(2)
        .map(|w| (w[1] - w[0]) * Interval::point((a.len() - 1) as f64))
        .collect()
}
fn positive_min(values: &[Interval]) -> Result<f64> {
    if values
        .iter()
        .any(|w| !w.lo().is_finite() || !w.hi().is_finite() || w.sign() != Some(Sign::Positive))
    {
        return Err(refused("multi-span homogeneous weights are unresolved"));
    }
    Ok(values.iter().map(|w| w.lo()).fold(f64::INFINITY, f64::min))
}
fn breaks(curve: &BCurve) -> Result<Vec<f64>> {
    let p = curve.degree();
    let r = curve.range();
    if p > 16
        || curve.knots()[..=p].iter().any(|&k| !same(k, r.lo))
        || curve.knots()[curve.points().len()..]
            .iter()
            .any(|&k| !same(k, r.hi))
    {
        return Err(refused("multi-span cut needs clamped degree <= 16"));
    }
    let mut out = vec![r.lo];
    for &k in &curve.knots()[p + 1..curve.points().len()] {
        if !same(k, *out.last().ok_or_else(|| refused("empty knot range"))?) {
            out.push(k);
        }
    }
    out.push(r.hi);
    if out.len() > 257
        || out
            .windows(2)
            .any(|w| (Interval::point(w[1]) - Interval::point(w[0])).sign() != Some(Sign::Positive))
        || out[1..out.len() - 1]
            .iter()
            .any(|&k| curve.knots().iter().filter(|&&v| same(k, v)).count() > p)
    {
        return Err(refused("multi-span cut needs at most 256 continuous spans"));
    }
    Ok(out)
}
fn c0_knots(degree: usize, count: usize) -> Vec<f64> {
    let mut knots = vec![0.0; degree + 1];
    for i in 1..count {
        knots.extend(std::iter::repeat_n(i as f64, degree));
    }
    knots.extend(std::iter::repeat_n(count as f64, degree + 1));
    knots
}
fn insertions(curve: &BCurve, source_breaks: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    for &k in &source_breaks[1..source_breaks.len() - 1] {
        let n = curve.knots().iter().filter(|&&v| same(k, v)).count();
        out.extend(std::iter::repeat_n(k, curve.degree() - n));
    }
    out
}

fn construct(
    a: &BCurve,
    b: &BCurve,
    source_breaks: &[f64],
    linear: ParamRange,
    linear_u: bool,
    plane: Plane,
) -> Result<(BCurve, BCurve)> {
    let p = a.degree();
    let q = p + 1;
    let normal = plane.frame().z();
    let origin = plane.frame().origin();
    let direction = b.points()[0] - a.points()[0];
    let denominator = direction.dot(normal);
    let projected = a
        .points()
        .iter()
        .map(|&point| {
            let fraction = (origin - point).dot(normal) / denominator;
            point + direction * fraction
        })
        .collect();
    let raw = BCurve::new(
        p,
        a.knots().to_vec(),
        projected,
        a.weights().map(<[f64]>::to_vec),
    )?;
    let additions = insertions(a, source_breaks);
    let refined = a.insert_knots(&additions)?;
    let projected = raw.insert_knots(&additions)?;
    let count = source_breaks.len() - 1;
    let curve = BCurve::new(
        p,
        c0_knots(p, count),
        projected.points().to_vec(),
        projected.weights().map(<[f64]>::to_vec),
    )?;
    let scale = (0..refined.points().len())
        .map(|i| weight(&refined, i))
        .fold(0.0, f64::max);
    let weights: Vec<_> = (0..refined.points().len())
        .map(|i| weight(&refined, i) / scale)
        .collect();
    let linear_coordinates: Vec<_> = refined
        .points()
        .iter()
        .map(|&point| {
            linear.lo + (linear.hi - linear.lo) * (origin - point).dot(normal) / denominator
        })
        .collect();
    let mut uv_points = Vec::new();
    let mut uv_weights = Vec::new();
    for span in 0..count {
        for i in usize::from(span > 0)..=q {
            let alpha = i as f64 / q as f64;
            let left = if i <= p {
                weights[span * p + i] * (1.0 - alpha)
            } else {
                0.0
            };
            let right = if i > 0 {
                weights[span * p + i - 1] * alpha
            } else {
                0.0
            };
            let w = left + right;
            let profile_t = if i == 0 {
                source_breaks[span]
            } else if i == q {
                source_breaks[span + 1]
            } else {
                source_breaks[span] + (source_breaks[span + 1] - source_breaks[span]) * right / w
            };
            let linear_t = if i == 0 {
                linear_coordinates[span * p]
            } else if i == q {
                linear_coordinates[(span + 1) * p]
            } else {
                (left * linear_coordinates[span * p + i]
                    + right * linear_coordinates[span * p + i - 1])
                    / w
            };
            uv_points.push(if linear_u {
                Point3::new(linear_t, profile_t, 0.0)
            } else {
                Point3::new(profile_t, linear_t, 0.0)
            });
            uv_weights.push(w);
        }
    }
    let uv = BCurve::new(q, c0_knots(q, count), uv_points, Some(uv_weights))?;
    // All shared controls/weights come from one global refinement and are stored
    // only once. Both curves have exactly the same integer span boundaries.
    for c in [&curve, &uv] {
        let reversed = c.reversed()?;
        let restored = reversed.reversed()?;
        if c.knots() != restored.knots()
            || c.points() != restored.points()
            || c.weights() != restored.weights()
            || reversed.range() != c.range()
        {
            return Err(refused("multi-span section reversal is not stable"));
        }
    }
    Ok((curve, uv))
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
    let (degree, count, knots, linear, linear_knots) = if linear_u {
        if pu != 1 || nu != 2 {
            return Err(refused("multi-span cut needs linear rulings"));
        }
        (pv, nv, vk, u, uk)
    } else {
        if pv != 1 || nv != 2 {
            return Err(refused("multi-span cut needs linear rulings"));
        }
        (pu, nu, uk, v, vk)
    };
    if linear_knots.len() != 4
        || !same(linear_knots[0], linear.lo)
        || !same(linear_knots[1], linear.lo)
        || !same(linear_knots[2], linear.hi)
        || !same(linear_knots[3], linear.hi)
        || tool.degrees() != (1, 1)
        || tool.counts() != (2, 2)
        || tool
            .weights()
            .is_some_and(|ws| ws.iter().any(|&w| !same(w, ws[0])))
    {
        return Err(refused(
            "multi-span cut needs clamped rulings and an affine plane tool",
        ));
    }
    let mut parallelogram_squared = Interval::point(0.0);
    for k in 0..3 {
        let coordinate = |i| Interval::point(coordinates(tool.points()[i])[k]);
        let residual = coordinate(0) - coordinate(1) + coordinate(3) - coordinate(2);
        parallelogram_squared = parallelogram_squared + residual.sqr();
    }
    if (parallelogram_squared - Interval::point(p.solver()).sqr()).sign() != Some(Sign::Negative) {
        return Err(refused("multi-span tool is not affine"));
    }
    let normal = plane.frame().z();
    let origin = plane.frame().origin();
    let mut a_points = Vec::with_capacity(count);
    let mut b_points = Vec::with_capacity(count);
    let mut source_weights = Vec::with_capacity(count);
    let mut orientation = None;
    let mut slope = f64::INFINITY;
    let mut ruling_bound = [0.0_f64; 3];
    let mut world_spacing = [0.0_f64; 3];
    for i in 0..count {
        let (ia, ib) = if linear_u {
            (i, nv + i)
        } else {
            (i * nv, i * nv + 1)
        };
        let (a, b) = (s.points()[ia], s.points()[ib]);
        let wa = s.weights().map_or(1.0, |w| w[ia]);
        if !same(wa, s.weights().map_or(1.0, |w| w[ib])) {
            return Err(refused("multi-span ruling weights must match"));
        }
        let (da, db) = (dot_delta(a, origin, normal), dot_delta(b, origin, normal));
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
            return Err(refused("multi-span section is not globally transverse"));
        };
        if orientation.is_some_and(|old| old != sign) {
            return Err(refused("multi-span ruling orientation changes"));
        }
        orientation = Some(sign);
        let delta = db - da;
        slope = slope.min(if sign == Sign::Positive {
            delta.lo()
        } else {
            -delta.hi()
        });
        for k in 0..3 {
            ruling_bound[k] = ruling_bound[k].max(absolute(
                Interval::point(coordinates(b)[k]) - Interval::point(coordinates(a)[k]),
            ));
            world_spacing[k] = world_spacing[k].max(spacing(coordinates(a)[k], coordinates(b)[k]));
        }
        a_points.push(a);
        b_points.push(b);
        source_weights.push(wa);
    }
    let a = BCurve::new(
        degree,
        knots.to_vec(),
        a_points,
        Some(source_weights.clone()),
    )?;
    let b = BCurve::new(degree, knots.to_vec(), b_points, Some(source_weights))?;
    let source_breaks = breaks(&a)?;
    if source_breaks.len() < 3 {
        return Err(refused("multi-span entry needs at least two source spans"));
    }
    crate::patch::certify_graph(s, p)?;
    let (curve, uv) = construct(&a, &b, &source_breaks, linear, linear_u, plane)?;
    certify(
        &a,
        &b,
        &curve,
        &uv,
        &source_breaks,
        linear,
        linear_u,
        plane,
        slope,
        ruling_bound,
        &mut world_spacing,
        p,
    )?;
    let linear_coordinate = |point: Point3| if linear_u { point.x } else { point.y };
    let fraction = |point: Point3| (linear_coordinate(point) - linear.lo) / (linear.hi - linear.lo);
    let fractions = [
        fraction(uv.points()[0]),
        fraction(uv.points()[uv.points().len() - 1]),
    ];
    Ok((
        BoundaryEdge {
            curve,
            uv: Some(uv),
        },
        fractions,
    ))
}

#[allow(clippy::too_many_arguments)]
fn certify(
    a: &BCurve,
    b: &BCurve,
    curve: &BCurve,
    uv: &BCurve,
    source_breaks: &[f64],
    linear: ParamRange,
    linear_u: bool,
    plane: Plane,
    slope: f64,
    ruling_bound: [f64; 3],
    world_spacing: &mut [f64; 3],
    p: &Precision,
) -> Result<()> {
    let degree = a.degree();
    let q = uv.degree();
    let actual_breaks = breaks(a)?;
    let n = actual_breaks.len() - 1;
    let identical =
        |a: &[f64], b: &[f64]| a.len() == b.len() && a.iter().zip(b).all(|(&a, &b)| same(a, b));
    if !identical(source_breaks, &actual_breaks)
        || b.degree() != degree
        || !identical(a.knots(), b.knots())
        || (0..a.points().len()).any(|i| !same(weight(a, i), weight(b, i)))
        || curve.degree() != degree
        || q != degree + 1
        || !identical(curve.knots(), &c0_knots(degree, n))
        || !identical(uv.knots(), &c0_knots(q, n))
    {
        return Err(refused(
            "multi-span certificate requires paired canonical stored bases",
        ));
    }
    let reference = a.points()[0];
    let aa = raw_spans(a, reference)?;
    let bb = raw_spans(b, reference)?;
    let xx = raw_spans(curve, reference)?;
    // UV is translated by the source's global profile and linear origins.
    let uv_reference = if linear_u {
        Point3::new(linear.lo, a.range().lo, 0.0)
    } else {
        Point3::new(a.range().lo, linear.lo, 0.0)
    };
    let vv = raw_spans(uv, uv_reference)?;
    if [aa.len(), bb.len(), xx.len(), vv.len()]
        .iter()
        .any(|&len| len != n)
    {
        return Err(refused("multi-span stored bases disagree"));
    }
    let (profile_axis, linear_axis) = if linear_u { (1, 0) } else { (0, 1) };
    let inv_linear = reciprocal(Interval::point(linear.hi) - Interval::point(linear.lo))?;
    let linear_quantization = Interval::point(spacing(linear.lo, linear.hi)) * inv_linear;
    let max_source_weight = (0..a.points().len())
        .map(|i| weight(a, i))
        .fold(0.0, f64::max);
    let min_source_weight = (0..a.points().len())
        .map(|i| weight(a, i))
        .fold(f64::INFINITY, f64::min);
    let ratio =
        Interval::point(max_source_weight) * reciprocal(Interval::point(min_source_weight))?;
    let mut min_span = f64::INFINITY;
    for pair in source_breaks.windows(2) {
        min_span = min_span.min((Interval::point(pair[1]) - Interval::point(pair[0])).lo());
    }
    let profile_quantization = Interval::point(spacing(a.range().lo, a.range().hi))
        * reciprocal(Interval::point(min_span))?;
    let common_quantization = Interval::point(spacing(0.0, n as f64));
    let mut source_derivative = [Interval::point(0.0); 3];
    let mut output_derivative = [Interval::point(0.0); 3];
    let max_output = (0..curve.points().len())
        .map(|i| weight(curve, i))
        .fold(0.0, f64::max);
    let min_output = (0..curve.points().len())
        .map(|i| weight(curve, i))
        .fold(f64::INFINITY, f64::min);
    let output_ratio = Interval::point(max_output) * reciprocal(Interval::point(min_output))?;
    for k in 0..3 {
        let lo = a
            .points()
            .iter()
            .chain(b.points())
            .map(|&v| coordinates(v)[k])
            .fold(f64::INFINITY, f64::min);
        let hi = a
            .points()
            .iter()
            .chain(b.points())
            .map(|&v| coordinates(v)[k])
            .fold(f64::NEG_INFINITY, f64::max);
        source_derivative[k] = Interval::point((2 * degree) as f64)
            * ratio
            * (Interval::point(hi) - Interval::point(lo));
        let lo = curve
            .points()
            .iter()
            .map(|&v| coordinates(v)[k])
            .fold(f64::INFINITY, f64::min);
        let hi = curve
            .points()
            .iter()
            .map(|&v| coordinates(v)[k])
            .fold(f64::NEG_INFINITY, f64::max);
        output_derivative[k] = Interval::point((2 * degree) as f64)
            * output_ratio
            * (Interval::point(hi) - Interval::point(lo));
        world_spacing[k] = world_spacing[k].max(spacing(lo, hi));
    }
    let normal = plane.frame().z();
    let plane_error = curve
        .points()
        .iter()
        .map(|&point| absolute(dot_delta(point, plane.frame().origin(), normal)))
        .fold(0.0, f64::max);
    for span in 0..n {
        let columns =
            |data: &[[Interval; 4]], axis: usize| data.iter().map(|c| c[axis]).collect::<Vec<_>>();
        let ws = columns(&aa[span], 3);
        let wx = columns(&xx[span], 3);
        let wv = columns(&vv[span], 3);
        let min_s = positive_min(&ws)?;
        let min_x = positive_min(&wx)?;
        let min_v = positive_min(&wv)?;
        let inverse_product =
            reciprocal(Interval::point(min_s) * Interval::point(min_x) * Interval::point(min_v))?;
        let inverse_uv = reciprocal(Interval::point(min_v))?;
        let source_width =
            Interval::point(source_breaks[span + 1]) - Interval::point(source_breaks[span]);
        let inverse_width = reciprocal(source_width)?;
        let offset = Interval::point(source_breaks[span]) - Interval::point(a.range().lo);
        let profile: Vec<_> = vv[span]
            .iter()
            .map(|v| (v[profile_axis] - offset * v[3]) * inverse_width)
            .collect();
        let linear_numerator: Vec<_> = vv[span]
            .iter()
            .map(|v| v[linear_axis] * inv_linear)
            .collect();
        // Exact stored endpoint ordinates plus positive derivative keep the
        // complete UV image inside this ORIGINAL source span (including seams).
        if !same(
            coordinates(uv.points()[span * q])[profile_axis],
            source_breaks[span],
        ) || !same(
            coordinates(uv.points()[(span + 1) * q])[profile_axis],
            source_breaks[span + 1],
        ) {
            return Err(refused(
                "multi-span UV does not preserve source seam endpoints",
            ));
        }
        let increasing = difference(
            &product(&derivative(&profile), &wv)?,
            &product(&profile, &derivative(&wv))?,
        );
        if increasing.iter().any(|v| v.sign() != Some(Sign::Positive)) {
            return Err(refused("multi-span UV monotonicity is unresolved"));
        }
        if linear_numerator.iter().zip(&wv).any(|(&v, &w)| {
            v.sign() != Some(Sign::Positive) || (w - v).sign() != Some(Sign::Positive)
        }) {
            return Err(refused("multi-span UV leaves a ruling"));
        }
        let identity = difference(
            &product(&profile, &[Interval::point(1.0); 2])?,
            &product(&wv, &[Interval::point(0.0), Interval::point(1.0)])?,
        );
        let profile_error = Interval::point(max_abs(&identity)) * inverse_uv + profile_quantization;
        let mut lift = [Interval::point(0.0); 3];
        for k in 0..3 {
            let ha = columns(&aa[span], k);
            let hb = columns(&bb[span], k);
            let hx = columns(&xx[span], k);
            // X*Ws*Wuv - [A*Wuv + (B-A)*V]*Wx, all relative to the same
            // spatial reference. None of the nominal extraction is trusted.
            let left = product(&product(&hx, &ws)?, &wv)?;
            let mut source = product(&ha, &wv)?;
            for (value, ruling) in source
                .iter_mut()
                .zip(product(&difference(&hb, &ha), &linear_numerator)?)
            {
                *value = *value + ruling;
            }
            let residual = difference(&left, &product(&source, &wx)?);
            lift[k] = Interval::point(max_abs(&residual)) * inverse_product
                + profile_error * source_derivative[k]
                + linear_quantization * Interval::point(ruling_bound[k])
                + common_quantization * output_derivative[k]
                + Interval::point(world_spacing[k]) * Interval::point(2.0);
        }
        let mut source_plane_error = Interval::point(plane_error) + Interval::point(p.solver());
        for k in 0..3 {
            source_plane_error = source_plane_error
                + lift[k] * Interval::point([normal.x, normal.y, normal.z][k].abs());
        }
        let root_error = source_plane_error * reciprocal(Interval::point(slope))?;
        let mut squared = Interval::point(0.0);
        for k in 0..3 {
            squared = squared + (lift[k] + root_error * Interval::point(ruling_bound[k])).sqr();
        }
        if (squared - Interval::point(p.default_local()).sqr()).sign() != Some(Sign::Negative) {
            return Err(refused(
                "multi-span section lacks a whole-branch spatial certificate",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]
    use super::*;

    fn fixture(domain: [f64; 2], scale: f64) -> (BCurve, BCurve, Vec<f64>, Plane) {
        let a = BCurve::new(
            2,
            vec![
                domain[0], domain[0], domain[0], domain[1], domain[1], domain[1],
            ],
            vec![
                Point3::ORIGIN,
                Point3::new(1.0, 0.0, 0.12),
                Point3::new(2.0, 0.0, 0.0),
            ],
            Some(vec![scale, scale * 0.8, scale]),
        )
        .unwrap();
        let span = domain[1] - domain[0];
        let a = a
            .insert_knots(&[
                domain[0] + span * 0.1,
                domain[0] + span * 0.4,
                domain[0] + span * 0.4,
                domain[0] + span * 0.75,
            ])
            .unwrap();
        // Restore a common scale AFTER nominal construction to exercise the
        // exact raw-input weight normalization used by the certificate.
        let a = BCurve::new(
            a.degree(),
            a.knots().to_vec(),
            a.points().to_vec(),
            Some(a.weights().unwrap().iter().map(|w| w * scale).collect()),
        )
        .unwrap();
        let b = BCurve::new(
            a.degree(),
            a.knots().to_vec(),
            a.points().iter().map(|&p| p + Vec3::Y * 2.0).collect(),
            a.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        let source_breaks = breaks(&a).unwrap();
        let frame = Frame::new(
            Point3::new(0.0, 0.7, 0.0),
            Vec3::new(-0.1, 1.0, -0.1),
            Vec3::X,
            &Precision::DEFAULT,
        )
        .unwrap();
        (a, b, source_breaks, Plane::new(frame))
    }

    fn check(
        a: &BCurve,
        b: &BCurve,
        curve: &BCurve,
        uv: &BCurve,
        knots: &[f64],
        plane: Plane,
    ) -> Result<()> {
        certify(
            a,
            b,
            curve,
            uv,
            knots,
            ParamRange::bounded(0.0, 1.0),
            false,
            plane,
            1.9,
            [0.0, 2.0, 0.0],
            &mut [0.0; 3],
            &Precision::DEFAULT,
        )
    }

    #[test]
    fn final_storage_certifies_nonuniform_source_and_common_weight_scales() {
        for scale in [2.0_f64.powi(-400), 1.0, 2.0_f64.powi(400)] {
            let (a, b, knots, plane) = fixture([-2.0, 5.0], scale);
            let (curve, uv) =
                construct(&a, &b, &knots, ParamRange::bounded(0.0, 1.0), false, plane).unwrap();
            check(&a, &b, &curve, &uv, &knots, plane).unwrap();
            assert_eq!(curve.range(), ParamRange::bounded(0.0, 4.0));
            for (i, &knot) in knots.iter().enumerate() {
                assert_eq!(uv.points()[i * uv.degree()].x.to_bits(), knot.to_bits());
            }
            assert_eq!(
                curve.knots(),
                curve.reversed().unwrap().reversed().unwrap().knots()
            );
            assert_eq!(
                uv.knots(),
                uv.reversed().unwrap().reversed().unwrap().knots()
            );
        }
    }

    #[test]
    fn independently_rejects_corrupt_xyz_uv_and_seam_storage() {
        let (a, b, knots, plane) = fixture([-2.0, 5.0], 1.0);
        let (curve, uv) =
            construct(&a, &b, &knots, ParamRange::bounded(0.0, 1.0), false, plane).unwrap();
        check(&a, &b, &curve, &uv, &knots, plane).unwrap();
        let mut points = curve.points().to_vec();
        points[3].x += 0.1;
        let bad = BCurve::new(
            curve.degree(),
            curve.knots().to_vec(),
            points,
            curve.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        assert!(check(&a, &b, &bad, &uv, &knots, plane).is_err());
        let mut points = uv.points().to_vec();
        points[2].y += 0.1;
        let bad = BCurve::new(
            uv.degree(),
            uv.knots().to_vec(),
            points,
            uv.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        assert!(check(&a, &b, &curve, &bad, &knots, plane).is_err());
        let shifted_knots = curve
            .knots()
            .iter()
            .map(|&t| if same(t, 1.0) { 0.75 } else { t })
            .collect();
        let bad = BCurve::new(
            curve.degree(),
            shifted_knots,
            curve.points().to_vec(),
            curve.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        assert!(check(&a, &b, &bad, &uv, &knots, plane).is_err());
        let mut points = uv.points().to_vec();
        points[uv.degree()].x = points[uv.degree()].x.next_up();
        let bad = BCurve::new(
            uv.degree(),
            uv.knots().to_vec(),
            points,
            uv.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        assert!(check(&a, &b, &curve, &bad, &knots, plane).is_err());
    }

    #[test]
    fn distinct_rational_spans_are_certified_without_a_single_bezier_ancestor() {
        let (_, _, _, plane) = fixture([-2.0, 5.0], 1.0);
        let a = BCurve::new(
            2,
            vec![-2.0, -2.0, -2.0, -1.3, -0.2, -0.2, 3.25, 5.0, 5.0, 5.0],
            [
                (0.0, 0.0),
                (0.2, 0.1),
                (0.7, 0.0),
                (1.1, -0.1),
                (1.4, 0.0),
                (1.7, 0.1),
                (2.0, 0.0),
            ]
            .into_iter()
            .map(|(x, z)| Point3::new(x, 0.0, z))
            .collect(),
            Some(vec![1.0, 0.8, 1.2, 0.9, 0.7, 1.1, 1.0]),
        )
        .unwrap();
        let b = BCurve::new(
            a.degree(),
            a.knots().to_vec(),
            a.points().iter().map(|&p| p + Vec3::Y * 2.0).collect(),
            a.weights().map(<[f64]>::to_vec),
        )
        .unwrap();
        let knots = breaks(&a).unwrap();
        let (curve, uv) =
            construct(&a, &b, &knots, ParamRange::bounded(0.0, 1.0), false, plane).unwrap();
        check(&a, &b, &curve, &uv, &knots, plane).unwrap();
    }

    #[test]
    fn huge_source_parameter_quantization_is_not_hidden_by_integer_edge_domain() {
        let (a, b, knots, plane) = fixture([1.0e15, 1.0e15 + 8.0], 1.0);
        let (curve, uv) =
            construct(&a, &b, &knots, ParamRange::bounded(0.0, 1.0), false, plane).unwrap();
        assert!(check(&a, &b, &curve, &uv, &knots, plane).is_err());
    }
}
