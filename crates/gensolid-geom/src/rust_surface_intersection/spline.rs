//! Original homogeneous algebra for rational spline extrusion / affine plane sections.
//! No foreign kernel source is used. These candidates do not certify B-rep topology.
use super::*;
use crate::ParamRange;

struct Extrusion {
    profile: BCurve,
    direction: Vec3,
    linear: ParamRange,
    linear_u: bool,
    error: f64,
}

fn extrusion(s: &BSurface, linear_u: bool, p: &Precision) -> Result<Extrusion> {
    let (nu, nv) = s.counts();
    let (pu, pv) = s.degrees();
    let (uk, vk) = s.knots();
    let (u, v) = s.ranges();
    let (degree, count, knots, linear, linear_knots) = if linear_u {
        if pu != 1 || nu != 2 {
            return Err(Error::unimplemented(
                "spline section needs one linear direction",
            ));
        }
        (pv, nv, vk, u, uk)
    } else {
        if pv != 1 || nv != 2 {
            return Err(Error::unimplemented(
                "spline section needs one linear direction",
            ));
        }
        (pu, nu, uk, v, vk)
    };
    if degree > 16 {
        return Err(Error::unimplemented("spline section degree exceeds 16"));
    }
    clamped(linear_knots, 1, linear.lo, linear.hi)?;
    let indices = |i| {
        if linear_u {
            (i, nv + i)
        } else {
            (i * nv, i * nv + 1)
        }
    };
    let (a, b) = indices(0);
    let direction = s.points()[b] - s.points()[a];
    let mut error = 0.0_f64;
    let mut points = Vec::with_capacity(count);
    let mut weights = Vec::with_capacity(count);
    for i in 0..count {
        let (a, b) = indices(i);
        let wa = s.weights().map_or(1.0, |w| w[a]);
        let wb = s.weights().map_or(1.0, |w| w[b]);
        // Equal paired weights are a representation identity: varying them
        // along a ruling changes its parameter and invalidates P(u) + v D.
        if !wa.total_cmp(&wb).is_eq() {
            return Err(Error::unimplemented(
                "ruling weights do not describe a translation",
            ));
        }
        error = error.max((s.points()[b] - s.points()[a] - direction).length());
        points.push(s.points()[a]);
        weights.push(wa);
    }
    if error > p.solver() {
        return Err(Error::unimplemented("spline rulings are not translations"));
    }
    let profile = BCurve::new(degree, knots.to_vec(), points, Some(weights))?;
    let r = profile.range();
    if knots[..=degree].iter().any(|k| !k.total_cmp(&r.lo).is_eq())
        || knots[count..].iter().any(|k| !k.total_cmp(&r.hi).is_eq())
        || knots
            .windows(degree + 1)
            .any(|ks| ks[0] > r.lo && ks[0] < r.hi && ks[0].total_cmp(&ks[degree]).is_eq())
    {
        return Err(Error::unimplemented(
            "spline section needs clamped, continuous profile spans",
        ));
    }
    Ok(Extrusion {
        profile,
        direction,
        linear,
        linear_u,
        error,
    })
}

fn strict_unit(x: Interval) -> bool {
    x.sign() == Some(Sign::Positive) && (Interval::point(1.0) - x).sign() == Some(Sign::Positive)
}

fn plane_uv(point: Point3, plane: &AffinePatch) -> Result<(Interval, Interval)> {
    let delta = interval_delta(point, plane.origin);
    let u = cross_dot(delta, plane.v, plane.normal) * plane.inverse_area;
    let v = cross_dot(plane.u, delta, plane.normal) * plane.inverse_area;
    if !strict_unit(u) || !strict_unit(v) {
        return Err(Error::unimplemented(
            "spline section is not wholly inside the bounded tool",
        ));
    }
    Ok((u, v))
}

fn bounds(profile: &BCurve) -> (f64, f64, f64) {
    let weights = profile.weights();
    let max = weights.map_or(1.0, |w| w.iter().copied().fold(0.0, f64::max));
    let min = weights.map_or(1.0, |w| w.iter().copied().fold(f64::INFINITY, f64::min));
    let reference = profile.points()[0];
    let diameter = profile
        .points()
        .iter()
        .map(|q| q.distance(reference))
        .fold(0.0, f64::max)
        * 2.0;
    (max, min, diameter)
}

/// Construct u(t)=t as (t W(t))/W(t), not Greville ordinates with rational weights.
/// Multiplication by the affine parameter raises the homogeneous degree by one.
fn first_uv(
    profile: &BCurve,
    fractions: &[Interval],
    extrusion: &Extrusion,
) -> Result<(BCurve, f64, f64)> {
    let degree = profile.degree();
    let r = profile.range();
    let (max_weight, _, _) = bounds(profile);
    let weights: Vec<_> = (0..=degree)
        .map(|i| profile.weights().map_or(1.0, |w| w[i]) / max_weight)
        .collect();
    let inv_degree = inverse(Interval::point((degree + 1) as f64))?;
    let inv_profile_span = inverse(Interval::point(r.hi) - Interval::point(r.lo))?;
    let linear = extrusion.linear;
    let inv_linear_span = inverse(Interval::point(linear.hi) - Interval::point(linear.lo))?;
    let mut points = Vec::with_capacity(degree + 2);
    let mut output_weights = Vec::with_capacity(degree + 2);
    let mut profile_error = 0.0_f64;
    let mut linear_error = 0.0_f64;
    let mut minimum_weight = f64::INFINITY;
    for i in 0..=degree + 1 {
        let alpha = Interval::point(i as f64) * inv_degree;
        let left = if i <= degree {
            Interval::point(weights[i]) * (Interval::point(1.0) - alpha)
        } else {
            Interval::point(0.0)
        };
        let right = if i > 0 {
            Interval::point(weights[i - 1]) * alpha
        } else {
            Interval::point(0.0)
        };
        let weight = left + right;
        let reciprocal = inverse(weight)?;
        let numerator = if i <= degree {
            left * fractions[i]
        } else {
            Interval::point(0.0)
        } + if i > 0 {
            right * fractions[i - 1]
        } else {
            Interval::point(0.0)
        };
        let profile_t = if i == 0 {
            r.lo
        } else if i == degree + 1 {
            r.hi
        } else {
            r.lo + (r.hi - r.lo) * midpoint(right * reciprocal)
        };
        let linear_t = linear.lo + (linear.hi - linear.lo) * midpoint(numerator * reciprocal);
        let stored_weight = midpoint(weight);
        if !profile_t.is_finite()
            || !(profile_t >= r.lo && profile_t <= r.hi)
            || !(linear_t > linear.lo && linear_t < linear.hi)
            || !stored_weight.is_finite()
            || stored_weight <= 0.0
        {
            return Err(Error::unimplemented(
                "rational section UV is not representable",
            ));
        }
        // Bound homogeneous coefficient storage errors and denominator elevation.
        // The desired normalized coordinates lie in [0,1], so denominator error
        // contributes at most its magnitude when comparing the two quotients.
        let w = Interval::point(stored_weight);
        let denominator_error = magnitude(w - weight);
        let represented_profile =
            (Interval::point(profile_t) - Interval::point(r.lo)) * inv_profile_span * w;
        let represented_linear =
            (Interval::point(linear_t) - Interval::point(linear.lo)) * inv_linear_span * w;
        profile_error =
            profile_error.max(magnitude(represented_profile - right) + denominator_error);
        linear_error =
            linear_error.max(magnitude(represented_linear - numerator) + denominator_error);
        minimum_weight = minimum_weight.min(stored_weight);
        points.push(if extrusion.linear_u {
            Point3::new(linear_t, profile_t, 0.0)
        } else {
            Point3::new(profile_t, linear_t, 0.0)
        });
        output_weights.push(stored_weight);
    }
    let reciprocal = magnitude(inverse(Interval::point(minimum_weight))?);
    Ok((
        BCurve::new(
            degree + 1,
            [vec![r.lo; degree + 2], vec![r.hi; degree + 2]].concat(),
            points,
            Some(output_weights),
        )?,
        profile_error * reciprocal,
        linear_error * reciprocal,
    ))
}

fn segment(
    profile: &BCurve,
    extrusion: &Extrusion,
    tool: &BSurface,
    plane: &AffinePatch,
    reciprocal: Interval,
    p: &Precision,
) -> Result<SurfaceIntersectionCandidate> {
    let (tool_u, tool_v) = tool.ranges();
    let mut xyz = Vec::with_capacity(profile.points().len());
    let mut tool_uv = Vec::with_capacity(profile.points().len());
    let mut fractions = Vec::with_capacity(profile.points().len());
    let mut error = extrusion.error;
    for &a in profile.points() {
        let fraction = -dot_delta(a, plane.origin, plane.normal) * reciprocal;
        if !strict_unit(fraction) {
            return Err(Error::unimplemented(
                "spline section requires a complete interior branch",
            ));
        }
        let point = a + extrusion.direction * midpoint(fraction);
        let (u, v) = plane_uv(point, plane)?;
        let uv = Point3::new(
            tool_u.lo + (tool_u.hi - tool_u.lo) * midpoint(u),
            tool_v.lo + (tool_v.hi - tool_v.lo) * midpoint(v),
            0.0,
        );
        if !(uv.x > tool_u.lo && uv.x < tool_u.hi && uv.y > tool_v.lo && uv.y < tool_v.hi) {
            return Err(Error::unimplemented(
                "rational section tool UV is not representable",
            ));
        }
        let lift = interval_delta(point, a);
        let d = extrusion.direction;
        let residual = [
            lift[0] - Interval::point(d.x) * fraction,
            lift[1] - Interval::point(d.y) * fraction,
            lift[2] - Interval::point(d.z) * fraction,
        ];
        error = error.max(
            magnitude(residual[0])
                .hypot(magnitude(residual[1]))
                .hypot(magnitude(residual[2])),
        );
        error = error.max(tool.eval(uv.x, uv.y).point.distance(point));
        fractions.push(fraction);
        xyz.push(point);
        tool_uv.push(uv);
    }
    let (uv, profile_error, linear_error) = first_uv(profile, &fractions, extrusion)?;
    let (max_weight, min_weight, diameter) = bounds(profile);
    let ratio = magnitude(Interval::point(max_weight) * inverse(Interval::point(min_weight))?);
    // For positive rational Bernstein weights, sum |N_i'| <= 2p and
    // |P_i-C(t)| <= control-hull diameter. This bounds the entire profile speed.
    let derivative = 2.0 * profile.degree() as f64 * ratio * diameter;
    let r = profile.range();
    let quantization = parameter_spacing(r.lo, r.hi)
        * magnitude(inverse(Interval::point(r.hi) - Interval::point(r.lo))?)
        * derivative
        + parameter_spacing(extrusion.linear.lo, extrusion.linear.hi)
            * magnitude(inverse(
                Interval::point(extrusion.linear.hi) - Interval::point(extrusion.linear.lo),
            )?)
            * extrusion.direction.length();
    let tool_quantization = parameter_spacing(tool_u.lo, tool_u.hi)
        * magnitude(inverse(
            Interval::point(tool_u.hi) - Interval::point(tool_u.lo),
        )?)
        * (tool.points()[2] - tool.points()[0]).length()
        + parameter_spacing(tool_v.lo, tool_v.hi)
            * magnitude(inverse(
                Interval::point(tool_v.hi) - Interval::point(tool_v.lo),
            )?)
            * (tool.points()[1] - tool.points()[0]).length();
    // Bezier extraction uses homogeneous knot insertion. Account conservatively
    // for its floating point workspace; this is a provider diagnostic, not the
    // modelling layer's independent whole-branch certificate.
    let scale = extrusion
        .profile
        .points()
        .iter()
        .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
        .fold(1.0, f64::max);
    let extraction_roundoff = (scale.next_up() - scale)
        * ((profile.degree() + 1) * (profile.degree() + 1) * 32) as f64
        * ratio
        * ratio;
    let deviation = error
        + extrusion.error
        + profile_error * derivative
        + linear_error * extrusion.direction.length()
        + quantization
        + tool_quantization
        + extraction_roundoff;
    if !deviation.is_finite() || deviation > p.linear() {
        return Err(Error::unimplemented(
            "rational section conditioning or UV quantization exceeds precision",
        ));
    }
    Ok(SurfaceIntersectionCandidate {
        curve: BCurve::new(
            profile.degree(),
            profile.knots().to_vec(),
            xyz,
            profile.weights().map(<[f64]>::to_vec),
        )?,
        first_uv: uv,
        second_uv: BCurve::new(
            profile.degree(),
            profile.knots().to_vec(),
            tool_uv,
            profile.weights().map(<[f64]>::to_vec),
        )?,
        reported_deviation: deviation,
    })
}

fn spline_section(
    first: &BSurface,
    second: &BSurface,
    linear_u: bool,
    p: &Precision,
) -> Result<Vec<SurfaceIntersectionCandidate>> {
    let plane = affine_patch(second, p)?;
    let extrusion = extrusion(first, linear_u, p)?;
    let slab = second
        .points()
        .iter()
        .map(|q| dot_delta(*q, plane.origin, plane.normal))
        .fold(Interval::point(0.0), Interval::hull);
    let distances: Vec<_> = first
        .points()
        .iter()
        .map(|q| dot_delta(*q, plane.origin, plane.normal))
        .collect();
    if distances
        .iter()
        .all(|d| (*d - Interval::point(slab.hi())).sign() == Some(Sign::Positive))
        || distances
            .iter()
            .all(|d| (*d - Interval::point(slab.lo())).sign() == Some(Sign::Negative))
    {
        return Ok(Vec::new());
    }
    let reciprocal = inverse(dot(extrusion.direction, plane.normal))?;
    // Require the whole source profile's control hull to cross internally. Do not
    // silently return selected spans of a partially clipped or tangent branch.
    for &a in extrusion.profile.points() {
        let fraction = -dot_delta(a, plane.origin, plane.normal) * reciprocal;
        if !strict_unit(fraction) {
            return Err(Error::unimplemented(
                "spline section needs unclipped complete crossing",
            ));
        }
        plane_uv(a + extrusion.direction * midpoint(fraction), &plane)?;
    }
    let range = extrusion.profile.range();
    let mut result = Vec::new();
    for span in extrusion
        .profile
        .knots()
        .windows(2)
        .filter(|k| k[0] < k[1] && k[0] >= range.lo && k[1] <= range.hi)
    {
        let profile = extrusion.profile.trim([span[0], span[1]])?;
        result.push(segment(
            &profile, &extrusion, second, &plane, reciprocal, p,
        )?);
    }
    Ok(result)
}

/// Paired algebraic candidates for positive-weight, continuous clamped spline
/// extrusions and bounded affine plane patches. One candidate is returned per
/// nonzero profile knot span, preserving that span's original parameter range.
/// Rational profile UV uses a degree-elevated homogeneous parameter product.
/// Both operand orders and extrusion axes are supported. Partial, tangent,
/// overlapping, discontinuous and unresolved ill-conditioned cases are refused.
/// These curves are not a certificate that a modelling operation is supported:
/// rational isoparametric cuts have an independent modelling certificate, while
/// general rational oblique cuts remain outside that certificate's scope.
pub fn intersect_spline_extrusion_plane(
    first: &BSurface,
    second: &BSurface,
    p: &Precision,
) -> Result<Vec<SurfaceIntersectionCandidate>> {
    if let Ok(result) = intersect_extrusion_plane(first, second, p) {
        return Ok(result);
    }
    for linear_u in [false, true] {
        if let Ok(result) = spline_section(first, second, linear_u, p) {
            return Ok(result);
        }
    }
    for linear_u in [false, true] {
        if let Ok(mut result) = spline_section(second, first, linear_u, p) {
            for branch in &mut result {
                std::mem::swap(&mut branch.first_uv, &mut branch.second_uv);
            }
            return Ok(result);
        }
    }
    Err(Error::unimplemented(
        "Rust spline sections need complete rational extrusion / affine plane branches",
    ))
}
