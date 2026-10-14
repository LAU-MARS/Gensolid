//! Polynomial extrusion/plane sections with whole-curve error bounds.
//! For S(t,v)=P(t)+v D, plane intersection is the affine projection of P along D.
//! Bernstein coefficients therefore construct both the spatial curve and its UV
//! curve without tracing or tessellation. Small extrusion representation errors
//! are bounded separately; interval transversality proves one complete branch.
use super::*;
use crate::boundary::BoundaryEdge;
use gensolid_geom::BCurve;

fn vector(p: Point3) -> [Interval; 3] {
    [
        Interval::point(p.x),
        Interval::point(p.y),
        Interval::point(p.z),
    ]
}
fn absolute(a: Interval) -> f64 {
    a.lo().abs().max(a.hi().abs())
}
fn reciprocal(a: Interval) -> Result<Interval> {
    if a.sign() != Some(Sign::Positive) && a.sign() != Some(Sign::Negative) {
        return Err(Error::unimplemented(
            "section denominator is not separated from zero",
        ));
    }
    Interval::new((1. / a.hi()).next_down(), (1. / a.lo()).next_up())
        .ok_or_else(|| Error::unimplemented("section denominator overflow"))
}

pub(super) fn section(
    s: &BSurface,
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
            return Err(Error::unimplemented(
                "section requires one linear direction",
            ));
        }
        (pv, nv, vk, v, u, uk)
    } else {
        if pv != 1 || nv != 2 {
            return Err(Error::unimplemented(
                "section requires one linear direction",
            ));
        }
        (pu, nu, uk, u, v, vk)
    };
    if count != degree + 1
        || degree > 16
        || knots[..count]
            .iter()
            .any(|&k| k.to_bits() != domain.lo.to_bits())
        || knots[count..]
            .iter()
            .any(|&k| k.to_bits() != domain.hi.to_bits())
        || linear_knots[0].to_bits() != linear.lo.to_bits()
        || linear_knots[3].to_bits() != linear.hi.to_bits()
        || s.weights()
            .is_some_and(|w| w.windows(2).any(|w| w[0].to_bits() != w[1].to_bits()))
    {
        return Err(Error::unimplemented(
            "non-isoparametric section requires a polynomial Bezier extrusion",
        ));
    }
    // The projection proof also rules out self-intersections of the target patch.
    crate::patch::certify_graph(s, p)?;
    let pair = |i| {
        if linear_u {
            (s.points()[i], s.points()[nv + i])
        } else {
            (s.points()[i * nv], s.points()[i * nv + 1])
        }
    };
    let (start, end) = pair(0);
    let direction = end - start;
    let d = [
        Interval::point(direction.x),
        Interval::point(direction.y),
        Interval::point(direction.z),
    ];
    let normal = plane.frame().z();
    let origin = plane.frame().origin();
    let n = [normal.x, normal.y, normal.z];
    let denominator =
        d[0] * Interval::point(n[0]) + d[1] * Interval::point(n[1]) + d[2] * Interval::point(n[2]);
    let inverse = reciprocal(denominator)?;
    let inverse_degree = reciprocal(Interval::point(degree as f64))?;
    let span = Interval::point(domain.hi) - Interval::point(domain.lo);
    let inverse_span = reciprocal(span)?;
    let linear_span = Interval::point(linear.hi) - Interval::point(linear.lo);
    let inverse_linear = reciprocal(linear_span)?;
    let mut extrusion_error = [0.0f64; 3];
    let mut lift_error = [0.0f64; 3];
    let mut derivative = [0.0f64; 3];
    let mut ruling_bound = [0.0f64; 3];
    let mut parameter_error = 0.0f64;
    let mut plane_error = 0.0f64;
    let mut slope = f64::INFINITY;
    let mut orientation = None;
    let mut points = Vec::new();
    let mut uv_points = Vec::new();
    let mut fractions = Vec::new();
    for i in 0..count {
        let (a, b) = pair(i);
        let a_vector = vector(a);
        let b_vector = vector(b);
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
                "section is not a full transverse branch",
            ));
        };
        if orientation.is_some_and(|old| old != sign) {
            return Err(Error::unimplemented("section changes crossing orientation"));
        }
        orientation = Some(sign);
        let difference = db - da;
        slope = slope.min(if sign == Sign::Positive {
            difference.lo()
        } else {
            -difference.hi()
        });
        let t_interval = (Interval::point(0.) - da) * inverse;
        if t_interval.sign() != Some(Sign::Positive)
            || (Interval::point(1.) - t_interval).sign() != Some(Sign::Positive)
        {
            return Err(Error::unimplemented(
                "section parameter leaves the target rectangle",
            ));
        }
        let t = t_interval.lo() + (t_interval.hi() - t_interval.lo()) * 0.5;
        let point = a + direction * t;
        let profile_parameter = if i == 0 {
            domain.lo
        } else if i + 1 == count {
            domain.hi
        } else {
            domain.lo + (domain.hi - domain.lo) * (i as f64 / degree as f64)
        };
        let linear_parameter = linear.lo + (linear.hi - linear.lo) * t;
        let uv = if linear_u {
            Point3::new(linear_parameter, profile_parameter, 0.)
        } else {
            Point3::new(profile_parameter, linear_parameter, 0.)
        };
        let ideal_parameter =
            Interval::point(domain.lo) + span * Interval::point(i as f64) * inverse_degree;
        parameter_error = parameter_error.max(absolute(
            Interval::point(profile_parameter) - ideal_parameter,
        ));
        let normalized =
            (Interval::point(linear_parameter) - Interval::point(linear.lo)) * inverse_linear;
        if normalized.sign() != Some(Sign::Positive)
            || (Interval::point(1.) - normalized).sign() != Some(Sign::Positive)
        {
            return Err(Error::unimplemented("rounded section UV leaves its domain"));
        }
        let point_vector = vector(point);
        for k in 0..3 {
            let ruling = b_vector[k] - a_vector[k];
            ruling_bound[k] = ruling_bound[k].max(absolute(ruling));
            extrusion_error[k] = extrusion_error[k].max(absolute(ruling - d[k]));
            lift_error[k] =
                lift_error[k].max(absolute(point_vector[k] - a_vector[k] - d[k] * normalized));
            if i > 0 {
                let previous = pair(i - 1);
                for (before, after) in [(previous.0, a), (previous.1, b)] {
                    derivative[k] = derivative[k].max(absolute(
                        (vector(after)[k] - vector(before)[k])
                            * Interval::point(degree as f64)
                            * inverse_span,
                    ));
                }
            }
        }
        plane_error = plane_error.max(absolute(dot_delta(point, origin, normal)));
        points.push(point);
        uv_points.push(uv);
        fractions.push(t);
    }
    let mut lifted = [Interval::point(0.); 3];
    let mut source_plane_error = Interval::point(plane_error) + Interval::point(p.solver());
    for k in 0..3 {
        lifted[k] = Interval::point(lift_error[k])
            + Interval::point(extrusion_error[k])
            + Interval::point(parameter_error) * Interval::point(derivative[k]);
        source_plane_error = source_plane_error + lifted[k] * Interval::point(n[k].abs());
    }
    // A complete Hausdorff tube bound: move the lifted UV curve to the unique
    // true root along each ruling. Strictly monotone profile UV spans the whole
    // interval, so this bounds every point of the true branch as well.
    let root_parameter_error = source_plane_error * reciprocal(Interval::point(slope))?;
    let mut error_squared = Interval::point(0.);
    for k in 0..3 {
        let bound = lifted[k] + root_parameter_error * Interval::point(ruling_bound[k]);
        error_squared = error_squared + bound.sqr();
    }
    if (error_squared - Interval::point(p.default_local()).sqr()).sign() != Some(Sign::Negative) {
        return Err(Error::unimplemented(
            "section error exceeds the whole-branch tolerance",
        ));
    }
    let profile_coordinate = |p: &Point3| if linear_u { p.y } else { p.x };
    for pair in uv_points.windows(2) {
        if (Interval::point(profile_coordinate(&pair[1]))
            - Interval::point(profile_coordinate(&pair[0])))
        .sign()
            != Some(Sign::Positive)
        {
            return Err(Error::unimplemented("section UV is not certified monotone"));
        }
    }
    Ok((
        BoundaryEdge {
            curve: BCurve::new(degree, knots.to_vec(), points, None)?,
            uv: Some(BCurve::new(degree, knots.to_vec(), uv_points, None)?),
        },
        [fractions[0], fractions[count - 1]],
    ))
}
