//! Whole-domain graph bounds for continuous rational ruled NURBS sheets.
//! Source coefficients are extracted with outward rounding, not floating trim.
use crate::spline_bounds::homogeneous_spans;
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{BCurve, BSurface};
use gensolid_math::{Frame, Interval, Point3, Sign, Vec3};

fn refused() -> Error {
    Error::unimplemented("ruled surface graph Jacobian is not certified")
}

fn reciprocal(value: Interval) -> Result<Interval> {
    if !value.lo().is_finite() || !value.hi().is_finite() || value.sign() != Some(Sign::Positive) {
        return Err(refused());
    }
    Interval::new((1.0 / value.hi()).next_down(), (1.0 / value.lo()).next_up())
        .filter(|x| x.lo().is_finite() && x.hi().is_finite())
        .ok_or_else(refused)
}

fn hull(values: impl IntoIterator<Item = Interval>) -> Result<Interval> {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for value in values {
        if !value.lo().is_finite() || !value.hi().is_finite() {
            return Err(refused());
        }
        lo = lo.min(value.lo());
        hi = hi.max(value.hi());
    }
    Interval::new(lo, hi).ok_or_else(refused)
}

fn derivative(values: &[Interval]) -> Vec<Interval> {
    values
        .windows(2)
        .map(|w| (w[1] - w[0]) * Interval::point((values.len() - 1) as f64))
        .collect()
}

fn binomial(n: usize, k: usize) -> f64 {
    // Profile degree <=16: all products here have degree <=31. Integer
    // recurrence intermediates and final binomials are exactly f64-representable.
    let mut value = 1_u64;
    for i in 1..=k.min(n - k) {
        value = value * (n + 1 - i) as u64 / i as u64;
    }
    value as f64
}

fn product(a: &[Interval], b: &[Interval]) -> Result<Vec<Interval>> {
    let p = a.len() - 1;
    let q = b.len() - 1;
    if p + q > 31 {
        return Err(refused());
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

fn project(span: &[[Interval; 4]], axis: Vec3) -> Vec<Interval> {
    span.iter()
        .map(|h| {
            h[0] * Interval::point(axis.x)
                + h[1] * Interval::point(axis.y)
                + h[2] * Interval::point(axis.z)
        })
        .collect()
}

fn derivative_numerator(h: &[Interval], w: &[Interval]) -> Result<Vec<Interval>> {
    Ok(product(&derivative(h), w)?
        .into_iter()
        .zip(product(h, &derivative(w))?)
        .map(|(a, b)| a - b)
        .collect())
}

fn certify_direction(s: &BSurface, frame: Frame, linear_u: bool, p: &Precision) -> Result<()> {
    let (pu, pv) = s.degrees();
    let (nu, nv) = s.counts();
    let (uk, vk) = s.knots();
    let (u, v) = s.ranges();
    let (degree, count, knots, domain, linear, linear_knots) = if linear_u {
        if pu != 1 || nu != 2 {
            return Err(refused());
        }
        (pv, nv, vk, v, u, uk)
    } else {
        if pv != 1 || nv != 2 {
            return Err(refused());
        }
        (pu, nu, uk, u, v, vk)
    };
    if degree > 16
        || linear_knots.len() != 4
        || linear_knots[..2]
            .iter()
            .any(|k| k.to_bits() != linear.lo.to_bits())
        || linear_knots[2..]
            .iter()
            .any(|k| k.to_bits() != linear.hi.to_bits())
    {
        return Err(refused());
    }
    let mut a = Vec::with_capacity(count);
    let mut b = Vec::with_capacity(count);
    let mut weights = Vec::with_capacity(count);
    for i in 0..count {
        let (ia, ib) = if linear_u {
            (i, nv + i)
        } else {
            (i * nv, i * nv + 1)
        };
        let wa = s.weights().map_or(1.0, |w| w[ia]);
        if wa.to_bits() != s.weights().map_or(1.0, |w| w[ib]).to_bits() {
            return Err(refused());
        }
        a.push(s.points()[ia]);
        b.push(s.points()[ib]);
        weights.push(wa);
    }
    let a = BCurve::new(degree, knots.to_vec(), a, Some(weights.clone()))?;
    let b = BCurve::new(degree, knots.to_vec(), b, Some(weights))?;
    let reference: Point3 = s.points()[0];
    let aa = homogeneous_spans(&a, reference)?;
    let bb = homogeneous_spans(&b, reference)?;
    let spans = knots
        .windows(2)
        .filter(|k| k[0] < k[1] && k[0] >= domain.lo && k[1] <= domain.hi)
        .collect::<Vec<_>>();
    if aa.len() != bb.len() || aa.len() != spans.len() {
        return Err(refused());
    }
    let global_width = Interval::point(domain.hi) - Interval::point(domain.lo);
    for ((a, b), span) in aa.iter().zip(&bb).zip(spans) {
        // Equal original row weights imply the same exact W in both rows.
        // The surface is affine in the normalized ruling parameter. Its profile
        // derivative numerator is thus affine too: the two row bounds enclose
        // the complete ruling, including both one-sided limits at profile seams.
        let w = a.iter().map(|h| h[3]).collect::<Vec<_>>();
        let inverse_weight = reciprocal(hull(w.iter().copied())?)?;
        let inverse_squared = inverse_weight.sqr();
        let profile_scale =
            global_width * reciprocal(Interval::point(span[1]) - Interval::point(span[0]))?;
        let bounds = |axis| -> Result<(Interval, Interval)> {
            let ha = project(a, axis);
            let hb = project(b, axis);
            let qa = derivative_numerator(&ha, &w)?;
            let qb = derivative_numerator(&hb, &w)?;
            let profile = hull(qa.into_iter().chain(qb))? * inverse_squared * profile_scale;
            let ruling = hull(ha.into_iter().zip(hb).map(|(a, b)| b - a))? * inverse_weight;
            Ok((profile, ruling))
        };
        let (profile_x, ruling_x) = bounds(frame.x())?;
        let (profile_y, ruling_y) = bounds(frame.y())?;
        let (j00, j11, j01, j10) = if linear_u {
            (ruling_x, profile_y, profile_x, ruling_y)
        } else {
            (profile_x, ruling_y, ruling_x, profile_y)
        };
        let off = (j01 + j10) * Interval::point(0.5);
        let margin = Interval::point(p.linear());
        if (j00 - margin).sign() != Some(Sign::Positive)
            || (j11 - margin).sign() != Some(Sign::Positive)
            || (j00 * j11 - off.sqr() - margin.sqr()).sign() != Some(Sign::Positive)
        {
            return Err(refused());
        }
    }
    Ok(())
}

/// Prove the unchanged global graph inequalities in the caller's single frame.
/// The caller keeps its size and continuity guards and general fallback budget.
pub(crate) fn certify(s: &BSurface, frame: Frame, p: &Precision) -> Result<()> {
    certify_direction(s, frame, false, p).or_else(|_| certify_direction(s, frame, true, p))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn ruled(height: f64, unequal_row_weights: bool) -> BSurface {
        let mut weights = vec![1.0, 1.0, 0.8, 0.8, 1.0, 1.0];
        if unequal_row_weights {
            weights[3] = 1.6;
        }
        BSurface::new(
            2,
            1,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            3,
            2,
            [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.1),
                Point3::new(2.0, 0.0, 0.0),
            ]
            .into_iter()
            .flat_map(|point| [point, point + Vec3::Y * height])
            .collect(),
            Some(weights),
        )
        .unwrap()
    }

    #[test]
    fn original_jacobian_margins_reject_singular_and_near_degenerate_rulings() {
        let p = Precision::DEFAULT;
        assert!(certify(&ruled(1.0, false), Frame::WORLD, &p).is_ok());
        assert!(certify(&ruled(0.0, false), Frame::WORLD, &p).is_err());
        assert!(certify(&ruled(p.linear() * 0.5, false), Frame::WORLD, &p).is_err());
    }

    #[test]
    fn unequal_row_weights_do_not_enter_the_affine_ruling_proof() {
        // General admission may use its unchanged fallback for this valid
        // NURBS family, but this shared-denominator proof must refuse it.
        assert!(certify(&ruled(1.0, true), Frame::WORLD, &Precision::DEFAULT).is_err());
    }

    #[test]
    fn shared_raw_extraction_rejects_degree_and_span_work_limits() {
        let high_degree = BCurve::new(
            18,
            [vec![0.0; 19], vec![1.0; 19]].concat(),
            (0..19)
                .map(|i| Point3::new(f64::from(i), 0.0, 0.0))
                .collect(),
            None,
        )
        .unwrap();
        assert!(homogeneous_spans(&high_degree, Point3::ORIGIN).is_err());

        let many_spans = BCurve::new(
            1,
            std::iter::once(0.0)
                .chain((0..=257).map(f64::from))
                .chain(std::iter::once(257.0))
                .collect(),
            (0..=257)
                .map(|i| Point3::new(f64::from(i), 0.0, 0.0))
                .collect(),
            None,
        )
        .unwrap();
        assert!(homogeneous_spans(&many_spans, Point3::ORIGIN).is_err());
    }
}
