//! Bounded numerical point projection.
//!
//! These searches produce converged candidates, not certificates of global optimality.
//! Comparisons in this module control numerical minimisation, never B-rep topology.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::Point3;

use crate::curve::CurveDerivs;
use crate::param::ParamRange;
use crate::store::{CurveProjection, ProjectionKind, SurfaceProjection};
use crate::surface::SurfaceDerivs;

const CURVE_SAMPLES: usize = 64;
const SURFACE_SAMPLES: usize = 16;
const MAX_ITERATIONS: usize = 96;
const MAX_HALVINGS: usize = 40;

fn require_bounded(range: ParamRange) -> Result<()> {
    if !range.is_bounded() {
        return Err(Error::unimplemented(
            "numerical projection needs a bounded parameter range",
        ));
    }
    if range.lo >= range.hi || !(range.hi - range.lo).is_finite() {
        return Err(Error::invalid_argument(
            "projection range must have finite positive width",
        ));
    }
    Ok(())
}

fn curve_sample(
    eval: &impl Fn(f64) -> CurveDerivs,
    t: f64,
    p: Point3,
) -> Result<(CurveProjection, f64)> {
    let value = eval(t);
    let distance = value.point.distance(p);
    let speed = value.d1.length();
    if !(value.point.is_finite() && value.d1.is_finite() && distance.is_finite()) {
        return Err(Error::degenerate("non-finite numerical curve evaluation"));
    }
    let gradient = if speed > 0.0 {
        (value.point - p).dot(value.d1 / speed)
    } else {
        0.0
    };
    if !gradient.is_finite() {
        return Err(Error::degenerate("non-finite curve projection gradient"));
    }
    Ok((
        CurveProjection {
            param: t,
            point: value.point,
            distance,
            kind: ProjectionKind::NumericalCandidate,
        },
        gradient,
    ))
}

fn keep_curve(best: &mut Option<CurveProjection>, candidate: CurveProjection) {
    if best.is_none_or(|current| candidate.distance < current.distance) {
        *best = Some(candidate);
    }
}

pub(crate) fn closest_on_curve(
    eval: impl Fn(f64) -> CurveDerivs,
    range: ParamRange,
    p: Point3,
    precision: &Precision,
) -> Result<CurveProjection> {
    require_bounded(range)?;
    let mut samples = Vec::with_capacity(CURVE_SAMPLES + 1);
    let mut best = None;
    let mut sampled_distance = f64::INFINITY;
    for i in 0..=CURVE_SAMPLES {
        let t = range.lo + (range.hi - range.lo) * (i as f64 / CURVE_SAMPLES as f64);
        let sample = curve_sample(&eval, t, p)?;
        sampled_distance = sampled_distance.min(sample.0.distance);
        if sample.1.abs() <= precision.solver()
            || (!range.periodic && (i == 0 || i == CURVE_SAMPLES))
        {
            keep_curve(&mut best, sample.0);
        }
        samples.push(sample);
    }
    // A negative-to-positive derivative transition brackets a stationary minimum.
    // Bisection cannot overshoot; its stopping criterion is the spatial tangential
    // residual rather than a parameter tolerance dependent on arbitrary knot scales.
    for pair in samples.windows(2) {
        if pair[0].1 >= 0.0 || pair[1].1 <= 0.0 {
            continue;
        }
        let (mut lo, mut hi) = (pair[0].0.param, pair[1].0.param);
        for _ in 0..MAX_ITERATIONS {
            let mid = lo + (hi - lo) * 0.5;
            let sample = curve_sample(&eval, mid, p)?;
            if sample.1.abs() <= precision.solver() || sample.0.distance <= precision.solver() {
                keep_curve(&mut best, sample.0);
                break;
            }
            if mid <= lo || mid >= hi {
                break;
            }
            if sample.1 < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
    }
    let best = best.ok_or_else(|| Error::unimplemented("curve projection did not converge"))?;
    if best.distance > sampled_distance + precision.linear() {
        return Err(Error::unimplemented(
            "curve projection did not converge near its best seed",
        ));
    }
    Ok(best)
}

/// Search every nonempty knot span, so narrow knot intervals cannot disappear
/// between uniform seeds. This still does not certify all minima inside a span.
pub(crate) fn closest_on_spline(
    eval: impl Fn(f64) -> CurveDerivs,
    range: ParamRange,
    knots: &[f64],
    p: Point3,
    precision: &Precision,
) -> Result<CurveProjection> {
    let mut best = None;
    for pair in knots.windows(2) {
        let lo = pair[0].max(range.lo);
        let hi = pair[1].min(range.hi);
        if lo >= hi {
            continue;
        }
        // At a fully repeated interior knot the two one-sided values can differ.
        // The upper domain endpoint is evaluated exactly; interior ends use the
        // last representable parameter on the left and the next span covers right.
        let end = if hi < range.hi { hi.next_down() } else { hi };
        if end <= lo {
            return Err(Error::unimplemented(
                "knot span has no representable interior",
            ));
        }
        keep_curve(
            &mut best,
            closest_on_curve(&eval, ParamRange::bounded(lo, end), p, precision)?,
        );
    }
    best.ok_or_else(|| Error::degenerate("spline has no nonempty knot spans"))
}

fn patch_sample(
    eval: &impl Fn(f64, f64) -> SurfaceDerivs,
    u: f64,
    v: f64,
    p: Point3,
) -> Result<(SurfaceProjection, SurfaceDerivs)> {
    let value = eval(u, v);
    let distance = value.point.distance(p);
    if !(value.point.is_finite()
        && value.du.is_finite()
        && value.dv.is_finite()
        && distance.is_finite())
    {
        return Err(Error::degenerate("non-finite numerical surface evaluation"));
    }
    Ok((
        SurfaceProjection {
            u,
            v,
            point: value.point,
            distance,
            kind: ProjectionKind::NumericalCandidate,
        },
        value,
    ))
}

fn blocked(range: ParamRange, value: f64, gradient: f64) -> bool {
    !range.periodic
        && ((value <= range.lo && gradient < 0.0) || (value >= range.hi && gradient > 0.0))
}

fn refine_patch(
    eval: &impl Fn(f64, f64) -> SurfaceDerivs,
    ranges: (ParamRange, ParamRange),
    seed: SurfaceProjection,
    p: Point3,
    precision: &Precision,
) -> Result<Option<SurfaceProjection>> {
    let (ur, vr) = ranges;
    let mut current = seed;
    for _ in 0..MAX_ITERATIONS {
        let (_, d) = patch_sample(eval, current.u, current.v, p)?;
        let r = p - d.point;
        let (su, sv) = (d.du.length(), d.dv.length());
        let tu = d.du.normalized(0.0);
        let tv = d.dv.normalized(0.0);
        let gu = tu.map_or(0.0, |t| r.dot(t));
        let gv = tv.map_or(0.0, |t| r.dot(t));
        let active_u = tu.is_some() && !blocked(ur, current.u, gu);
        let active_v = tv.is_some() && !blocked(vr, current.v, gv);
        let residual =
            (if active_u { gu.abs() } else { 0.0 }).max(if active_v { gv.abs() } else { 0.0 });
        if current.distance <= precision.solver() || residual <= precision.solver() {
            return Ok(Some(current));
        }
        let (mut step_u, mut step_v) = (
            if active_u { gu / su } else { 0.0 },
            if active_v { gv / sv } else { 0.0 },
        );
        if let (Some(tu), Some(tv)) = (tu, tv) {
            let cosine = tu.dot(tv);
            let determinant = 1.0 - cosine * cosine;
            if active_u && active_v && determinant > precision.angular() {
                step_u = (gu - cosine * gv) / (determinant * su);
                step_v = (gv - cosine * gu) / (determinant * sv);
            }
        }
        if !(step_u.is_finite() && step_v.is_finite()) {
            return Ok(None);
        }
        let mut accepted = None;
        let mut scale = 1.0;
        for _ in 0..MAX_HALVINGS {
            let u = ur.normalize(current.u + step_u * scale);
            let v = vr.normalize(current.v + step_v * scale);
            let (candidate, next) = patch_sample(eval, u, v, p)?;
            let next_r = p - next.point;
            let next_gu = next.du.normalized(0.0).map_or(0.0, |t| next_r.dot(t));
            let next_gv = next.dv.normalized(0.0).map_or(0.0, |t| next_r.dot(t));
            let next_residual = (if blocked(ur, u, next_gu) {
                0.0
            } else {
                next_gu.abs()
            })
            .max(if blocked(vr, v, next_gv) {
                0.0
            } else {
                next_gv.abs()
            });
            // Near the minimum, the quadratic distance improvement can round away
            // while the tangential residual still improves. Permit a solver-sized
            // distance plateau only together with better first-order optimality.
            if candidate.distance < current.distance
                || (candidate.distance <= current.distance + precision.solver()
                    && next_residual < residual)
            {
                accepted = Some(candidate);
                break;
            }
            scale *= 0.5;
        }
        let Some(next) = accepted else {
            return Ok(None);
        };
        current = next;
    }
    Ok(None)
}

pub(crate) fn closest_on_patch(
    eval: impl Fn(f64, f64) -> SurfaceDerivs,
    ur: ParamRange,
    vr: ParamRange,
    p: Point3,
    precision: &Precision,
) -> Result<SurfaceProjection> {
    require_bounded(ur)?;
    require_bounded(vr)?;
    let count = SURFACE_SAMPLES + 1;
    let mut samples = Vec::with_capacity(count * count);
    let mut sampled_distance = f64::INFINITY;
    for i in 0..count {
        for j in 0..count {
            let u = ur.lo + (ur.hi - ur.lo) * (i as f64 / SURFACE_SAMPLES as f64);
            let v = vr.lo + (vr.hi - vr.lo) * (j as f64 / SURFACE_SAMPLES as f64);
            let (sample, _) = patch_sample(&eval, u, v, p)?;
            sampled_distance = sampled_distance.min(sample.distance);
            samples.push(sample);
        }
    }
    let mut best: Option<SurfaceProjection> = None;
    for i in 0..count {
        for j in 0..count {
            let seed = samples[i * count + j];
            let local_minimum = (i.saturating_sub(1)..=(i + 1).min(SURFACE_SAMPLES)).all(|ni| {
                (j.saturating_sub(1)..=(j + 1).min(SURFACE_SAMPLES))
                    .all(|nj| samples[ni * count + nj].distance >= seed.distance)
            });
            if local_minimum
                && let Some(candidate) = refine_patch(&eval, (ur, vr), seed, p, precision)?
                && best.is_none_or(|current| candidate.distance < current.distance)
            {
                best = Some(candidate);
            }
        }
    }
    let best = best.ok_or_else(|| Error::unimplemented("surface projection did not converge"))?;
    if best.distance > sampled_distance + precision.linear() {
        return Err(Error::unimplemented(
            "surface projection did not converge near its best seed",
        ));
    }
    Ok(best)
}

pub(crate) fn closest_on_spline_patch(
    eval: impl Fn(f64, f64) -> SurfaceDerivs,
    ranges: (ParamRange, ParamRange),
    knots: (&[f64], &[f64]),
    p: Point3,
    precision: &Precision,
) -> Result<SurfaceProjection> {
    let span = |pair: &[f64], range: ParamRange| {
        let lo = pair[0].max(range.lo);
        let hi = pair[1].min(range.hi);
        (lo < hi).then(|| ParamRange::bounded(lo, if hi < range.hi { hi.next_down() } else { hi }))
    };
    let mut best: Option<SurfaceProjection> = None;
    for ur in knots.0.windows(2).filter_map(|pair| span(pair, ranges.0)) {
        for vr in knots.1.windows(2).filter_map(|pair| span(pair, ranges.1)) {
            let candidate = closest_on_patch(&eval, ur, vr, p, precision)?;
            if best.is_none_or(|current| candidate.distance < current.distance) {
                best = Some(candidate);
            }
        }
    }
    best.ok_or_else(|| Error::degenerate("surface has no nonempty knot spans"))
}
