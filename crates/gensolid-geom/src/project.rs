//! Numerical closest-point search for curves and surfaces without closed forms.

use gensolid_core::{Error, Precision, Result};
use gensolid_math::Point3;

use crate::curve::CurveDerivs;
use crate::param::ParamRange;
use crate::store::CurveProjection;
use crate::surface::SurfaceDerivs;

const CURVE_SAMPLES: usize = 64;
const SURFACE_SAMPLES: usize = 16;
const MAX_ITERATIONS: usize = 64;
const MAX_HALVINGS: usize = 40;

fn require_bounded(range: ParamRange) -> Result<()> {
    if range.is_bounded() {
        Ok(())
    } else {
        Err(Error::unsupported(
            "numerical projection needs a bounded parameter range",
        ))
    }
}

pub(crate) fn closest_on_curve(
    eval: impl Fn(f64) -> CurveDerivs,
    range: ParamRange,
    p: Point3,
    precision: &Precision,
) -> Result<CurveProjection> {
    require_bounded(range)?;
    let at = |t: f64| {
        let point = eval(t).point;
        CurveProjection {
            param: t,
            point,
            distance: point.distance(p),
        }
    };
    let step = (range.hi - range.lo) / CURVE_SAMPLES as f64;
    let mut best = (0..=CURVE_SAMPLES)
        .map(|i| at(range.lo + step * i as f64))
        .min_by(|a, b| a.distance.total_cmp(&b.distance))
        .ok_or_else(|| Error::internal("no samples"))?;

    let tiny = precision.degenerate_length() * precision.degenerate_length();
    for _ in 0..MAX_ITERATIONS {
        let d = eval(best.param);
        let speed2 = d.d1.length_squared();
        if speed2 <= tiny {
            break;
        }
        let step = (p - d.point).dot(d.d1) / speed2;
        let accepted = descend(
            |scale| at(range.normalize(best.param + step * scale)),
            best.distance,
        );
        let Some(next) = accepted else { break };
        let moved = (next.param - best.param).abs() * speed2.sqrt();
        best = next;
        if moved <= precision.solver() {
            break;
        }
    }
    Ok(best)
}

/// Backtracking line search: the first candidate at scale 1, 1/2, 1/4, ... that is closer
/// than `current`. Plain Gauss-Newton overshoots on strongly curved geometry far from the
/// query point, so a full step is never accepted blindly.
fn descend<T>(candidate: impl Fn(f64) -> T, current: f64) -> Option<T>
where
    T: Distance,
{
    let mut scale = 1.0;
    for _ in 0..MAX_HALVINGS {
        let c = candidate(scale);
        if c.distance() < current {
            return Some(c);
        }
        scale *= 0.5;
    }
    None
}

trait Distance {
    fn distance(&self) -> f64;
}

impl Distance for CurveProjection {
    fn distance(&self) -> f64 {
        self.distance
    }
}

impl Distance for (f64, f64, f64) {
    fn distance(&self) -> f64 {
        self.0
    }
}

pub(crate) fn distance_to_patch(
    eval: impl Fn(f64, f64) -> SurfaceDerivs,
    ur: ParamRange,
    vr: ParamRange,
    p: Point3,
    precision: &Precision,
) -> Result<f64> {
    require_bounded(ur)?;
    require_bounded(vr)?;
    let (du, dv) = (
        (ur.hi - ur.lo) / SURFACE_SAMPLES as f64,
        (vr.hi - vr.lo) / SURFACE_SAMPLES as f64,
    );
    let mut best = (f64::INFINITY, ur.lo, vr.lo);
    for i in 0..=SURFACE_SAMPLES {
        for j in 0..=SURFACE_SAMPLES {
            let (u, v) = (ur.lo + du * i as f64, vr.lo + dv * j as f64);
            let dist = eval(u, v).point.distance(p);
            if dist < best.0 {
                best = (dist, u, v);
            }
        }
    }

    for _ in 0..MAX_ITERATIONS {
        let (dist, u, v) = best;
        let d = eval(u, v);
        let r = p - d.point;
        let (a, b, c) = (d.du.dot(d.du), d.du.dot(d.dv), d.dv.dot(d.dv));
        let det = a * c - b * b;
        if !(det.is_finite() && det.abs() > precision.angular() * a * c) {
            break;
        }
        let (gu, gv) = (r.dot(d.du), r.dot(d.dv));
        let (step_u, step_v) = ((c * gu - b * gv) / det, (a * gv - b * gu) / det);
        let candidate = |scale: f64| {
            let (nu, nv) = (
                ur.normalize(u + step_u * scale),
                vr.normalize(v + step_v * scale),
            );
            (eval(nu, nv).point.distance(p), nu, nv)
        };
        let Some(next) = descend(candidate, dist) else {
            break;
        };
        let moved = (d.du * (next.1 - u) + d.dv * (next.2 - v)).length();
        best = next;
        if moved <= precision.solver() {
            break;
        }
    }
    Ok(best.0)
}
