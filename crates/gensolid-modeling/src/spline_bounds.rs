//! Outward-rounded homogeneous extraction from raw stored NURBS inputs.
//! Original Boehm knot insertion; no floating trimmed curve is used as proof.
use gensolid_core::{Error, Result};
use gensolid_geom::BCurve;
use gensolid_math::{Interval, Point3, Sign};

fn refused(message: &str) -> Error {
    Error::unimplemented(message)
}
fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}
fn coordinates(point: Point3) -> [f64; 3] {
    [point.x, point.y, point.z]
}
fn weight(curve: &BCurve, i: usize) -> f64 {
    curve.weights().map_or(1.0, |w| w[i])
}
fn reciprocal(value: Interval) -> Result<Interval> {
    if value.sign() != Some(Sign::Positive) {
        return Err(refused("homogeneous denominator is not positive"));
    }
    Interval::new((1.0 / value.hi()).next_down(), (1.0 / value.lo()).next_up())
        .filter(|v| v.lo().is_finite() && v.hi().is_finite())
        .ok_or_else(|| refused("homogeneous denominator overflow"))
}
fn positive_min(values: &[Interval]) -> Result<f64> {
    if values
        .iter()
        .any(|v| v.sign() != Some(Sign::Positive) || !v.lo().is_finite() || !v.hi().is_finite())
    {
        return Err(refused("homogeneous denominator is unresolved"));
    }
    Ok(values.iter().map(|v| v.lo()).fold(f64::INFINITY, f64::min))
}
/// Homogeneous Bernstein coefficients enclosed directly from the stored basis.
/// A single global weight scale cancels in the quotient and keeps common scales
/// such as 2^+/-400 harmless. The division itself is outward-rounded.
pub(crate) fn homogeneous_spans(
    curve: &BCurve,
    reference: Point3,
) -> Result<Vec<Vec<[Interval; 4]>>> {
    let p = curve.degree();
    let range = curve.range();
    if p > 17
        || curve.knots()[..=p].iter().any(|&k| !same(k, range.lo))
        || curve.knots()[curve.points().len()..]
            .iter()
            .any(|&k| !same(k, range.hi))
    {
        return Err(refused("homogeneous extraction needs clamped degree <= 17"));
    }
    let scale = (0..curve.points().len())
        .map(|i| weight(curve, i))
        .fold(0.0, f64::max);
    let inv_scale = reciprocal(Interval::point(scale))?;
    let mut controls: Vec<[Interval; 4]> = curve
        .points()
        .iter()
        .enumerate()
        .map(|(i, &point)| {
            let w = Interval::point(weight(curve, i)) * inv_scale;
            std::array::from_fn(|j| {
                if j == 3 {
                    w
                } else {
                    (Interval::point(coordinates(point)[j])
                        - Interval::point(coordinates(reference)[j]))
                        * w
                }
            })
        })
        .collect();
    let mut knots = curve.knots().to_vec();
    let r = curve.range();
    let mut unique = Vec::new();
    for &k in &knots[p + 1..curve.points().len()] {
        if !same(k, r.lo) && !same(k, r.hi) && unique.last().is_none_or(|&v| !same(k, v)) {
            unique.push(k);
        }
    }
    if unique.len() >= 256 {
        return Err(refused("homogeneous extraction exceeds 256 spans"));
    }
    for knot in unique {
        let mut multiplicity = knots.iter().filter(|&&k| same(k, knot)).count();
        if multiplicity > p {
            return Err(refused("homogeneous extraction found a full break"));
        }
        while multiplicity < p {
            // Knot ordering is representation identity, not a topology decision.
            let k = knots.partition_point(|&v| v <= knot) - 1;
            let mut next = Vec::with_capacity(controls.len() + 1);
            next.extend_from_slice(&controls[..=k - p]);
            for i in k - p + 1..=k - multiplicity {
                let denominator = Interval::point(knots[i + p]) - Interval::point(knots[i]);
                if denominator.sign() != Some(Sign::Positive) {
                    return Err(refused("homogeneous knot interval is unresolved"));
                }
                let alpha =
                    (Interval::point(knot) - Interval::point(knots[i])) * reciprocal(denominator)?;
                next.push(std::array::from_fn(|j| {
                    (Interval::point(1.0) - alpha) * controls[i - 1][j] + alpha * controls[i][j]
                }));
            }
            next.extend_from_slice(&controls[k - multiplicity..]);
            knots.insert(k + 1, knot);
            controls = next;
            multiplicity += 1;
        }
    }
    let mut out = Vec::new();
    for k in p..controls.len() {
        if knots[k] < knots[k + 1] {
            let span = controls[k - p..=k].to_vec();
            positive_min(&span.iter().map(|q| q[3]).collect::<Vec<_>>())?;
            if span
                .iter()
                .flatten()
                .any(|v| !v.lo().is_finite() || !v.hi().is_finite())
            {
                return Err(refused("homogeneous extraction overflow"));
            }
            out.push(span);
        }
    }
    Ok(out)
}
