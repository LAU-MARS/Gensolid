//! General NURBS surface intersection exploration in parameter space.
//! Interval exclusions preserve a cover of every possible intersection; continuation
//! supplies residual-checked numerical curves. Curve connectivity near singularities
//! is not certified and is never returned as an exact IntersectionSet.
use crate::nurbs_bounds::{bound_surface, knot_spans};
use crate::nurbs_intersect::{IntersectionOptions, inverse};
use crate::{BSurface, ParamRange};
use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Point3, Vec3};

/// Residual-checked numerical samples on both original surfaces.
#[derive(Clone, Debug)]
pub struct IntersectionSample {
    /// U,V on the first surface followed by U,V on the second.
    pub parameters: [f64; 4],
    /// Point evaluated on the first surface.
    pub point: Point3,
    /// Distance to the second surface evaluated at its paired parameters.
    pub residual: f64,
}
/// Numerical continuation branch; connectivity is a candidate, not a B-rep edge.
#[derive(Clone, Debug)]
pub struct NumericalIntersectionCurve {
    /// Ordered samples. Display interpolation is not used by kernel operations.
    pub samples: Vec<IntersectionSample>,
}
/// Honest surface/surface result: global enclosure plus numerical candidate curves.
#[derive(Clone, Debug)]
pub struct NurbsSurfaceIntersectionReport {
    /// Direct parameter-space numerical continuation curves.
    pub curves: Vec<NumericalIntersectionCurve>,
    /// Parameter boxes covering ALL roots not excluded by outward-rounded bounds.
    /// These may contain false positives; they do not assert branch completeness.
    pub possible_cells: Vec<[ParamRange; 4]>,
    /// Number of visited interval cells.
    pub visited_cells: usize,
    /// Whether the subdivision budget was exhausted.
    pub budget_exhausted: bool,
    /// True when the two original NURBS representations are identical.
    pub identical_overlap: bool,
}
impl NurbsSurfaceIntersectionReport {
    /// A non-overlapping empty cover is a certified empty intersection.
    pub fn certified_empty(&self) -> bool {
        !self.identical_overlap && self.possible_cells.is_empty()
    }
}
fn c(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn sample(a: &BSurface, b: &BSurface, x: [f64; 4]) -> IntersectionSample {
    let p = a.eval(x[0], x[1]).point;
    let q = b.eval(x[2], x[3]).point;
    IntersectionSample {
        parameters: x,
        point: p,
        residual: p.distance(q),
    }
}
fn system(a: &BSurface, b: &BSurface, x: [f64; 4]) -> ([f64; 3], [[f64; 4]; 3]) {
    let p = a.eval(x[0], x[1]);
    let q = b.eval(x[2], x[3]);
    let columns = [c(p.du), c(p.dv), c(-q.du), c(-q.dv)];
    (
        c(p.point - q.point),
        std::array::from_fn(|i| std::array::from_fn(|j| columns[j][i])),
    )
}
fn domain(a: &BSurface, b: &BSurface) -> [ParamRange; 4] {
    let (au, av) = a.ranges();
    let (bu, bv) = b.ranges();
    [au, av, bu, bv]
}
fn correct(
    a: &BSurface,
    b: &BSurface,
    mut x: [f64; 4],
    bounds: [ParamRange; 4],
    p: &Precision,
) -> Option<IntersectionSample> {
    for _ in 0..32 {
        let point = sample(a, b, x);
        if point.residual <= p.solver() {
            return Some(point);
        }
        let (f, j) = system(a, b, x);
        let gram = std::array::from_fn(|i| {
            std::array::from_fn(|k| j[i].iter().zip(j[k]).map(|(a, b)| a * b).sum())
        });
        let inv = inverse(gram)?;
        let lambda: [f64; 3] =
            std::array::from_fn(|i| inv[i].iter().zip(f).map(|(a, b)| a * b).sum());
        let delta: [f64; 4] = std::array::from_fn(|i| (0..3).map(|k| j[k][i] * lambda[k]).sum());
        let mut accepted = None;
        for k in 0..12 {
            let alpha = 0.5_f64.powi(k);
            let next = std::array::from_fn(|i| {
                (x[i] - alpha * delta[i]).clamp(bounds[i].lo, bounds[i].hi)
            });
            if next.iter().all(|v| v.is_finite()) && sample(a, b, next).residual < point.residual {
                accepted = Some(next);
                break;
            }
        }
        x = accepted?;
    }
    None
}
fn determinant(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}
fn tangent(
    a: &BSurface,
    b: &BSurface,
    x: [f64; 4],
    d: [ParamRange; 4],
    p: &Precision,
) -> Option<[f64; 4]> {
    let (_, j) = system(a, b, x);
    let mut t = std::array::from_fn(|skip| {
        let indices = (0..4).filter(|&i| i != skip).collect::<Vec<_>>();
        let m = std::array::from_fn(|i| {
            std::array::from_fn(|k| j[i][indices[k]] * (d[indices[k]].hi - d[indices[k]].lo))
        });
        determinant(m) * if skip % 2 == 0 { 1.0 } else { -1.0 }
    });
    let scale = t.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !scale.is_finite()
        || scale
            <= p.angular()
                * j.iter()
                    .map(|row| row.iter().map(|v| v * v).sum::<f64>().sqrt())
                    .product::<f64>()
    {
        return None;
    }
    for v in &mut t {
        *v /= scale;
    }
    Some(t)
}
fn distance(a: [f64; 4], b: [f64; 4], d: [ParamRange; 4]) -> f64 {
    (0..4)
        .map(|i| ((a[i] - b[i]) / (d[i].hi - d[i].lo)).powi(2))
        .sum::<f64>()
        .sqrt()
}
fn trace(
    a: &BSurface,
    b: &BSurface,
    seed: &IntersectionSample,
    sign: f64,
    p: &Precision,
) -> Vec<IntersectionSample> {
    let d = domain(a, b);
    let mut current = seed.clone();
    let mut old = None;
    let mut out = Vec::new();
    // Fixed continuation sampling resolution is a display/work budget, not a
    // topological predicate. Every sample still satisfies the original residual.
    let step: f64 = 1.0 / 64.0;
    for _ in 0..1024 {
        let Some(mut direction) = tangent(a, b, current.parameters, d, p) else {
            break;
        };
        if let Some(previous) = old {
            if direction
                .iter()
                .zip(previous)
                .map(|(a, b)| a * b)
                .sum::<f64>()
                < 0.0
            {
                direction = direction.map(|x| -x);
            }
        } else {
            direction = direction.map(|x| x * sign);
        }
        let mut length = step;
        for i in 0..4 {
            let delta = direction[i] * (d[i].hi - d[i].lo);
            if delta > 0.0 {
                length = length.min((d[i].hi - current.parameters[i]) / delta);
            } else if delta < 0.0 {
                length = length.min((d[i].lo - current.parameters[i]) / delta);
            }
        }
        if length <= p.angular() {
            break;
        }
        let target = std::array::from_fn(|i| {
            current.parameters[i] + length * direction[i] * (d[i].hi - d[i].lo)
        });
        let Some(next) = correct(a, b, target, d, p) else {
            break;
        };
        if distance(current.parameters, next.parameters, d) <= p.angular() {
            break;
        }
        let closes = out.len() > 8 && distance(next.parameters, seed.parameters, d) < step * 0.5;
        old = Some(direction);
        current = next;
        out.push(current.clone());
        if closes {
            break;
        }
    }
    out
}
/// Intersect arbitrary-degree positive-weight NURBS surfaces directly in their
/// bounded knot domains. This provides certified exclusion/overlap plus numerical
/// curves and an explicit global possible-root cover. It never asserts certified
/// branch connectivity for nonempty regular, singular or tangential intersections.
pub fn intersect_nurbs_surfaces(
    a: &BSurface,
    b: &BSurface,
    p: &Precision,
    options: IntersectionOptions,
) -> Result<NurbsSurfaceIntersectionReport> {
    if options.max_cells == 0 || options.max_depth == 0 {
        return Err(Error::invalid_argument(
            "surface intersection budgets must be nonzero",
        ));
    }
    let d = domain(a, b);
    let mut result = NurbsSurfaceIntersectionReport {
        curves: Vec::new(),
        possible_cells: Vec::new(),
        visited_cells: 0,
        budget_exhausted: false,
        identical_overlap: a == b,
    };
    if result.identical_overlap {
        return Ok(result);
    }
    let (pu, pv) = a.degrees();
    let (qu, qv) = b.degrees();
    let (n, m) = a.counts();
    let (k, l) = b.counts();
    let (au, av) = a.knots();
    let (bu, bv) = b.knots();
    let spans = [
        knot_spans(au, pu, n),
        knot_spans(av, pv, m),
        knot_spans(bu, qu, k),
        knot_spans(bv, qv, l),
    ];
    let count = spans.iter().try_fold(1usize, |n, s| n.checked_mul(s.len()));
    if count.is_none_or(|n| n > options.max_cells) {
        result.possible_cells.push(d);
        result.budget_exhausted = true;
        return Ok(result);
    }
    let mut pending = Vec::new();
    for &u in &spans[0] {
        for &v in &spans[1] {
            for &s in &spans[2] {
                for &t in &spans[3] {
                    pending.push(([u, v, s, t], 0usize));
                }
            }
        }
    }
    while let Some((cell, depth)) = pending.pop() {
        if result.visited_cells >= options.max_cells {
            result.budget_exhausted = true;
            result.possible_cells.push(cell);
            result
                .possible_cells
                .extend(pending.into_iter().map(|x| x.0));
            break;
        }
        result.visited_cells += 1;
        let aa = bound_surface(a, cell[0], cell[1])?;
        let bb = bound_surface(b, cell[2], cell[3])?;
        if (0..3).any(|i| !(aa.point[i] - bb.point[i]).contains(0.0)) {
            continue;
        }
        if depth >= options.max_depth.min(16) {
            result.possible_cells.push(cell);
            continue;
        }
        let axis = (0..4)
            .max_by(|&i, &j| {
                ((cell[i].hi - cell[i].lo) / (d[i].hi - d[i].lo))
                    .total_cmp(&((cell[j].hi - cell[j].lo) / (d[j].hi - d[j].lo)))
            })
            .ok_or_else(|| Error::internal("no parameter axis"))?;
        let mid = cell[axis].lo + (cell[axis].hi - cell[axis].lo) * 0.5;
        if mid <= cell[axis].lo || mid >= cell[axis].hi {
            result.possible_cells.push(cell);
            continue;
        }
        let mut left = cell;
        let mut right = cell;
        left[axis].hi = mid;
        right[axis].lo = mid;
        pending.push((right, depth + 1));
        pending.push((left, depth + 1));
    }
    // Numerical candidate association is only for display, never exact B-rep topology.
    for cell in &result.possible_cells {
        let center = cell.map(|r| r.lo + (r.hi - r.lo) * 0.5);
        if result.curves.iter().any(|curve| {
            curve
                .samples
                .iter()
                .any(|s| distance(center, s.parameters, d) < 1.0 / 32.0)
        }) {
            continue;
        }
        let Some(seed) = correct(a, b, center, d, p) else {
            continue;
        };
        if result.curves.iter().any(|curve| {
            curve
                .samples
                .iter()
                .any(|s| distance(seed.parameters, s.parameters, d) < 1.0 / 32.0)
        }) {
            continue;
        }
        let mut samples = trace(a, b, &seed, -1.0, p);
        samples.reverse();
        samples.push(seed.clone());
        samples.extend(trace(a, b, &seed, 1.0, p));
        result.curves.push(NumericalIntersectionCurve { samples });
        if result.curves.len() >= 128 {
            break;
        }
    }
    Ok(result)
}
#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    fn patch(vertical: bool, z: f64) -> BSurface {
        let points = (0..3)
            .flat_map(|i| {
                (0..3).map(move |j| {
                    let u = i as f64 / 2.0;
                    let v = j as f64 / 2.0;
                    if vertical {
                        Point3::new(u, 0.45, v - 0.5 + z)
                    } else {
                        Point3::new(u, v, z)
                    }
                })
            })
            .collect();
        BSurface::new(
            2,
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            3,
            3,
            points,
            None,
        )
        .unwrap()
    }
    #[test]
    fn general_surface_trace_preserves_original_residual_and_global_cover() {
        let a = patch(false, 0.0);
        let b = patch(true, 0.0);
        let p = Precision::DEFAULT;
        let r = intersect_nurbs_surfaces(
            &a,
            &b,
            &p,
            IntersectionOptions {
                max_cells: 5000,
                max_depth: 12,
            },
        )
        .unwrap();
        assert!(!r.certified_empty());
        assert!(!r.curves.is_empty());
        assert!(!r.possible_cells.is_empty());
        for curve in r.curves {
            assert!(curve.samples.len() > 2);
            for sample in curve.samples {
                assert!(sample.residual <= p.solver());
                assert!((sample.point.y - 0.45).abs() <= p.linear());
            }
        }
    }
    #[test]
    fn disjoint_and_identical_surfaces_have_definitive_classifications() {
        let a = patch(false, 0.0);
        let b = patch(false, 2.0);
        let p = Precision::DEFAULT;
        assert!(
            intersect_nurbs_surfaces(&a, &b, &p, Default::default())
                .unwrap()
                .certified_empty()
        );
        let r = intersect_nurbs_surfaces(&a, &a, &p, Default::default()).unwrap();
        assert!(r.identical_overlap);
        assert!(!r.certified_empty());
    }
}
