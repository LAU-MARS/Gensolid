//! Homogeneous global interpolation of compatible NURBS section curves.
use crate::{BCurve, BSurface};
use gensolid_core::{Error, Precision, Result};
use gensolid_math::Point3;

/// Interpolate two or more compatible positive-weight NURBS section curves.
/// Curves share degree, knot vector and control count (explicit correspondence).
/// The result passes through every section at its supplied parameter. Degree one
/// is ruled; higher degrees produce a smooth interpolating surface. Invalid weights
/// or ill-conditioned interpolation are rejected rather than approximated silently.
pub fn loft_nurbs_sections(
    curves: &[BCurve],
    parameters: &[f64],
    degree: usize,
    p: &Precision,
) -> Result<BSurface> {
    let n = curves.len();
    if n < 2 || n != parameters.len() || degree == 0 || degree >= n || n > 64 {
        return Err(Error::invalid_argument(
            "loft requires 2..64 sections and degree in 1..section_count",
        ));
    }
    let first = &curves[0];
    if curves.iter().any(|c| {
        c.degree() != first.degree()
            || c.knots() != first.knots()
            || c.points().len() != first.points().len()
    }) {
        return Err(Error::invalid_argument(
            "NURBS sections need compatible knots, degree and control correspondence",
        ));
    }
    if parameters.iter().any(|x| !x.is_finite())
        || parameters.windows(2).any(|r| r[1] - r[0] <= p.angular())
    {
        return Err(Error::invalid_argument(
            "loft parameters must be finite and strictly increasing",
        ));
    }
    let mut knots = vec![parameters[0]; degree + 1];
    for j in 1..n - degree {
        knots.push(parameters[j..j + degree].iter().sum::<f64>() / degree as f64);
    }
    knots.extend(vec![parameters[n - 1]; degree + 1]);
    let basis = (0..n)
        .map(|i| {
            BCurve::new(
                degree,
                knots.clone(),
                (0..n)
                    .map(|j| Point3::new(if i == j { 1.0 } else { 0.0 }, 0.0, 0.0))
                    .collect(),
                None,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let matrix = parameters
        .iter()
        .map(|&t| basis.iter().map(|b| b.eval(t).point.x).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let max = curves
        .iter()
        .flat_map(|c| {
            c.weights()
                .map_or_else(|| vec![1.0; c.points().len()], |w| w.to_vec())
        })
        .fold(0.0, f64::max);
    let reference = first.points()[0];
    let mut points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..first.points().len() {
        let values = curves
            .iter()
            .map(|curve| {
                let w = curve.weights().map_or(1.0, |ws| ws[i]) / max;
                let q = curve.points()[i] - reference;
                [q.x * w, q.y * w, q.z * w, w]
            })
            .collect::<Vec<_>>();
        let controls = solve(matrix.clone(), values, p)?;
        for c in controls {
            if !c[3].is_finite() || c[3] <= 0.0 {
                return Err(Error::unimplemented(
                    "interpolating loft would need non-positive rational weights",
                ));
            }
            points
                .push(reference + gensolid_math::Vec3::new(c[0] / c[3], c[1] / c[3], c[2] / c[3]));
            weights.push(c[3]);
        }
    }
    BSurface::new(
        first.degree(),
        degree,
        first.knots().to_vec(),
        knots,
        first.points().len(),
        n,
        points,
        Some(weights),
    )
}
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<[f64; 4]>, p: &Precision) -> Result<Vec<[f64; 4]>> {
    let n = a.len();
    for k in 0..n {
        let pivot = (k..n)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .ok_or_else(|| Error::internal("empty interpolation matrix"))?;
        a.swap(k, pivot);
        b.swap(k, pivot);
        let d = a[k][k];
        if !d.is_finite() || d.abs() <= p.angular() {
            return Err(Error::degenerate("ill-conditioned loft interpolation"));
        }
        for value in &mut a[k][k..] {
            *value /= d;
        }
        for v in &mut b[k] {
            *v /= d;
        }
        let pivot_row = a[k].clone();
        let pivot_rhs = b[k];
        for i in 0..n {
            if i == k {
                continue;
            }
            let factor = a[i][k];
            for (value, pivot) in a[i][k..].iter_mut().zip(&pivot_row[k..]) {
                *value -= factor * pivot;
            }
            for (value, pivot) in b[i].iter_mut().zip(pivot_rhs) {
                *value -= factor * pivot;
            }
        }
    }
    Ok(b)
}
#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    #[test]
    fn smooth_rational_loft_interpolates_every_input_section() {
        let curves = (0..4)
            .map(|j| {
                BCurve::new(
                    2,
                    vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                    (0..3)
                        .map(|i| Point3::new(i as f64, j as f64, ((i + j) % 3) as f64))
                        .collect(),
                    Some(vec![1.0, 0.8, 1.0]),
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let p = Precision::DEFAULT;
        let parameters = [0.0, 0.2, 0.7, 1.0];
        let surface = loft_nurbs_sections(&curves, &parameters, 3, &p).unwrap();
        for (curve, t) in curves.iter().zip(parameters) {
            for i in 0..31 {
                let u = i as f64 / 30.0;
                assert!(curve.eval(u).point.distance(surface.eval(u, t).point) <= p.linear());
            }
        }
    }
}
