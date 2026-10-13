//! Non-uniform rational B-spline curves and surfaces.
//!
//! Evaluation follows Piegl & Tiller, *The NURBS Book* (2nd ed.): `find_span` is A2.1 and
//! `basis_derivatives` is A2.3; rational derivatives use the quotient rule on the
//! homogeneous form.

use gensolid_core::{Error, Result};
use gensolid_math::{Point3, Vec3};

use crate::curve::CurveDerivs;
use crate::param::ParamRange;
use crate::surface::SurfaceDerivs;

/// Index `i` of the knot span `[U_i, U_{i+1})` containing `u`, for `n + 1` control points
/// of degree `p`. `u` is clamped to the domain `[U_p, U_{n+1}]`.
fn find_span(n: usize, p: usize, u: f64, knots: &[f64]) -> usize {
    if u >= knots[n + 1] {
        // Last non-empty span, so the domain end evaluates like its interior.
        let mut span = n;
        while span > p && knots[span] >= knots[n + 1] {
            span -= 1;
        }
        return span;
    }
    if u <= knots[p] {
        let mut span = p;
        while span < n && knots[span + 1] <= knots[p] {
            span += 1;
        }
        return span;
    }
    let (mut low, mut high) = (p, n + 1);
    let mut mid = (low + high) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

/// Non-zero basis functions and their derivatives up to order `n` at `u` in span `span`.
/// Row `k` holds the `k`-th derivatives of `N_{span-p..=span, p}`.
fn basis_derivatives(span: usize, u: f64, p: usize, n: usize, knots: &[f64]) -> Vec<Vec<f64>> {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    let mut ders = vec![vec![0.0; p + 1]; n + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let pi = p as isize;
    let mut a = [vec![0.0; p + 1], vec![0.0; p + 1]];
    for r in 0..=p {
        let ri = r as isize;
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=n.min(p) {
            let ki = k as isize;
            let rk = ri - ki;
            let pk = (pi - ki) as usize;
            let mut d = 0.0;
            if rk >= 0 {
                let rk = rk as usize;
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if ri - 1 <= pi - ki { k - 1 } else { p - r };
            for j in j1..=j2 {
                let idx = (rk + j as isize) as usize;
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx];
                d += a[s2][j] * ndu[idx][pk];
            }
            if ri <= pi - ki {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as f64;
    for (k, row) in ders.iter_mut().enumerate().skip(1).take(n.min(p)) {
        for value in row {
            *value *= factor;
        }
        factor *= (p - k) as f64;
    }
    ders
}

fn validate_knots(degree: usize, count: usize, knots: &[f64], direction: &str) -> Result<()> {
    if degree == 0 {
        return Err(Error::invalid_argument(format!(
            "{direction} degree must be at least 1"
        )));
    }
    if count <= degree {
        return Err(Error::invalid_argument(format!(
            "{direction}: need at least degree + 1 control points"
        )));
    }
    if count.checked_add(degree).and_then(|n| n.checked_add(1)) != Some(knots.len()) {
        return Err(Error::invalid_argument(format!(
            "{direction}: knot count must equal control point count + degree + 1"
        )));
    }
    if knots.iter().any(|k| !k.is_finite()) || knots.windows(2).any(|w| w[0] > w[1]) {
        return Err(Error::invalid_argument(format!(
            "{direction}: knots must be finite and non-decreasing"
        )));
    }
    if knots[degree] >= knots[count] {
        return Err(Error::degenerate(format!(
            "{direction}: parameter domain is empty"
        )));
    }
    if !(knots[count] - knots[degree]).is_finite() {
        return Err(Error::invalid_argument(format!(
            "{direction}: parameter domain width must be finite"
        )));
    }
    // Knots are non-decreasing here, so `>=` means equal.
    if knots.windows(degree + 2).any(|w| w[0] >= w[degree + 1]) {
        return Err(Error::invalid_argument(format!(
            "{direction}: knot multiplicity exceeds degree + 1"
        )));
    }
    Ok(())
}

fn validate_points_and_weights(points: &[Point3], weights: Option<&[f64]>) -> Result<()> {
    if points.iter().any(|p| !p.is_finite()) {
        return Err(Error::invalid_argument("control points must be finite"));
    }
    if let Some(weights) = weights {
        if weights.len() != points.len() {
            return Err(Error::invalid_argument(
                "weight count must equal control point count",
            ));
        }
        if weights.iter().any(|w| !(w.is_finite() && *w > 0.0)) {
            return Err(Error::invalid_argument(
                "weights must be finite and positive",
            ));
        }
    }
    Ok(())
}

/// NURBS curve. Polynomial when `weights` is `None`.
#[derive(Clone, Debug, PartialEq)]
pub struct BCurve {
    degree: usize,
    knots: Vec<f64>,
    points: Vec<Point3>,
    weights: Option<Vec<f64>>,
}

impl BCurve {
    /// Creates a curve after validating degree, knots, points and weights.
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        points: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self> {
        validate_knots(degree, points.len(), &knots, "curve")?;
        validate_points_and_weights(&points, weights.as_deref())?;
        Ok(Self {
            degree,
            knots,
            points,
            weights,
        })
    }

    /// Degree.
    pub const fn degree(&self) -> usize {
        self.degree
    }

    /// Knot vector.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// Control points.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// Weights, if rational.
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// Parameter domain `[U_p, U_{n+1}]`.
    pub fn range(&self) -> ParamRange {
        ParamRange::bounded(self.knots[self.degree], self.knots[self.points.len()])
    }

    /// Point and first derivative at `t` (clamped to the domain).
    pub fn eval(&self, t: f64) -> CurveDerivs {
        let p = self.degree;
        let n = self.points.len() - 1;
        let t = self.range().normalize(t);
        let span = find_span(n, p, t, &self.knots);
        let ders = basis_derivatives(span, t, p, 1, &self.knots);
        // Rational geometry is invariant under a common weight scale. Normalising
        // the active weights prevents harmless common scales from overflowing the
        // homogeneous coordinates (or underflowing their denominator).
        let weight_scale = self.weights.as_ref().map_or(1.0, |weights| {
            weights[span - p..=span].iter().copied().fold(0.0, f64::max)
        });
        // Evaluate relative to an active control point. This avoids subtracting
        // large translated homogeneous coordinates in the quotient derivative.
        let reference = self.points[span - p];
        let (mut a, mut da) = (Vec3::ZERO, Vec3::ZERO);
        let (mut w, mut dw) = (0.0, 0.0);
        for (k, (&n0, &n1)) in ders[0].iter().zip(&ders[1]).enumerate() {
            let i = span - p + k;
            let wi = self.weights.as_ref().map_or(1.0, |ws| ws[i] / weight_scale);
            let pi = (self.points[i] - reference) * wi;
            a += pi * n0;
            da += pi * n1;
            w += wi * n0;
            dw += wi * n1;
        }
        let c = a / w;
        CurveDerivs {
            point: reference + c,
            d1: (da - c * dw) / w,
        }
    }
}

/// NURBS surface. Control points are stored row-major: index `i * v_count + j` where `i`
/// runs along `u` and `j` along `v`.
#[derive(Clone, Debug, PartialEq)]
pub struct BSurface {
    u_degree: usize,
    v_degree: usize,
    u_knots: Vec<f64>,
    v_knots: Vec<f64>,
    u_count: usize,
    v_count: usize,
    points: Vec<Point3>,
    weights: Option<Vec<f64>>,
}

impl BSurface {
    /// Creates a surface after validating degrees, knots, points and weights.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        u_degree: usize,
        v_degree: usize,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        u_count: usize,
        v_count: usize,
        points: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self> {
        if u_count.checked_mul(v_count) != Some(points.len()) {
            return Err(Error::invalid_argument(
                "control net size must equal u_count * v_count",
            ));
        }
        validate_knots(u_degree, u_count, &u_knots, "u")?;
        validate_knots(v_degree, v_count, &v_knots, "v")?;
        validate_points_and_weights(&points, weights.as_deref())?;
        Ok(Self {
            u_degree,
            v_degree,
            u_knots,
            v_knots,
            u_count,
            v_count,
            points,
            weights,
        })
    }

    /// Degrees in `u` and `v`.
    pub const fn degrees(&self) -> (usize, usize) {
        (self.u_degree, self.v_degree)
    }

    /// Control net dimensions in `u` and `v`.
    pub const fn counts(&self) -> (usize, usize) {
        (self.u_count, self.v_count)
    }

    /// Knot vectors in the `u` and `v` directions.
    pub fn knots(&self) -> (&[f64], &[f64]) {
        (&self.u_knots, &self.v_knots)
    }

    /// Row-major control points.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// Row-major weights, if rational.
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// Parameter domains in `u` and `v`.
    pub fn ranges(&self) -> (ParamRange, ParamRange) {
        (
            ParamRange::bounded(self.u_knots[self.u_degree], self.u_knots[self.u_count]),
            ParamRange::bounded(self.v_knots[self.v_degree], self.v_knots[self.v_count]),
        )
    }

    /// Point and first partial derivatives at `(u, v)` (clamped to the domain).
    pub fn eval(&self, u: f64, v: f64) -> SurfaceDerivs {
        let (ur, vr) = self.ranges();
        let (u, v) = (ur.normalize(u), vr.normalize(v));
        let (pu, pv) = (self.u_degree, self.v_degree);
        let su = find_span(self.u_count - 1, pu, u, &self.u_knots);
        let sv = find_span(self.v_count - 1, pv, v, &self.v_knots);
        let nu = basis_derivatives(su, u, pu, 1, &self.u_knots);
        let nv = basis_derivatives(sv, v, pv, 1, &self.v_knots);
        let mut weight_scale: f64 = 1.0;
        if let Some(weights) = &self.weights {
            weight_scale = 0.0;
            for i in su - pu..=su {
                for j in sv - pv..=sv {
                    weight_scale = weight_scale.max(weights[i * self.v_count + j]);
                }
            }
        }
        let reference = self.points[(su - pu) * self.v_count + sv - pv];
        let (mut a, mut au, mut av) = (Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
        let (mut w, mut wu, mut wv) = (0.0, 0.0, 0.0);
        for (k, (&u0, &u1)) in nu[0].iter().zip(&nu[1]).enumerate() {
            for (l, (&v0, &v1)) in nv[0].iter().zip(&nv[1]).enumerate() {
                let idx = (su - pu + k) * self.v_count + (sv - pv + l);
                let wi = self
                    .weights
                    .as_ref()
                    .map_or(1.0, |ws| ws[idx] / weight_scale);
                let pw = (self.points[idx] - reference) * wi;
                let (b, bu, bv) = (u0 * v0, u1 * v0, u0 * v1);
                a += pw * b;
                au += pw * bu;
                av += pw * bv;
                w += wi * b;
                wu += wi * bu;
                wv += wi * bv;
            }
        }
        let s = a / w;
        SurfaceDerivs {
            point: reference + s,
            du: (au - s * wu) / w,
            dv: (av - s * wv) / w,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_1_SQRT_2;

    use gensolid_core::Precision;
    use proptest::prelude::*;

    use super::*;

    /// Quarter circle of radius 1 as a rational quadratic.
    fn quarter_circle() -> BCurve {
        BCurve::new(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            Some(vec![1.0, FRAC_1_SQRT_2, 1.0]),
        )
        .unwrap()
    }

    #[test]
    fn rejects_bad_knot_vectors() {
        let pts = vec![Point3::ORIGIN, Point3::new(1.0, 0.0, 0.0)];
        assert!(BCurve::new(1, vec![0.0, 0.0, 1.0], pts.clone(), None).is_err());
        assert!(BCurve::new(1, vec![0.0, 1.0, 0.5, 1.0], pts.clone(), None).is_err());
        assert!(BCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], pts, None).is_ok());
    }

    #[test]
    fn linear_bcurve_matches_segment() {
        let c = BCurve::new(
            1,
            vec![0.0, 0.0, 2.0, 2.0],
            vec![Point3::ORIGIN, Point3::new(4.0, 2.0, 0.0)],
            None,
        )
        .unwrap();
        let d = c.eval(0.5);
        assert_eq!(d.point, Point3::new(1.0, 0.5, 0.0));
        assert_eq!(d.d1, Vec3::new(2.0, 1.0, 0.0));
    }

    proptest! {
        #[test]
        fn rational_quarter_circle_is_exact(t in 0.0f64..=1.0) {
            let tol = Precision::DEFAULT.linear();
            let d = quarter_circle().eval(t);
            prop_assert!((d.point.to_vec().length() - 1.0).abs() <= tol);
            // The tangent of a circle is perpendicular to its radius.
            prop_assert!(d.point.to_vec().dot(d.d1).abs() <= tol * d.d1.length());
        }

        #[test]
        fn curve_derivative_matches_finite_difference(t in 0.05f64..0.95) {
            let c = BCurve::new(
                3,
                vec![0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0],
                vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 2.0, 0.0),
                    Point3::new(2.0, -1.0, 1.0),
                    Point3::new(3.0, 1.0, 2.0),
                    Point3::new(4.0, 0.0, 0.0),
                ],
                Some(vec![1.0, 2.0, 0.5, 1.5, 1.0]),
            )
            .unwrap();
            let h = Precision::DEFAULT.linear().sqrt();
            let fd = (c.eval(t + h).point - c.eval(t - h).point) / (2.0 * h);
            let d1 = c.eval(t).d1;
            prop_assert!((fd - d1).length() <= h * (1.0 + d1.length()));
        }
    }

    #[test]
    fn bilinear_surface_eval() {
        let s = BSurface::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            2,
            2,
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 1.0),
            ],
            None,
        )
        .unwrap();
        let d = s.eval(0.5, 0.5);
        assert_eq!(d.point, Point3::new(0.5, 0.5, 0.25));
        assert_eq!(d.du, Vec3::new(1.0, 0.0, 0.5));
        assert_eq!(d.dv, Vec3::new(0.0, 1.0, 0.5));
    }
}
