//! Pure Rust rational curve construction and exact spline editing.
//!
//! These algorithms are original implementations of homogeneous B-spline knot
//! insertion and Bernstein degree elevation; no SMLib or OCCT source is copied.
//! Edits preserve the parameter interval. Degree elevation uses Bezier extraction
//! and may retain extra interior knots instead of performing knot removal.

use crate::BCurve;
use crate::nurbs::{basis_derivatives, find_span};
use gensolid_core::{Error, ErrorCode, Precision, Result};
use gensolid_math::{Point3, Vec3};

/// A rational curve's position and its first two parameter derivatives.
#[derive(Clone, Copy, Debug)]
pub struct CurveDerivs2 {
    /// Evaluated position.
    pub point: Point3,
    /// First parameter derivative.
    pub d1: Vec3,
    /// Second parameter derivative.
    pub d2: Vec3,
}

impl BCurve {
    /// Evaluate two derivatives, clamping finite parameters to the natural range.
    /// Non-finite inputs or unrepresentable results return an error.
    pub fn eval2(&self, t: f64) -> Result<CurveDerivs2> {
        if !t.is_finite() {
            return Err(Error::invalid_argument("curve parameter must be finite"));
        }
        let p = self.degree();
        let t = self.range().normalize(t);
        let span = find_span(self.points().len() - 1, p, t, self.knots());
        let ders = basis_derivatives(span, t, p, 2, self.knots());
        let scale = self.weights().map_or(1.0, |ws| {
            ws[span - p..=span].iter().copied().fold(0.0, f64::max)
        });
        let reference = self.points()[span - p];
        let mut a = [Vec3::ZERO; 3];
        let mut w = [0.0; 3];
        for (j, row) in ders.iter().enumerate() {
            for (k, basis) in row.iter().enumerate() {
                let i = span - p + k;
                let wi = self.weights().map_or(1.0, |ws| ws[i] / scale);
                let coefficient = basis * wi;
                a[j] += (self.points()[i] - reference) * coefficient;
                w[j] += coefficient;
            }
        }
        let c = a[0] / w[0];
        let d1 = (a[1] - c * w[1]) / w[0];
        let out = CurveDerivs2 {
            point: reference + c,
            d1,
            d2: (a[2] - c * w[2] - d1 * (2.0 * w[1])) / w[0],
        };
        if !out.point.is_finite() || !out.d1.is_finite() || !out.d2.is_finite() {
            return Err(Error::degenerate("curve derivatives are not representable"));
        }
        Ok(out)
    }

    /// Insert exact parameter values into a copy of a clamped NURBS curve.
    /// Interior multiplicity may reach degree + 1 (a full parametric break).
    /// Knot order and equality describe the representation, not topology.
    pub fn insert_knots(&self, knots: &[f64]) -> Result<Self> {
        require_clamped(self)?;
        let mut result = self.clone();
        for &knot in knots {
            let range = self.range();
            if !knot.is_finite() || knot <= range.lo || knot >= range.hi {
                return Err(Error::invalid_argument(
                    "inserted knots must be inside the domain",
                ));
            }
            result = insert_one(&result, knot)?;
        }
        Ok(result)
    }

    /// Extract a closed, nonempty subinterval without changing its parameters.
    /// Clamped input is required; no sampling or approximation is performed.
    /// At a discontinuous full break, each new endpoint uses the limit from
    /// inside the retained interval (right at the start, left at the end).
    pub fn trim(&self, interval: [f64; 2]) -> Result<Self> {
        require_clamped(self)?;
        let [lo, hi] = interval;
        let range = self.range();
        if !lo.is_finite() || !hi.is_finite() || lo >= hi || lo < range.lo || hi > range.hi {
            return Err(Error::invalid_argument(
                "trim interval is outside the curve domain",
            ));
        }
        let mut curve = self.clone();
        for boundary in [lo, hi] {
            while multiplicity(curve.knots(), boundary) < curve.degree() + 1 {
                curve = insert_one(&curve, boundary)?;
            }
        }
        let first = curve.knots().partition_point(|&k| k < lo);
        let last = curve.knots().partition_point(|&k| k < hi);
        Self::new(
            curve.degree(),
            curve.knots()[first..last + curve.degree() + 1].to_vec(),
            curve.points()[first..last].to_vec(),
            curve.weights().map(|w| w[first..last].to_vec()),
        )
    }

    /// Exactly raise degree using homogeneous Bezier extraction and elevation.
    /// Continuous joins retain multiplicity `degree`; pre-existing full breaks
    /// remain full breaks. Extra knots are deliberately not removed. Degrees
    /// above 64 return `CapacityExceeded` to bound the extraction workspace.
    pub fn elevate_degree(&self, degree: usize) -> Result<Self> {
        require_clamped(self)?;
        if degree < self.degree() {
            return Err(Error::invalid_argument(
                "elevated degree is below the original degree",
            ));
        }
        if degree == self.degree() {
            return Ok(self.clone());
        }
        if degree > 64 {
            return Err(Error::new(
                ErrorCode::CapacityExceeded,
                "degree elevation is limited to 64",
            ));
        }
        let range = self.range();
        let mut breaks = Vec::new();
        for &k in self.knots() {
            if k > range.lo && k < range.hi && breaks.last().is_none_or(|&b| !same(b, k)) {
                breaks.push(k);
            }
        }
        let mut extracted = self.clone();
        for &k in &breaks {
            while multiplicity(extracted.knots(), k) < self.degree() + 1 {
                extracted = insert_one(&extracted, k)?;
            }
        }
        let (reference, controls) = homogeneous(&extracted)?;
        let mut elevated = Vec::new();
        let mut knots = vec![range.lo; degree + 1];
        for (index, segment) in controls.chunks_exact(self.degree() + 1).enumerate() {
            let mut current = segment.to_vec();
            for d in self.degree()..degree {
                let mut next = Vec::with_capacity(d + 2);
                next.push(current[0]);
                for i in 1..=d {
                    next.push(blend(current[i], current[i - 1], i as f64 / (d + 1) as f64));
                }
                next.push(current[d]);
                current = next;
            }
            let connected =
                index > 0 && multiplicity(self.knots(), breaks[index - 1]) <= self.degree();
            if index > 0 {
                knots.extend(std::iter::repeat_n(
                    breaks[index - 1],
                    degree + usize::from(!connected),
                ));
            }
            elevated.extend(current.into_iter().skip(usize::from(connected)));
        }
        knots.extend(std::iter::repeat_n(range.hi, degree + 1));
        from_homogeneous(degree, knots, reference, elevated, self.weights().is_some())
    }
}

fn same(a: f64, b: f64) -> bool {
    // IEEE equality of knot values is representation identity, not a geometric
    // coincidence decision; signed zero denotes the same knot.
    a.partial_cmp(&b) == Some(std::cmp::Ordering::Equal)
}

fn multiplicity(knots: &[f64], value: f64) -> usize {
    knots.iter().filter(|&&k| same(k, value)).count()
}

fn require_clamped(curve: &BCurve) -> Result<()> {
    let r = curve.range();
    if curve.knots()[..=curve.degree()]
        .iter()
        .any(|&k| !same(k, r.lo))
        || curve.knots()[curve.points().len()..]
            .iter()
            .any(|&k| !same(k, r.hi))
    {
        return Err(Error::unimplemented(
            "spline edits currently require clamped end knots",
        ));
    }
    Ok(())
}

fn homogeneous(curve: &BCurve) -> Result<(Point3, Vec<[f64; 4]>)> {
    let reference = curve.points()[0];
    let scale = curve
        .weights()
        .map_or(1.0, |w| w.iter().copied().fold(0.0, f64::max));
    let mut controls = Vec::with_capacity(curve.points().len());
    for (i, &point) in curve.points().iter().enumerate() {
        let w = curve.weights().map_or(1.0, |ws| ws[i] / scale);
        let q = (point - reference) * w;
        if !q.is_finite() || w <= 0.0 {
            return Err(Error::degenerate(
                "homogeneous curve coordinates are not representable",
            ));
        }
        controls.push([q.x, q.y, q.z, w]);
    }
    Ok((reference, controls))
}

fn blend(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    std::array::from_fn(|i| a[i] * (1.0 - t) + b[i] * t)
}

fn from_homogeneous(
    degree: usize,
    knots: Vec<f64>,
    reference: Point3,
    controls: Vec<[f64; 4]>,
    rational: bool,
) -> Result<BCurve> {
    if controls
        .iter()
        .any(|q| q.iter().any(|v| !v.is_finite()) || q[3] <= 0.0)
    {
        return Err(Error::degenerate(
            "edited homogeneous coordinates are not representable",
        ));
    }
    let points = controls
        .iter()
        .map(|q| reference + Vec3::new(q[0] / q[3], q[1] / q[3], q[2] / q[3]))
        .collect();
    let weights = rational.then(|| controls.iter().map(|q| q[3]).collect());
    BCurve::new(degree, knots, points, weights)
}

fn insert_one(curve: &BCurve, knot: f64) -> Result<BCurve> {
    let p = curve.degree();
    let s = multiplicity(curve.knots(), knot);
    if s > p {
        return Err(Error::invalid_argument(
            "inserted knot exceeds degree + 1 multiplicity",
        ));
    }
    let k = curve.knots().partition_point(|&x| x <= knot) - 1;
    let (reference, controls) = homogeneous(curve)?;
    let mut out = Vec::with_capacity(controls.len() + 1);
    out.extend_from_slice(&controls[..=k - p]);
    for i in k - p + 1..=k - s {
        let alpha = (knot - curve.knots()[i]) / (curve.knots()[i + p] - curve.knots()[i]);
        out.push(blend(controls[i - 1], controls[i], alpha));
    }
    out.extend_from_slice(&controls[k - s..]);
    let mut knots = curve.knots().to_vec();
    knots.insert(k + 1, knot);
    from_homogeneous(p, knots, reference, out, curve.weights().is_some())
}

/// Create a bounded degree-one line with parameter interval `[0, 1]`.
pub fn line_segment(start: Point3, end: Point3, precision: &Precision) -> Result<BCurve> {
    if !start.is_finite() || !end.is_finite() {
        return Err(Error::invalid_argument("line endpoints must be finite"));
    }
    let length = start.distance(end);
    if !length.is_finite() {
        return Err(Error::invalid_argument(
            "line displacement must be representable",
        ));
    }
    if length <= precision.degenerate_length() {
        return Err(Error::degenerate(
            "line segment is shorter than session precision",
        ));
    }
    BCurve::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![start, end], None)
}

/// Exact rational circle in the XY plane, parameter interval `[0, TAU]`.
/// NURBS parameters are not uniform angles between the quadrant endpoints.
pub fn circle(center: Point3, radius: f64, precision: &Precision) -> Result<BCurve> {
    ellipse(center, [radius, radius], precision)
}

/// Exact rational ellipse in the XY plane, parameter interval `[0, TAU]`.
pub fn ellipse(center: Point3, radii: [f64; 2], precision: &Precision) -> Result<BCurve> {
    conic(center, radii, [0.0, std::f64::consts::TAU], precision)
}

/// Exact counterclockwise XY circular arc, at most one full revolution.
/// Angles set the end points and knot interval; the interior rational parameter
/// is not uniform angular speed.
pub fn arc(center: Point3, radius: f64, angles: [f64; 2], precision: &Precision) -> Result<BCurve> {
    conic(center, [radius, radius], angles, precision)
}

fn conic(center: Point3, radii: [f64; 2], angles: [f64; 2], p: &Precision) -> Result<BCurve> {
    if !center.is_finite()
        || radii.iter().any(|r| !r.is_finite())
        || angles.iter().any(|a| !a.is_finite())
    {
        return Err(Error::invalid_argument("conic parameters must be finite"));
    }
    if radii.iter().any(|&r| r <= p.degenerate_length()) {
        return Err(Error::degenerate(
            "conic radii must exceed session precision",
        ));
    }
    let sweep = angles[1] - angles[0];
    if !sweep.is_finite() || sweep <= p.angular() || sweep > std::f64::consts::TAU {
        return Err(Error::invalid_argument(
            "arc angle span must be in (angular precision, TAU]",
        ));
    }
    let segments = (sweep / std::f64::consts::FRAC_PI_2).ceil() as usize;
    let mut points = Vec::with_capacity(2 * segments + 1);
    let mut weights = Vec::with_capacity(2 * segments + 1);
    let mut knots = vec![angles[0]; 3];
    // Reduce the original angle through sin_cos, not remainder by rounded TAU:
    // the latter loses the phase of large, otherwise valid starting angles.
    let (start_sin, start_cos) = angles[0].sin_cos();
    for i in 0..segments {
        let a = angles[0] + sweep * (i as f64 / segments as f64);
        let b = angles[0] + sweep * ((i + 1) as f64 / segments as f64);
        let geometric_a = sweep * (i as f64 / segments as f64);
        let geometric_b = sweep * ((i + 1) as f64 / segments as f64);
        let middle = geometric_a + (geometric_b - geometric_a) * 0.5;
        let w = ((geometric_b - geometric_a) * 0.5).cos();
        let pos = |t: f64, weight: f64| {
            let (s, c) = t.sin_cos();
            center
                + Vec3::new(
                    radii[0] * (start_cos * c - start_sin * s) / weight,
                    radii[1] * (start_sin * c + start_cos * s) / weight,
                    0.0,
                )
        };
        if i == 0 {
            points.push(pos(geometric_a, 1.0));
            weights.push(1.0);
        }
        points.push(pos(middle, w));
        weights.push(w);
        points.push(pos(geometric_b, 1.0));
        weights.push(1.0);
        if b <= a {
            return Err(Error::invalid_argument(
                "arc knots cannot represent the subdivisions",
            ));
        }
        knots.extend(std::iter::repeat_n(
            b,
            if i + 1 == segments { 3 } else { 2 },
        ));
    }
    BCurve::new(2, knots, points, Some(weights))
}
