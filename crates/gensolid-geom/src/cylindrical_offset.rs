//! Native Rust offsets of a single rational quadratic cylindrical patch.
//!
//! Construction proposes a changed radius. Acceptance is separate: outward-rounded
//! Bernstein coefficients bound the actual source/output displacement, its length,
//! its orthogonality to both source tangents, and both Jacobian orientations on the
//! entire parameter rectangle. No point samples are used as an error certificate.

use crate::BSurface;
use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Frame, Interval, Point3, Sign, Vec3};

fn unsupported(message: &str) -> Error {
    Error::unimplemented(message)
}

fn reciprocal(value: Interval) -> Result<Interval> {
    if value.sign() != Some(Sign::Positive) && value.sign() != Some(Sign::Negative) {
        return Err(unsupported("offset certificate denominator contains zero"));
    }
    Interval::new((1.0 / value.hi()).next_down(), (1.0 / value.lo()).next_up())
        .ok_or_else(|| unsupported("offset certificate quotient is unbounded"))
}

fn coordinates(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn shape(surface: &BSurface) -> Result<[f64; 3]> {
    if surface.degrees() != (2, 1) || surface.counts() != (3, 2) {
        return Err(unsupported(
            "cylindrical offset requires a 3 by 2 quadratic/linear patch",
        ));
    }
    let (u, v) = surface.ranges();
    let (uk, vk) = surface.knots();
    if uk[..3].iter().any(|k| k.to_bits() != u.lo.to_bits())
        || uk[3..].iter().any(|k| k.to_bits() != u.hi.to_bits())
        || vk[..2].iter().any(|k| k.to_bits() != v.lo.to_bits())
        || vk[2..].iter().any(|k| k.to_bits() != v.hi.to_bits())
    {
        return Err(unsupported("cylindrical offset requires one clamped span"));
    }
    let weights = surface
        .weights()
        .ok_or_else(|| unsupported("polynomial curved offsets are not supported"))?;
    if (0..3).any(|i| weights[2 * i].to_bits() != weights[2 * i + 1].to_bits())
        || weights[0].to_bits() != weights[4].to_bits()
    {
        return Err(unsupported(
            "cylindrical offset requires repeated arc weights along the axis",
        ));
    }
    Ok([weights[0], weights[2], weights[4]])
}

struct Cylinder {
    center: Point3,
    axis: Vec3,
    radius: f64,
    orientation: f64,
}

// Recognition selects the supported family and a proposal. The independent
// polynomial certificate below, not these numerical residuals, admits output.
fn cylinder(surface: &BSurface, p: &Precision) -> Result<Cylinder> {
    let weights = shape(surface)?;
    let w = weights[1] / weights[0];
    let squared_sine = 1.0 - w * w;
    if !w.is_finite() || w < std::f64::consts::FRAC_1_SQRT_2 || squared_sine <= p.angular() {
        return Err(unsupported(
            "cylindrical offset supports nondegenerate arcs of at most a quarter turn",
        ));
    }
    let points = surface.points();
    let direction = points[1] - points[0];
    let axis = direction
        .normalized(p.degenerate_length())
        .ok_or_else(|| Error::degenerate("cylinder extrusion is degenerate"))?;
    for i in 1..3 {
        if ((points[2 * i + 1] - points[2 * i]) - direction).length() > p.solver() {
            return Err(unsupported("surface is not a constant linear extrusion"));
        }
    }
    let middle = points[0] + (points[4] - points[0]) * 0.5;
    let center = middle - (points[2] - middle) * (w * w / squared_sine);
    let radius = (points[0] - center).length();
    if !radius.is_finite() || radius <= p.degenerate_length() {
        return Err(Error::degenerate("cylinder radius is degenerate"));
    }
    if ((points[4] - center).length() - radius).abs() > p.solver()
        || ((points[2] - center).length() * w - radius).abs() > p.solver()
        || [0, 2, 4]
            .iter()
            .any(|&i| (points[i] - center).dot(axis).abs() > p.solver())
    {
        return Err(unsupported(
            "extruded quadratic is not a perpendicular circular arc",
        ));
    }
    let (u, v) = surface.ranges();
    let mid = surface.eval(u.lo + (u.hi - u.lo) * 0.5, v.lo + (v.hi - v.lo) * 0.5);
    let relative = mid.point - center;
    let radial = relative - axis * relative.dot(axis);
    let normal = mid.du.cross(mid.dv);
    let sign = Interval::point(normal.x) * Interval::point(radial.x)
        + Interval::point(normal.y) * Interval::point(radial.y)
        + Interval::point(normal.z) * Interval::point(radial.z);
    let orientation = match sign.sign() {
        Some(Sign::Positive) => 1.0,
        Some(Sign::Negative) => -1.0,
        _ => return Err(Error::degenerate("cylinder orientation is unresolved")),
    };
    Ok(Cylinder {
        center,
        axis,
        radius,
        orientation,
    })
}

/// Offset a single clamped circular-arc extrusion, with arc angle at most π/2.
/// The signed distance follows `du × dv`. The exact rational basis is retained.
/// Radius collapse, orientation reversal and uncertified results are rejected.
pub fn offset_cylindrical(source: &BSurface, distance: f64, p: &Precision) -> Result<BSurface> {
    if !distance.is_finite() {
        return Err(Error::invalid_argument("offset distance must be finite"));
    }
    let c = cylinder(source, p)?;
    let radius = c.radius + c.orientation * distance;
    if !radius.is_finite() || radius <= p.degenerate_length() {
        return Err(Error::degenerate(
            "cylindrical offset collapses or crosses the axis",
        ));
    }
    let scale = radius / c.radius;
    let (uk, vk) = source.knots();
    let output = BSurface::new(
        2,
        1,
        uk.to_vec(),
        vk.to_vec(),
        3,
        2,
        source
            .points()
            .iter()
            .map(|&point| {
                let relative = point - c.center;
                let axial = c.axis * relative.dot(c.axis);
                c.center + axial + (relative - axial) * scale
            })
            .collect(),
        source.weights().map(<[f64]>::to_vec),
    )?;
    certify_offset(source, &output, distance, p)?;
    Ok(output)
}

#[derive(Clone)]
struct Polynomial {
    u: usize,
    v: usize,
    coefficients: Vec<Interval>,
}

fn binomial(n: usize, k: usize) -> usize {
    let k = k.min(n - k);
    (0..k).fold(1, |value, i| value * (n - i) / (i + 1))
}

impl Polynomial {
    fn elevated(&self, u: usize, v: usize) -> Result<Self> {
        if u < self.u || v < self.v {
            return Err(Error::internal(
                "cannot lower certificate polynomial degree",
            ));
        }
        let mut result = self.clone();
        for along_u in [true, false] {
            while if along_u { result.u < u } else { result.v < v } {
                let old = result;
                let du = old.u + usize::from(along_u);
                let dv = old.v + usize::from(!along_u);
                let degree = if along_u { du } else { dv };
                let inverse = reciprocal(Interval::point(degree as f64))?;
                let mut coefficients = Vec::new();
                for i in 0..=du {
                    for j in 0..=dv {
                        let k = if along_u { i } else { j };
                        if k == 0 {
                            coefficients.push(old.coefficients[i * (old.v + 1) + j]);
                        } else if k == degree {
                            let a = i - usize::from(along_u);
                            let b = j - usize::from(!along_u);
                            coefficients.push(old.coefficients[a * (old.v + 1) + b]);
                        } else {
                            let alpha = Interval::point(k as f64) * inverse;
                            let previous = (i - usize::from(along_u)) * (old.v + 1) + j
                                - usize::from(!along_u);
                            let next = i * (old.v + 1) + j;
                            coefficients.push(
                                old.coefficients[previous] * alpha
                                    + old.coefficients[next] * (Interval::point(1.0) - alpha),
                            );
                        }
                    }
                }
                result = Self {
                    u: du,
                    v: dv,
                    coefficients,
                };
            }
        }
        Ok(result)
    }
    fn scaled(&self, factor: Interval) -> Self {
        Self {
            u: self.u,
            v: self.v,
            coefficients: self.coefficients.iter().map(|&c| c * factor).collect(),
        }
    }
    fn add(&self, other: &Self) -> Result<Self> {
        if self.u != other.u || self.v != other.v {
            return Err(Error::internal("offset polynomial degree mismatch"));
        }
        Ok(Self {
            u: self.u,
            v: self.v,
            coefficients: self
                .coefficients
                .iter()
                .zip(&other.coefficients)
                .map(|(&a, &b)| a + b)
                .collect(),
        })
    }
    fn sub(&self, other: &Self) -> Result<Self> {
        self.add(&other.scaled(Interval::point(-1.0)))
    }
    fn derivative(&self, along_u: bool) -> Self {
        let degree = if along_u { self.u } else { self.v };
        if degree == 0 {
            return Self {
                u: self.u,
                v: self.v,
                coefficients: vec![Interval::point(0.0); self.coefficients.len()],
            };
        }
        let u = self.u - usize::from(along_u);
        let v = self.v - usize::from(!along_u);
        let mut coefficients = Vec::new();
        for i in 0..=u {
            for j in 0..=v {
                let start = i * (self.v + 1) + j;
                let end = start + if along_u { self.v + 1 } else { 1 };
                coefficients.push(
                    (self.coefficients[end] - self.coefficients[start])
                        * Interval::point(degree as f64),
                );
            }
        }
        Self { u, v, coefficients }
    }
    fn mul(&self, other: &Self) -> Result<Self> {
        let (u, v) = (self.u + other.u, self.v + other.v);
        let mut coefficients = vec![Interval::point(0.0); (u + 1) * (v + 1)];
        for i in 0..=self.u {
            for j in 0..=self.v {
                for k in 0..=other.u {
                    for l in 0..=other.v {
                        // Degrees are bounded by (14,6) in this certificate;
                        // these integer binomials and products are exactly f64-representable.
                        let numerator = binomial(self.u, i)
                            * binomial(other.u, k)
                            * binomial(self.v, j)
                            * binomial(other.v, l);
                        let denominator = binomial(u, i + k) * binomial(v, j + l);
                        let factor = Interval::point(numerator as f64)
                            * reciprocal(Interval::point(denominator as f64))?;
                        let index = (i + k) * (v + 1) + j + l;
                        coefficients[index] = coefficients[index]
                            + self.coefficients[i * (self.v + 1) + j]
                                * other.coefficients[k * (other.v + 1) + l]
                                * factor;
                    }
                }
            }
        }
        Ok(Self { u, v, coefficients })
    }
    fn lower(&self) -> f64 {
        self.coefficients
            .iter()
            .map(|c| c.lo())
            .fold(f64::INFINITY, f64::min)
    }
    fn absolute_upper(&self) -> f64 {
        self.coefficients
            .iter()
            .map(|c| c.lo().abs().max(c.hi().abs()))
            .fold(0.0, f64::max)
    }
}

type VectorPolynomial = [Polynomial; 3];

fn dot(a: &VectorPolynomial, b: &VectorPolynomial) -> Result<Polynomial> {
    a[0].mul(&b[0])?
        .add(&a[1].mul(&b[1])?)?
        .add(&a[2].mul(&b[2])?)
}

fn cross(a: &VectorPolynomial, b: &VectorPolynomial) -> Result<VectorPolynomial> {
    Ok([
        a[1].mul(&b[2])?.sub(&a[2].mul(&b[1])?)?,
        a[2].mul(&b[0])?.sub(&a[0].mul(&b[2])?)?,
        a[0].mul(&b[1])?.sub(&a[1].mul(&b[0])?)?,
    ])
}

fn derivative_numerator(
    h: &VectorPolynomial,
    w: &Polynomial,
    along_u: bool,
) -> Result<VectorPolynomial> {
    let dw = w.derivative(along_u);
    let component = |i: usize| h[i].derivative(along_u).mul(w)?.sub(&h[i].mul(&dw)?);
    Ok([component(0)?, component(1)?, component(2)?])
}

fn sqrt_upper(value: Interval) -> Result<f64> {
    let result = value.hi().sqrt().next_up();
    if !result.is_finite() {
        return Err(unsupported("offset certificate norm is unbounded"));
    }
    Ok(result)
}

fn norm_upper(vector: &VectorPolynomial) -> Result<f64> {
    sqrt_upper(
        vector
            .iter()
            .map(|p| Interval::point(p.absolute_upper()).sqr())
            .fold(Interval::point(0.0), |a, b| a + b),
    )
}

fn positive(polynomial: &Polynomial, message: &str) -> Result<()> {
    if polynomial
        .coefficients
        .iter()
        .all(|c| c.sign() == Some(Sign::Positive))
    {
        Ok(())
    } else {
        Err(unsupported(message))
    }
}

fn weight_polynomial(surface: &BSurface) -> Result<Polynomial> {
    let weights = surface
        .weights()
        .ok_or_else(|| unsupported("cylindrical offset needs rational weights"))?;
    let maximum = weights.iter().copied().fold(0.0, f64::max);
    let scale = reciprocal(Interval::point(maximum))?;
    Ok(Polynomial {
        u: 2,
        v: 1,
        coefficients: weights
            .iter()
            .map(|&w| Interval::point(w) * scale)
            .collect(),
    })
}

fn homogeneous(points: &[Point3], origin: Point3, w: &Polynomial) -> VectorPolynomial {
    let origin = coordinates(origin);
    std::array::from_fn(|k| Polynomial {
        u: 2,
        v: 1,
        coefficients: points
            .iter()
            .enumerate()
            .map(|(i, &point)| {
                (Interval::point(coordinates(point)[k]) - Interval::point(origin[k]))
                    * w.coefficients[i]
            })
            .collect(),
    })
}

fn project(vector: &VectorPolynomial, direction: Vec3) -> Result<Polynomial> {
    vector[0]
        .scaled(Interval::point(direction.x))
        .add(&vector[1].scaled(Interval::point(direction.y)))?
        .add(&vector[2].scaled(Interval::point(direction.z)))
}

/// A tighter whole-domain proof of the existing regular graph criterion for a
/// supported circular extrusion. Uses the same midpoint frame and normalized
/// parameters as the modelling layer: `j00 > linear`, `j11 > linear`, and
/// `j00*j11 - ((j01+j10)/2)^2 > linear^2`. Positive rational denominators are
/// cleared algebraically; all numerator coefficients carry outward rounding.
/// The caller retains its control-hull/size and continuity checks.
pub fn certify_cylindrical_graph(source: &BSurface, p: &Precision) -> Result<()> {
    let _ = cylinder(source, p)?;
    let (u, v) = source.ranges();
    let mid = source.eval((u.lo + u.hi) * 0.5, (v.lo + v.hi) * 0.5);
    let frame = Frame::new(mid.point, mid.du.cross(mid.dv), mid.du, p)?;
    let w = weight_polynomial(source)?;
    let h = homogeneous(source.points(), source.points()[0], &w);
    // Bernstein differentiation is in normalized [0,1] coordinates. This is
    // exactly the natural derivative multiplied by its positive knot span.
    let qu = derivative_numerator(&h, &w, true)?;
    let qv = derivative_numerator(&h, &w, false)?;
    let j00 = project(&qu, frame.x())?.elevated(4, 2)?;
    let j11 = project(&qv, frame.y())?.elevated(4, 2)?;
    let off = project(&qv, frame.x())?
        .elevated(4, 2)?
        .add(&project(&qu, frame.y())?.elevated(4, 2)?)?
        .scaled(Interval::point(0.5));
    let w2 = w.mul(&w)?;
    let linear = Interval::point(p.linear());
    positive(
        &j00.sub(&w2.scaled(linear))?,
        "cylindrical graph U margin is not certified",
    )?;
    positive(
        &j11.sub(&w2.scaled(linear))?,
        "cylindrical graph V margin is not certified",
    )?;
    let determinant = j00
        .mul(&j11)?
        .sub(&off.mul(&off)?)?
        .sub(&w2.mul(&w2)?.scaled(linear.sqr()))?;
    positive(
        &determinant,
        "cylindrical graph determinant margin is not certified",
    )
}

/// Independently certify the actual provider output over the entire rectangle.
/// Matching degrees, knots and unscaled weights are required. The full Bernstein
/// hull bounds homogeneous displacement length and tangent dot products. A strict
/// cross-product projection bound converts tangent residuals into position error;
/// signed triple products exclude a reversed offset and a reversed output surface.
/// Output is accepted only below `default_local - solver` positional error.
pub fn certify_offset(
    source: &BSurface,
    output: &BSurface,
    distance: f64,
    p: &Precision,
) -> Result<()> {
    if !distance.is_finite() {
        return Err(Error::invalid_argument("offset distance must be finite"));
    }
    let _ = cylinder(source, p)?;
    if source.degrees() != output.degrees()
        || source.counts() != output.counts()
        || source.knots() != output.knots()
        || source.weights() != output.weights()
    {
        return Err(unsupported(
            "cylindrical offset provider changed the rational basis",
        ));
    }
    let w = weight_polynomial(source)?;
    let h = homogeneous(source.points(), source.points()[0], &w);
    let out = homogeneous(output.points(), source.points()[0], &w);
    let d: VectorPolynomial = std::array::from_fn(|k| Polynomial {
        u: 2,
        v: 1,
        coefficients: source
            .points()
            .iter()
            .zip(output.points())
            .enumerate()
            .map(|(i, (&a, &b))| {
                (Interval::point(coordinates(b)[k]) - Interval::point(coordinates(a)[k]))
                    * w.coefficients[i]
            })
            .collect(),
    });
    let qu = derivative_numerator(&h, &w, true)?;
    let qv = derivative_numerator(&h, &w, false)?;
    let normal = cross(&qu, &qv)?;
    let out_normal = cross(
        &derivative_numerator(&out, &w, true)?,
        &derivative_numerator(&out, &w, false)?,
    )?;
    positive(
        &dot(&normal, &out_normal)?,
        "cylindrical offset reverses or loses surface regularity",
    )?;
    let (u, v) = source.ranges();
    let mid = source.eval(u.lo + (u.hi - u.lo) * 0.5, v.lo + (v.hi - v.lo) * 0.5);
    // This evaluation only selects a direction. Its suitability over the whole
    // face is proved from original-coordinate interval coefficients below.
    let direction = mid
        .du
        .cross(mid.dv)
        .normalized(p.angular() * mid.du.length() * mid.dv.length())
        .ok_or_else(|| Error::degenerate("cylinder has no reference normal"))?;
    let projection = normal[0]
        .scaled(Interval::point(direction.x))
        .add(&normal[1].scaled(Interval::point(direction.y)))?
        .add(&normal[2].scaled(Interval::point(direction.z)))?;
    positive(
        &projection,
        "whole cylindrical patch regularity is not certified",
    )?;
    let direction_norm = sqrt_upper(
        Interval::point(direction.x).sqr()
            + Interval::point(direction.y).sqr()
            + Interval::point(direction.z).sqr(),
    )?;
    let cross_lower =
        Interval::point(projection.lower()) * reciprocal(Interval::point(direction_norm))?;
    let weight_lower = Interval::point(w.lower());
    let absolute_distance = Interval::point(distance.abs());
    let error = if distance.abs() <= p.solver() {
        Interval::point(norm_upper(&d)?) * reciprocal(weight_lower)? + absolute_distance
    } else {
        positive(
            &dot(&d, &normal)?.scaled(Interval::point(distance.signum())),
            "provider displacement has an uncertified normal direction",
        )?;
        let length_residual =
            dot(&d, &d)?.sub(&w.mul(&w)?.scaled(Interval::point(distance).sqr()))?;
        // ||D/W|-|distance|| <= |D.D-distance^2 W^2|/(Wmin^2 |distance|).
        let length_error = Interval::point(length_residual.absolute_upper())
            * reciprocal(weight_lower.sqr() * absolute_distance)?;
        // Tangent component <= (|D.Qu| |Qv| + |D.Qv| |Qu|)/(Wmin |Qu x Qv|).
        let tangent_error = (Interval::point(dot(&d, &qu)?.absolute_upper())
            * Interval::point(norm_upper(&qv)?)
            + Interval::point(dot(&d, &qv)?.absolute_upper()) * Interval::point(norm_upper(&qu)?))
            * reciprocal(weight_lower * cross_lower)?;
        length_error + tangent_error * Interval::point(2.0)
    };
    if (error - Interval::point(p.default_local() - p.solver())).sign() != Some(Sign::Negative) {
        return Err(unsupported(
            "cylindrical offset exceeds the whole-face normal-distance error budget",
        ));
    }
    Ok(())
}
