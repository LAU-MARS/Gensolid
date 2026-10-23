//! Original homogeneous knot insertion with outward-rounded interval controls.
//! Floating-point curve trimming is used only to propose a cylinder, never to
//! certify the original spline or the actual provider output.

use super::*;
use crate::BCurve;

// Exact knot multiplicity and stored-weight equality describe the spline basis,
// not a geometric coincidence decision. Signed zero denotes the same knot.
#[allow(clippy::float_cmp)]
fn same(a: f64, b: f64) -> bool {
    a == b
}

pub(super) struct BezierSpan {
    pub range: [f64; 2],
    pub w: Polynomial,
    pub h: VectorPolynomial,
    pub out: VectorPolynomial,
    pub d: VectorPolynomial,
}

fn shape(surface: &BSurface) -> Result<()> {
    if surface.degrees() != (2, 1) || surface.counts().1 != 2 {
        return Err(unsupported(
            "cylindrical offset requires quadratic U and one linear V span",
        ));
    }
    let (u, v) = surface.ranges();
    let (uk, vk) = surface.knots();
    if uk[..3].iter().any(|&k| !same(k, u.lo))
        || uk[uk.len() - 3..].iter().any(|&k| !same(k, u.hi))
        || vk[..2].iter().any(|&k| !same(k, v.lo))
        || vk[2..].iter().any(|&k| !same(k, v.hi))
        || uk
            .windows(3)
            .any(|k| k[0] > u.lo && k[0] < u.hi && same(k[0], k[2]))
    {
        return Err(unsupported(
            "cylindrical offset requires clamped, continuous spans",
        ));
    }
    let weights = surface
        .weights()
        .ok_or_else(|| unsupported("polynomial curved offsets are not supported"))?;
    if weights.as_chunks::<2>().0.iter().any(|w| !same(w[0], w[1])) {
        return Err(unsupported(
            "cylindrical offset requires repeated weights along the axis",
        ));
    }
    Ok(())
}

/// Extract both surfaces and their directly formed displacement in one common
/// basis. Every interval starts at a raw input f64, including knot differences.
pub(super) fn extract(source: &BSurface, output: Option<&BSurface>) -> Result<Vec<BezierSpan>> {
    shape(source)?;
    let weights = source
        .weights()
        .ok_or_else(|| unsupported("missing rational weights"))?;
    let scale = reciprocal(Interval::point(weights.iter().copied().fold(0.0, f64::max)))?;
    let origin = coordinates(source.points()[0]);
    let mut controls: Vec<[Interval; 10]> = source
        .points()
        .iter()
        .enumerate()
        .map(|(i, &point)| {
            let a = coordinates(point);
            let b = coordinates(output.map_or(point, |surface| surface.points()[i]));
            let w = Interval::point(weights[i]) * scale;
            std::array::from_fn(|j| match j {
                0..=2 => (Interval::point(a[j]) - Interval::point(origin[j])) * w,
                3 => w,
                4..=6 => (Interval::point(b[j - 4]) - Interval::point(origin[j - 4])) * w,
                _ => (Interval::point(b[j - 7]) - Interval::point(a[j - 7])) * w,
            })
        })
        .collect();
    let mut knots = source.knots().0.to_vec();
    let domain = source.ranges().0;
    let simple: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&t| {
            t > domain.lo && t < domain.hi && knots.iter().filter(|&&k| same(k, t)).count() == 1
        })
        .collect();
    for t in simple {
        // Degree two, inserting an existing simple knot exactly once. Convex
        // homogeneous interpolation is enclosed before any f64 rounding can
        // escape: alpha=(t-U_i)/(U_{i+2}-U_i).
        let k = knots.partition_point(|&knot| knot <= t) - 1;
        let count = controls.len() / 2;
        let mut next = vec![[Interval::point(0.0); 10]; controls.len() + 2];
        for i in 0..=k - 2 {
            next[2 * i] = controls[2 * i];
            next[2 * i + 1] = controls[2 * i + 1];
        }
        for i in k - 1..count {
            next[2 * (i + 1)] = controls[2 * i];
            next[2 * (i + 1) + 1] = controls[2 * i + 1];
        }
        let i = k - 1;
        let alpha = (Interval::point(t) - Interval::point(knots[i]))
            * reciprocal(Interval::point(knots[i + 2]) - Interval::point(knots[i]))?;
        for v in 0..2 {
            for j in 0..10 {
                next[2 * i + v][j] = alpha * controls[2 * i + v][j]
                    + (Interval::point(1.0) - alpha) * controls[2 * (i - 1) + v][j];
            }
        }
        knots.insert(k + 1, t);
        controls = next;
    }
    let mut spans = Vec::new();
    for k in 2..controls.len() / 2 {
        if same(knots[k], knots[k + 1]) {
            continue;
        }
        let polynomial = |component: usize| Polynomial {
            u: 2,
            v: 1,
            coefficients: controls[2 * (k - 2)..2 * (k + 1)]
                .iter()
                .map(|c| c[component])
                .collect(),
        };
        let w = polynomial(3);
        positive(
            &w,
            "extracted cylinder denominator is not certified positive",
        )?;
        spans.push(BezierSpan {
            range: [knots[k], knots[k + 1]],
            w,
            h: std::array::from_fn(polynomial),
            out: std::array::from_fn(|i| polynomial(i + 4)),
            d: std::array::from_fn(|i| polynomial(i + 7)),
        });
    }
    Ok(spans)
}

// A numerical proposal only. Unequal endpoint weights are valid (in particular
// after non-midpoint knot refinement), so locate the center by intersecting the
// normals to the two endpoint tangent lines instead of assuming symmetric weights.
pub(super) fn cylinder(surface: &BSurface, p: &Precision) -> Result<Cylinder> {
    // Denominator positivity is a certificate precondition. Check it before
    // any floating-point proposal work (including homogeneous curve trimming).
    let spans = extract(surface, None)?;
    let points = surface.points();
    let direction = points[1] - points[0];
    let axis = direction
        .normalized(p.degenerate_length())
        .ok_or_else(|| Error::degenerate("cylinder extrusion is degenerate"))?;
    let weights = surface
        .weights()
        .ok_or_else(|| unsupported("missing rational weights"))?;
    let profile = BCurve::new(
        2,
        surface.knots().0.to_vec(),
        points.iter().step_by(2).copied().collect(),
        Some(weights.iter().step_by(2).copied().collect()),
    )?;
    let first_end = profile
        .knots()
        .iter()
        .copied()
        .find(|&t| t > profile.range().lo)
        .ok_or_else(|| unsupported("cylinder has no nonempty span"))?;
    let first = profile.trim([profile.range().lo, first_end])?;
    let a = first.points()[0];
    let b = first.points()[1];
    let c = first.points()[2];
    let t0 = b - a;
    let t1 = c - b;
    let inward = t0.cross(t1).cross(t0);
    let denominator = inward.dot(t1);
    if !denominator.is_finite() || denominator == 0.0 {
        return Err(unsupported("circular arc center cannot be proposed"));
    }
    let center = a + inward * ((c - a).dot(t1) / denominator);
    let radius = (a - center).length();
    if !center.is_finite() || !radius.is_finite() || radius <= p.degenerate_length() {
        return Err(Error::degenerate("cylinder radius is degenerate"));
    }
    let mid = surface.eval(
        (surface.ranges().0.lo + first_end) * 0.5,
        (surface.ranges().1.lo + surface.ranges().1.hi) * 0.5,
    );
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
    let candidate = Cylinder {
        center,
        axis,
        radius,
        orientation,
    };
    certify_family(surface, &spans, &candidate, p)?;
    Ok(candidate)
}

fn certify_family(
    source: &BSurface,
    spans: &[BezierSpan],
    c: &Cylinder,
    p: &Precision,
) -> Result<()> {
    let zero = Interval::point(0.0);
    let points = source.points();
    let coordinates_delta = |a: Point3, b: Point3| {
        let a = coordinates(a);
        let b = coordinates(b);
        std::array::from_fn::<_, 3, _>(|i| Interval::point(b[i]) - Interval::point(a[i]))
    };
    let direction = coordinates_delta(points[0], points[1]);
    // The entire source, not merely each piece's own axis, must be a common
    // translation extrusion. Positive repeated V weights make this control bound
    // a full-domain bound on the axial translation residual.
    for pair in points.as_chunks::<2>().0 {
        let delta = coordinates_delta(pair[0], pair[1]);
        let residual = (0..3)
            .map(|i| (delta[i] - direction[i]).sqr())
            .fold(zero, |a, b| a + b);
        if (residual - Interval::point(p.solver()).sqr()).sign() != Some(Sign::Negative) {
            return Err(unsupported(
                "cylinder spans do not share a constant extrusion direction",
            ));
        }
    }
    let axis_norm = Interval::point(c.axis.x).sqr()
        + Interval::point(c.axis.y).sqr()
        + Interval::point(c.axis.z).sqr();
    let inverse_axis_norm = reciprocal(axis_norm)?;
    let origin = coordinates(points[0]);
    let center = coordinates(c.center);
    let axis = [c.axis.x, c.axis.y, c.axis.z];
    for span in spans {
        let relative: VectorPolynomial = [0, 1, 2]
            .map(|i| {
                span.h[i].add(
                    &span
                        .w
                        .scaled(Interval::point(origin[i]) - Interval::point(center[i])),
                )
            })
            .into_iter()
            .collect::<Result<Vec<_>>>()?
            .try_into()
            .map_err(|_| Error::internal("cylinder coordinate count"))?;
        let axial = project(&relative, c.axis)?.scaled(inverse_axis_norm);
        let radial: VectorPolynomial = [0, 1, 2]
            .map(|i| relative[i].sub(&axial.scaled(Interval::point(axis[i]))))
            .into_iter()
            .collect::<Result<Vec<_>>>()?
            .try_into()
            .map_err(|_| Error::internal("cylinder coordinate count"))?;
        let residual = dot(&radial, &radial)?
            .sub(&span.w.mul(&span.w)?.scaled(Interval::point(c.radius).sqr()))?;
        let radial_error = Interval::point(residual.absolute_upper())
            * reciprocal(Interval::point(span.w.lower()).sqr() * Interval::point(c.radius))?;
        if (radial_error - Interval::point(p.solver())).sign() != Some(Sign::Negative) {
            return Err(unsupported(
                "spans are not certified on one common circular cylinder",
            ));
        }
        let normal = cross(
            &derivative_numerator(&span.h, &span.w, true)?,
            &derivative_numerator(&span.h, &span.w, false)?,
        )?;
        positive(
            &dot(&normal, &radial)?.scaled(Interval::point(c.orientation)),
            "cylinder spans do not share a regular radial orientation",
        )?;
    }
    Ok(())
}
