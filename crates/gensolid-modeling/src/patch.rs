//! Regular NURBS sheets and exact tensor-product subdivision.
use crate::{
    body_check_complete,
    boundary::{self, BoundaryEdge, BoundaryFace},
};
use gensolid_core::{Error, Precision, Result};
use gensolid_geom::{
    BCurve, BSurface, ParamRange, Surface,
    backend::{GeometryBackend, SurfaceConstruction},
    nurbs_bounds::bound_surface,
};
use gensolid_math::{Frame, Interval, Point3, Sign, Vec3};
use gensolid_topo::{BodyId, Model, Sense};

fn dot(a: [Interval; 3], b: Vec3) -> Interval {
    a[0] * Interval::point(b.x) + a[1] * Interval::point(b.y) + a[2] * Interval::point(b.z)
}

/// Sufficient global injectivity certificate: the projected parameter map has a
/// positive definite symmetric Jacobian on the complete rectangular domain.
/// Interval subdivision operates on NURBS equations, never tessellated geometry.
pub(crate) fn certify_graph(s: &BSurface, p: &Precision) -> Result<()> {
    if s.points()
        .iter()
        .any(|v| [v.x, v.y, v.z].iter().any(|x| x.abs() > p.size() * 0.5))
    {
        return Err(Error::invalid_argument(
            "surface control hull outside session size",
        ));
    }
    let (u, v) = s.ranges();
    let mid = s.eval((u.lo + u.hi) * 0.5, (v.lo + v.hi) * 0.5);
    let frame = Frame::new(mid.point, mid.du.cross(mid.dv), mid.du, p)?;
    let (uk, vk) = s.knots();
    for (knots, degree, domain) in [(uk, s.degrees().0, u), (vk, s.degrees().1, v)] {
        if knots
            .windows(degree + 1)
            .any(|w| w[0] > domain.lo && w[0] < domain.hi && w[0].to_bits() == w[degree].to_bits())
        {
            return Err(Error::unimplemented(
                "discontinuous interior surface knots cannot bound a sheet",
            ));
        }
    }
    // This proves the same positive-definite normalized projection Jacobian
    // with tighter rational Bernstein bounds for a supported circular extrusion.
    // Size/continuity guards above and the general interval fallback stay intact.
    if gensolid_geom::cylindrical_offset::certify_cylindrical_graph(s, p).is_ok() {
        return Ok(());
    }
    let mut pending = Vec::new();
    for a in uk
        .windows(2)
        .filter(|a| a[0] < a[1] && a[0] >= u.lo && a[1] <= u.hi)
    {
        for b in vk
            .windows(2)
            .filter(|a| a[0] < a[1] && a[0] >= v.lo && a[1] <= v.hi)
        {
            pending.push((
                ParamRange::bounded(a[0], a[1]),
                ParamRange::bounded(b[0], b[1]),
            ));
        }
    }
    let mut visited = 0;
    while let Some((a, b)) = pending.pop() {
        visited += 1;
        if visited > 4096 {
            return Err(Error::unimplemented(
                "surface graph certification budget exhausted",
            ));
        }
        let bounds = bound_surface(s, a, b)?;
        let j00 = dot(bounds.du, frame.x()) * Interval::point(u.hi - u.lo);
        let j11 = dot(bounds.dv, frame.y()) * Interval::point(v.hi - v.lo);
        let off = (dot(bounds.dv, frame.x()) * Interval::point(v.hi - v.lo)
            + dot(bounds.du, frame.y()) * Interval::point(u.hi - u.lo))
            * Interval::point(0.5);
        let margin = Interval::point(p.linear()).sqr();
        if (j00 - Interval::point(p.linear())).sign() == Some(Sign::Positive)
            && (j11 - Interval::point(p.linear())).sign() == Some(Sign::Positive)
            && (j00 * j11 - off.sqr() - margin).sign() == Some(Sign::Positive)
        {
            continue;
        }
        if j00.sign() == Some(Sign::Negative) || j11.sign() == Some(Sign::Negative) {
            return Err(Error::unimplemented(
                "surface is outside the regular graph family",
            ));
        }
        if (a.hi - a.lo) / (u.hi - u.lo) > (b.hi - b.lo) / (v.hi - v.lo) {
            let m = (a.lo + a.hi) * 0.5;
            if m <= a.lo || m >= a.hi {
                return Err(Error::degenerate("unresolvable U interval"));
            }
            pending.extend([
                (ParamRange::bounded(a.lo, m), b),
                (ParamRange::bounded(m, a.hi), b),
            ]);
        } else {
            let m = (b.lo + b.hi) * 0.5;
            if m <= b.lo || m >= b.hi {
                return Err(Error::degenerate("unresolvable V interval"));
            }
            pending.extend([
                (a, ParamRange::bounded(b.lo, m)),
                (a, ParamRange::bounded(m, b.hi)),
            ]);
        }
    }
    Ok(())
}

pub(crate) fn patch_face(s: &BSurface) -> Result<BoundaryFace> {
    let (u, v) = s.ranges();
    let uv = [
        Point3::new(u.lo, v.lo, 0.0),
        Point3::new(u.hi, v.lo, 0.0),
        Point3::new(u.hi, v.hi, 0.0),
        Point3::new(u.lo, v.hi, 0.0),
    ];
    let curves = [
        s.iso_curve(false, v.lo)?,
        s.iso_curve(true, u.hi)?,
        s.iso_curve(false, v.hi)?.reversed()?,
        s.iso_curve(true, u.lo)?.reversed()?,
    ];
    Ok(BoundaryFace {
        surface: Surface::BSurface(s.clone()),
        sense: Sense::Positive,
        edges: curves
            .into_iter()
            .enumerate()
            .map(|(i, curve)| {
                Ok(BoundaryEdge {
                    curve,
                    uv: Some(boundary::line(uv[i], uv[(i + 1) % 4])?),
                })
            })
            .collect::<Result<_>>()?,
    })
}

/// Create a regular graph sheet with four exact isoparametric boundaries and UV curves.
/// Periodic, folded and singular patches are refused if the graph certificate fails.
pub fn body_create_sheet_surface(model: &mut Model, surface: &BSurface) -> Result<BodyId> {
    certify_graph(surface, model.precision())?;
    let face = patch_face(surface)?;
    model.transaction(|m| {
        let body = boundary::assemble(m, &[face], false)?;
        body_check_complete(m, body)?;
        Ok(body)
    })
}

/// Invoke an injected backend, copy its NURBS output into the kernel, and create a checked sheet.
pub fn body_create_sheet_backend(
    model: &mut Model,
    backend: &dyn GeometryBackend,
    request: &SurfaceConstruction,
) -> Result<BodyId> {
    if backend.precision() != *model.precision() {
        return Err(Error::invalid_argument(
            "backend and model precisions differ",
        ));
    }
    body_create_sheet_surface(model, &backend.construct_surface(request)?)
}

pub(crate) fn split_curve(c: &BCurve, t: f64) -> Result<(BCurve, BCurve)> {
    let r = c.range();
    let n = c.points().len();
    if n != c.degree() + 1
        || c.knots()[..n].iter().any(|x| x.to_bits() != r.lo.to_bits())
        || c.knots()[n..].iter().any(|x| x.to_bits() != r.hi.to_bits())
    {
        return Err(Error::unimplemented(
            "subdivision currently requires a clamped Bezier span",
        ));
    }
    if !t.is_finite() || t <= 0.0 || t >= 1.0 {
        return Err(Error::invalid_argument("split fraction must be interior"));
    }
    let mut work: Vec<[f64; 4]> = c
        .points()
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let w = c.weights().map_or(1.0, |w| w[i]);
            [p.x * w, p.y * w, p.z * w, w]
        })
        .collect();
    let mut left = vec![work[0]];
    let mut right = vec![work[n - 1]];
    for level in 1..n {
        for i in 0..n - level {
            let next = work[i + 1];
            for (value, next) in work[i].iter_mut().zip(next) {
                *value = (1.0 - t) * *value + t * next;
            }
        }
        left.push(work[0]);
        right.push(work[n - level - 1]);
    }
    right.reverse();
    let m = r.lo + (r.hi - r.lo) * t;
    let make = |h: Vec<[f64; 4]>, a: f64, b: f64| {
        BCurve::new(
            c.degree(),
            [vec![a; n], vec![b; n]].concat(),
            h.iter()
                .map(|v| Point3::new(v[0] / v[3], v[1] / v[3], v[2] / v[3]))
                .collect(),
            Some(h.iter().map(|v| v[3]).collect()),
        )
    };
    Ok((make(left, r.lo, m)?, make(right, m, r.hi)?))
}

pub(crate) fn split_surface(
    s: &BSurface,
    constant_u: bool,
    t: f64,
) -> Result<(BSurface, BSurface)> {
    let (nu, nv) = s.counts();
    let (pu, pv) = s.degrees();
    let (uk, vk) = s.knots();
    let count = if constant_u { nv } else { nu };
    let mut left = Vec::new();
    let mut right = Vec::new();
    for k in 0..count {
        let indices: Vec<_> = if constant_u {
            (0..nu).map(|i| i * nv + k).collect()
        } else {
            (0..nv).map(|j| k * nv + j).collect()
        };
        let c = BCurve::new(
            if constant_u { pu } else { pv },
            if constant_u { uk } else { vk }.to_vec(),
            indices.iter().map(|&i| s.points()[i]).collect(),
            s.weights().map(|w| indices.iter().map(|&i| w[i]).collect()),
        )?;
        let (a, b) = split_curve(&c, t)?;
        left.push(a);
        right.push(b);
    }
    let make = |curves: Vec<BCurve>| {
        let mut points = Vec::new();
        let mut weights = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                let (k, index) = if constant_u { (j, i) } else { (i, j) };
                points.push(curves[k].points()[index]);
                weights.push(curves[k].weights().map_or(1.0, |w| w[index]));
            }
        }
        BSurface::new(
            pu,
            pv,
            if constant_u {
                curves[0].knots().to_vec()
            } else {
                uk.to_vec()
            },
            if constant_u {
                vk.to_vec()
            } else {
                curves[0].knots().to_vec()
            },
            nu,
            nv,
            points,
            Some(weights),
        )
    };
    Ok((make(left)?, make(right)?))
}

/// Offset a planar or supported cylindrical NURBS sheet through an injected backend.
/// Independently certify the provider's complete result before creating topology.
/// Planes allow a matching rational basis or a polynomial affine source; circular
/// quadratic extrusions require the same rational basis and a full normal-distance
/// certificate. Other curved offsets and arbitrary reparameterization are refused.
pub fn body_offset_sheet_backend(
    model: &mut Model,
    body: BodyId,
    backend: &dyn GeometryBackend,
    distance: f64,
) -> Result<BodyId> {
    if !distance.is_finite() || backend.precision() != *model.precision() {
        return Err(Error::invalid_argument(
            "invalid offset or backend precision",
        ));
    }
    body_check_complete(model, body)?;
    if model.body(body)?.kind != gensolid_topo::BodyKind::Sheet {
        return Err(Error::invalid_argument("offset input must be a sheet"));
    }
    let faces = model.body_faces(body)?;
    let [id] = faces.as_slice() else {
        return Err(Error::unimplemented(
            "backend sheet offset requires one face",
        ));
    };
    let face = boundary::extract(model, *id)?;
    let Surface::BSurface(source) = face.surface else {
        return Err(Error::unimplemented(
            "backend sheet offset requires a NURBS carrier",
        ));
    };
    let canonical = patch_face(&source)?;
    if face.sense != Sense::Positive || face.edges.len() != 4 {
        return Err(Error::unimplemented(
            "offset requires a natural full parameter rectangle",
        ));
    }
    for (a, b) in face.edges.iter().zip(&canonical.edges) {
        if !boundary::same_curve(&a.curve, &b.curve, model.precision())? {
            return Err(Error::unimplemented("offset of trimmed sheet"));
        }
    }
    let (u, v) = source.ranges();
    let mid = source.eval((u.lo + u.hi) * 0.5, (v.lo + v.hi) * 0.5);
    let p = model.precision();
    let normal = mid
        .du
        .cross(mid.dv)
        .normalized(p.angular())
        .ok_or_else(|| Error::degenerate("singular offset normal"))?;
    let planar = source.points().iter().all(|q| {
        let delta = [
            Interval::point(q.x) - Interval::point(mid.point.x),
            Interval::point(q.y) - Interval::point(mid.point.y),
            Interval::point(q.z) - Interval::point(mid.point.z),
        ];
        let distance = dot(delta, normal);
        distance.lo() >= -p.solver() && distance.hi() <= p.solver()
    });
    if !planar {
        let output = backend.offset_surface(&source, distance)?;
        gensolid_geom::cylindrical_offset::certify_offset(&source, &output, distance, p)?;
        return body_create_sheet_surface(model, &output);
    }
    let output = backend.offset_surface(&source, distance)?;
    let references = offset_reference(&source, &output, p)?;
    for (a, &b) in references.into_iter().zip(output.points()) {
        let delta = [
            (Interval::point(b.x) - a[0]) - Interval::point(normal.x) * Interval::point(distance),
            (Interval::point(b.y) - a[1]) - Interval::point(normal.y) * Interval::point(distance),
            (Interval::point(b.z) - a[2]) - Interval::point(normal.z) * Interval::point(distance),
        ];
        let squared = delta
            .into_iter()
            .map(|v| v.sqr())
            .fold(Interval::point(0.0), |a, b| a + b);
        if (squared - Interval::point(p.default_local() - p.solver()).sqr()).sign()
            != Some(Sign::Negative)
        {
            return Err(Error::unimplemented(
                "offset result is outside the whole-surface error budget",
            ));
        }
    }
    body_create_sheet_surface(model, &output)
}

fn offset_reference(
    source: &BSurface,
    output: &BSurface,
    p: &Precision,
) -> Result<Vec<[Interval; 3]>> {
    let point = |p: Point3| {
        [
            Interval::point(p.x),
            Interval::point(p.y),
            Interval::point(p.z),
        ]
    };
    let weights = |s: &BSurface| {
        (0..s.points().len())
            .map(|i| s.weights().map_or(1.0, |w| w[i]))
            .collect::<Vec<_>>()
    };
    let a = weights(source);
    let b = weights(output);
    if source.degrees() == output.degrees()
        && source.counts() == output.counts()
        && source.knots() == output.knots()
        && a.iter().map(|w| w / a[0]).collect::<Vec<_>>()
            == b.iter().map(|w| w / b[0]).collect::<Vec<_>>()
    {
        return Ok(source.points().iter().copied().map(point).collect());
    }
    // Linear precision of a B-spline basis: affine functions have Greville
    // coefficients. This also certifies degree-elevated approximations of planes.
    if source.degrees() != (1, 1)
        || source.counts() != (2, 2)
        || source.ranges() != output.ranges()
        || !a.windows(2).all(|w| w[0].to_bits() == w[1].to_bits())
        || !b.windows(2).all(|w| w[0].to_bits() == w[1].to_bits())
    {
        return Err(Error::unimplemented(
            "offset backend changed the rational basis",
        ));
    }
    let (u, v) = source.ranges();
    let (uk, vk) = source.knots();
    if uk[0].to_bits() != u.lo.to_bits()
        || uk[3].to_bits() != u.hi.to_bits()
        || vk[0].to_bits() != v.lo.to_bits()
        || vk[3].to_bits() != v.hi.to_bits()
    {
        return Err(Error::unimplemented(
            "affine offset reference requires clamped source",
        ));
    }
    let origin = point(source.points()[0]);
    let x = point(source.points()[2]);
    let y = point(source.points()[1]);
    let end = point(source.points()[3]);
    let residual = (0..3)
        .map(|k| (end[k] - x[k] - y[k] + origin[k]).sqr())
        .fold(Interval::point(0.0), |a, b| a + b);
    if (residual - Interval::point(p.solver()).sqr()).sign() != Some(Sign::Negative) {
        return Err(Error::unimplemented(
            "source parameterization is not affine",
        ));
    }
    let greville =
        |knots: &[f64], degree: usize, count: usize, range: ParamRange| -> Result<Vec<Interval>> {
            let denominator = Interval::point(range.hi) - Interval::point(range.lo);
            if denominator.contains(0.0) {
                return Err(Error::degenerate("offset parameter domain collapsed"));
            }
            let reciprocal = Interval::new(
                (1.0 / denominator.hi()).next_down(),
                (1.0 / denominator.lo()).next_up(),
            )
            .ok_or_else(|| Error::degenerate("offset domain quotient"))?;
            let divisor = Interval::new(
                (1.0 / degree as f64).next_down(),
                (1.0 / degree as f64).next_up(),
            )
            .ok_or_else(|| Error::degenerate("invalid Greville degree"))?;
            Ok((0..count)
                .map(|i| {
                    let average = knots[i + 1..=i + degree]
                        .iter()
                        .map(|&k| Interval::point(k))
                        .fold(Interval::point(0.0), |a, b| a + b)
                        * divisor;
                    (average - Interval::point(range.lo)) * reciprocal
                })
                .collect())
        };
    let (uk, vk) = output.knots();
    let (pu, pv) = output.degrees();
    let (nu, nv) = output.counts();
    let us = greville(uk, pu, nu, u)?;
    let vs = greville(vk, pv, nv, v)?;
    let mut reference = Vec::new();
    for a in us {
        for &b in &vs {
            reference.push(std::array::from_fn(|k| {
                origin[k] + (x[k] - origin[k]) * a + (y[k] - origin[k]) * b
            }));
        }
    }
    Ok(reference)
}
