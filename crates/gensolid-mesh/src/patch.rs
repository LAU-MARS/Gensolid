//! Analytic patch interpolation and interval-enclosed vertex evaluation.

use gensolid_core::{Error, ErrorCode, Precision, Result, Tolerance};
use gensolid_geom::{ParamRange, Surface};
use gensolid_math::{Frame, Interval, Point3};

use crate::Tessellation;

fn enclosed_ratio(numerator: f64, denominator: f64) -> Result<Interval> {
    let value = numerator / denominator;
    Interval::new(value.next_down(), value.next_up())
        .ok_or_else(|| Error::degenerate("invalid interval coefficient"))
}

/// Enclose mathematical sine/cosine without assuming platform libm rounding.
/// Reduce by enclosed π/2 and bound the omitted alternating Taylor term.
fn trig(angle: f64) -> Result<(Interval, Interval)> {
    if !angle.is_finite() {
        return Err(Error::invalid_argument("patch angle is not finite"));
    }
    let quotient = (angle / std::f64::consts::FRAC_PI_2).round();
    // Integer range-reduction resource limit; this is not a geometric tolerance.
    if quotient.abs() > (1_u64 << 40) as f64 {
        return Err(Error::unsupported(
            "patch angle exceeds certified range reduction",
        ));
    }
    let q = quotient as i64;
    let half_pi = std::f64::consts::FRAC_PI_2;
    let pi = Interval::new(half_pi.next_down(), half_pi.next_up())
        .ok_or_else(|| Error::internal("invalid pi enclosure"))?;
    let x = Interval::point(angle) - pi * Interval::point(quotient);
    if x.lo().abs().max(x.hi().abs()) > 1.0 {
        return Err(Error::unsupported(
            "angle range reduction cannot be bounded",
        ));
    }
    let x2 = x.sqr();
    let mut sin = x;
    let mut cos = Interval::point(1.0);
    let mut sin_term = x;
    let mut cos_term = Interval::point(1.0);
    for k in 1..=16 {
        sin_term = -(sin_term * x2 * enclosed_ratio(1.0, (2 * k * (2 * k + 1)) as f64)?);
        cos_term = -(cos_term * x2 * enclosed_ratio(1.0, ((2 * k - 1) * 2 * k) as f64)?);
        sin = sin + sin_term;
        cos = cos + cos_term;
    }
    let next_sin = sin_term * x2 * enclosed_ratio(1.0, (34 * 35) as f64)?;
    let next_cos = cos_term * x2 * enclosed_ratio(1.0, (33 * 34) as f64)?;
    let sin_error = next_sin.lo().abs().max(next_sin.hi().abs());
    let cos_error = next_cos.lo().abs().max(next_cos.hi().abs());
    sin = sin
        + Interval::new(-sin_error, sin_error)
            .ok_or_else(|| Error::degenerate("unbounded sine remainder"))?;
    cos = cos
        + Interval::new(-cos_error, cos_error)
            .ok_or_else(|| Error::degenerate("unbounded cosine remainder"))?;
    Ok(match q.rem_euclid(4) {
        0 => (sin, cos),
        1 => (cos, -sin),
        2 => (-sin, -cos),
        _ => (-cos, sin),
    })
}

fn transform(frame: &Frame, local: [Interval; 3]) -> [Interval; 3] {
    let origin = frame.origin();
    let axes = [frame.x(), frame.y(), frame.z()];
    let coordinate = |v: gensolid_math::Vec3, i: usize| match i {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    };
    let origins = [origin.x, origin.y, origin.z];
    std::array::from_fn(|i| {
        Interval::point(origins[i])
            + local[0] * Interval::point(coordinate(axes[0], i))
            + local[1] * Interval::point(coordinate(axes[1], i))
            + local[2] * Interval::point(coordinate(axes[2], i))
    })
}

fn point(surface: &Surface, u: f64, v: f64) -> Result<(Point3, f64)> {
    let (si, co) = trig(u)?;
    let iv = Interval::point(v);
    let (value, enclosure) = match surface {
        Surface::Plane(p) => (
            p.eval(u, v).point,
            transform(p.frame(), [Interval::point(u), iv, Interval::point(0.0)]),
        ),
        Surface::Cylinder(s) => {
            let r = Interval::point(s.radius());
            (
                s.eval(u, v).point,
                transform(s.frame(), [r * co, r * si, iv]),
            )
        }
        Surface::Sphere(s) => {
            let (sv, cv) = trig(v)?;
            let r = Interval::point(s.radius());
            (
                s.eval(u, v).point,
                transform(s.frame(), [r * cv * co, r * cv * si, r * sv]),
            )
        }
        Surface::Torus(s) => {
            let (sv, cv) = trig(v)?;
            let (major, minor) = s.radii();
            let r = Interval::point(minor);
            let rho = Interval::point(major) + r * cv;
            (
                s.eval(u, v).point,
                transform(s.frame(), [rho * co, rho * si, r * sv]),
            )
        }
        Surface::Cone(s) => {
            let (r, t) = s.shape();
            let (st, ct) = trig(t)?;
            let tan_interval = divide(st, ct)?;
            let rho = Interval::point(r) + iv * tan_interval;
            (
                s.eval(u, v).point,
                transform(s.frame(), [rho * co, rho * si, iv]),
            )
        }
        _ => {
            return Err(Error::unimplemented(
                "certified patch tessellation supports analytic carriers only",
            ));
        }
    };
    if !value.is_finite() {
        return Err(Error::degenerate("patch vertex is not finite"));
    }
    let coordinates = [value.x, value.y, value.z];
    let mut squared = Interval::point(0.0);
    for i in 0..3 {
        let delta = enclosure[i] - Interval::point(coordinates[i]);
        let bound = delta.lo().abs().max(delta.hi().abs());
        squared = squared + Interval::point(bound).sqr();
    }
    let error = squared.hi().sqrt().next_up();
    if !error.is_finite() {
        return Err(Error::degenerate("vertex enclosure is not finite"));
    }
    Ok((value, error))
}

fn divide(a: Interval, b: Interval) -> Result<Interval> {
    if b.contains(0.0) {
        return Err(Error::unsupported("analytic denominator includes zero"));
    }
    let x = 1.0 / b.lo();
    let y = 1.0 / b.hi();
    let reciprocal = Interval::new(x.min(y).next_down(), x.max(y).next_up())
        .ok_or_else(|| Error::degenerate("unbounded reciprocal"))?;
    Ok(a * reciprocal)
}

fn hessian(surface: &Surface, v: ParamRange) -> Result<[f64; 3]> {
    let local = match surface {
        Surface::Plane(_) => [0.0, 0.0, 0.0],
        Surface::Cylinder(s) => [s.radius(), 0.0, 0.0],
        Surface::Sphere(s) => [s.radius(), s.radius(), s.radius()],
        Surface::Torus(s) => {
            let (major, minor) = s.radii();
            [(major + minor).next_up(), minor, minor]
        }
        Surface::Cone(s) => {
            let (r, a) = s.shape();
            let (si, co) = trig(a)?;
            let t = divide(si, co)?;
            let radius = Interval::point(r)
                + Interval::new(v.lo, v.hi)
                    .ok_or_else(|| Error::invalid_argument("invalid cone range"))?
                    * t;
            [
                radius.lo().abs().max(radius.hi().abs()),
                0.0,
                t.lo().abs().max(t.hi().abs()),
            ]
        }
        _ => {
            return Err(Error::unimplemented(
                "NURBS/procedural patch curvature bounds are not implemented",
            ));
        }
    };
    let frame = match surface {
        Surface::Plane(s) => s.frame(),
        Surface::Cylinder(s) => s.frame(),
        Surface::Sphere(s) => s.frame(),
        Surface::Torus(s) => s.frame(),
        Surface::Cone(s) => s.frame(),
        _ => return Err(Error::unimplemented("patch frame is unavailable")),
    };
    let mut norm = Interval::point(0.0);
    for axis in [frame.x(), frame.y(), frame.z()] {
        for component in [axis.x, axis.y, axis.z] {
            norm = norm + Interval::point(component).sqr();
        }
    }
    // Frobenius norm bounds the exact affine-frame operator, including small
    // departures from ideal orthonormality in the stored floating frame.
    let scale = norm.hi().sqrt().next_up();
    Ok(local.map(|value| {
        (Interval::point(value) * Interval::point(scale))
            .hi()
            .max(0.0)
    }))
}

fn bound(h: [f64; 3], du: f64, dv: f64) -> f64 {
    let du = Interval::point(du);
    let dv = Interval::point(dv);
    ((Interval::point(h[0]) * du.sqr()
        + Interval::point(h[1]) * dv.sqr()
        + Interval::point(2.0 * h[2]) * du * dv)
        * Interval::point(0.125))
    .hi()
    .max(0.0)
}

/// Tessellate a finite analytic parameter rectangle. For each triangle, Taylor's
/// theorem and bounded parameter variance give `(Muu du² + 2 Muv du dv + Mvv dv²)/8`.
/// Vertex intervals enclose the exact transcendental/affine evaluation as well.
/// This is a surface-patch API; arbitrary trimmed curved faces remain unimplemented.
pub fn surface_tessellate_patch(
    surface: &Surface,
    u: ParamRange,
    v: ParamRange,
    chord: Tolerance,
    precision: &Precision,
) -> Result<Tessellation> {
    let chord = Tolerance::new(chord.value(), precision)?.value();
    if !u.is_bounded() || !v.is_bounded() || u.hi <= u.lo || v.hi <= v.lo {
        return Err(Error::invalid_argument(
            "patch parameter ranges must be finite and increasing",
        ));
    }
    let width_u = (u.hi - u.lo).next_up();
    let width_v = (v.hi - v.lo).next_up();
    if !width_u.is_finite() || !width_v.is_finite() {
        return Err(Error::invalid_argument("patch range width is not finite"));
    }
    let h = hessian(surface, v)?;
    let (mut nu, mut nv) = (1usize, 1usize);
    let budget = chord * 0.5;
    loop {
        let du = (width_u / nu as f64).next_up();
        let dv = (width_v / nv as f64).next_up();
        if bound(h, du, dv) <= budget {
            break;
        }
        let u_score = h[0] * du * du + h[2] * du * dv;
        let v_score = h[1] * dv * dv + h[2] * du * dv;
        // Numeric refinement chooses cells; it does not decide B-rep topology.
        if u_score >= v_score {
            nu = nu.checked_mul(2).ok_or_else(|| {
                Error::new(ErrorCode::CapacityExceeded, "patch subdivision overflow")
            })?;
        } else {
            nv = nv.checked_mul(2).ok_or_else(|| {
                Error::new(ErrorCode::CapacityExceeded, "patch subdivision overflow")
            })?;
        }
        if (nu as u128 + 1) * (nv as u128 + 1) > 1_000_000 {
            return Err(Error::new(
                ErrorCode::CapacityExceeded,
                "certified patch exceeds mesh vertex budget",
            ));
        }
    }
    let us: Vec<_> = (0..=nu)
        .map(|i| {
            if i == 0 {
                u.lo
            } else if i == nu {
                u.hi
            } else {
                u.lo + (u.hi - u.lo) * (i as f64 / nu as f64)
            }
        })
        .collect();
    let vs: Vec<_> = (0..=nv)
        .map(|i| {
            if i == 0 {
                v.lo
            } else if i == nv {
                v.hi
            } else {
                v.lo + (v.hi - v.lo) * (i as f64 / nv as f64)
            }
        })
        .collect();
    let max_u = us
        .windows(2)
        .map(|p| (p[1] - p[0]).next_up())
        .fold(0.0, f64::max);
    let max_v = vs
        .windows(2)
        .map(|p| (p[1] - p[0]).next_up())
        .fold(0.0, f64::max);
    if us.windows(2).chain(vs.windows(2)).any(|p| p[1] <= p[0]) {
        return Err(Error::unsupported(
            "parameter cells collapse under rounding",
        ));
    }
    let mut mesh = Tessellation {
        vertices: Vec::with_capacity((nu + 1) * (nv + 1)),
        triangles: Vec::with_capacity(2 * nu * nv),
        interpolation_bound: bound(h, max_u, max_v),
        vertex_roundoff_bound: 0.0,
    };
    for &u in &us {
        for &v in &vs {
            let (point, error) = point(surface, u, v)?;
            mesh.vertices.push(point);
            mesh.vertex_roundoff_bound = mesh.vertex_roundoff_bound.max(error);
        }
    }
    for i in 0..nu {
        for j in 0..nv {
            let a = i * (nv + 1) + j;
            let b = (i + 1) * (nv + 1) + j;
            mesh.triangles.push([a, b, b + 1]);
            mesh.triangles.push([a, b + 1, a + 1]);
        }
    }
    if mesh.deviation_bound() > chord {
        return Err(Error::unsupported(
            "floating vertex error exceeds the chord budget",
        ));
    }
    Ok(mesh)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use gensolid_geom::{Cylinder, Sphere, Torus};
    use gensolid_math::{Frame, Vec3};
    use proptest::prelude::*;

    #[test]
    fn sine_cosine_enclosures_cover_independent_platform_evaluations() {
        for i in -400..=400 {
            let a = i as f64 / 7.0;
            let (s, c) = trig(a).unwrap();
            assert!(s.contains(a.sin()));
            assert!(c.contains(a.cos()));
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(48))]
        #[test]
        fn curved_patch_barycentric_residual_is_below_the_certified_bound(
            r in 0.1f64..10.0,ox in -20.0f64..20.0,oy in -20.0f64..20.0,oz in -20.0f64..20.0,
            ax in -1.0f64..1.0,ay in -1.0f64..1.0,az in 0.1f64..1.0,kind in 0usize..3,
        ) {
            let p=Precision::DEFAULT;
            let f=Frame::from_axis(Point3::new(ox,oy,oz),Vec3::new(ax,ay,az),&p).unwrap();
            let surface=match kind {0=>Surface::Cylinder(Cylinder::new(f,r,&p).unwrap()),1=>Surface::Sphere(Sphere::new(f,r,&p).unwrap()),_=>Surface::Torus(Torus::new(f,3.0*r,r,&p).unwrap())};
            let u=ParamRange::bounded(0.0,1.0);let v=ParamRange::bounded(0.0,0.75);
            let chord=Tolerance::new(p.linear()*100_000.0,&p).unwrap();
            let mesh=surface_tessellate_patch(&surface,u,v,chord,&p).unwrap();
            prop_assert!(mesh.deviation_bound()<=chord.value());
            let cells=mesh.triangles.len()/2;
            // Use each actual mesh triangle's vertex projection to recover its
            // surface parameters, then test affine points against independent eval.
            let mut store=gensolid_geom::GeomStore::new();let id=store.add_surface(surface).unwrap();
            for triangle in mesh.triangles.iter().step_by((mesh.triangles.len()/17).max(1)) {
                let a=mesh.vertices[triangle[0]];let b=mesh.vertices[triangle[1]];let c=mesh.vertices[triangle[2]];
                let linear=a+(b-a)*0.25+(c-a)*0.375;
                let distance=store.distance_to_surface(id,linear,&p).unwrap();
                prop_assert!(distance<=mesh.deviation_bound(),"{distance} > {}",mesh.deviation_bound());
            }
            prop_assert!(cells>0);
        }
    }
}
