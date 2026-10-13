//! Outward-rounded homogeneous de Boor bounds. No polygonization is involved.
use crate::{BCurve, BSurface, ParamRange};
use gensolid_core::{Error, Result};
use gensolid_math::{Interval, Point3};

/// Enclosures for position and first derivatives on a single knot rectangle.
#[derive(Clone, Copy, Debug)]
pub struct NurbsBounds {
    /// Cartesian position enclosure.
    pub point: [Interval; 3],
    /// First parameter derivative enclosure.
    pub du: [Interval; 3],
    /// Second parameter derivative enclosure; zero for curves.
    pub dv: [Interval; 3],
}

pub(crate) fn divide(a: Interval, b: Interval) -> Result<Interval> {
    if b.contains(0.0) {
        return Err(Error::unimplemented("interval denominator contains zero"));
    }
    let x = 1.0 / b.lo();
    let y = 1.0 / b.hi();
    Ok(a * Interval::new(x.min(y).next_down(), x.max(y).next_up())
        .ok_or_else(|| Error::degenerate("unbounded quotient"))?)
}
#[derive(Clone, Copy)]
struct Jet {
    v: Interval,
    d: [Interval; 2],
}
impl Jet {
    fn constant(v: Interval) -> Self {
        Self {
            v,
            d: [Interval::point(0.0); 2],
        }
    }
    fn sub(self, b: Self) -> Self {
        Self {
            v: self.v - b.v,
            d: std::array::from_fn(|i| self.d[i] - b.d[i]),
        }
    }
    fn add(self, b: Self) -> Self {
        Self {
            v: self.v + b.v,
            d: std::array::from_fn(|i| self.d[i] + b.d[i]),
        }
    }
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            d: std::array::from_fn(|i| self.d[i] * b.v + self.v * b.d[i]),
        }
    }
    fn scale(self, b: Interval) -> Self {
        Self {
            v: self.v * b,
            d: self.d.map(|x| x * b),
        }
    }
}
fn span(knots: &[f64], degree: usize, count: usize, r: ParamRange) -> Result<usize> {
    if !r.lo.is_finite() || !r.hi.is_finite() || r.lo > r.hi {
        return Err(Error::invalid_argument("invalid interval parameter domain"));
    }
    (degree..count)
        .find(|&i| knots[i] < knots[i + 1] && r.lo >= knots[i] && r.hi <= knots[i + 1])
        .ok_or_else(|| Error::invalid_argument("bound domain must lie in one nonempty knot span"))
}
fn deboor(
    mut work: Vec<[Jet; 4]>,
    knots: &[f64],
    degree: usize,
    span: usize,
    r: ParamRange,
    axis: usize,
) -> Result<[Jet; 4]> {
    let mut t = Jet::constant(
        Interval::new(r.lo, r.hi)
            .ok_or_else(|| Error::invalid_argument("reversed parameter domain"))?,
    );
    t.d[axis] = Interval::point(1.0);
    let one = Jet::constant(Interval::point(1.0));
    for level in 1..=degree {
        for j in (level..=degree).rev() {
            let left = Interval::point(knots[span - degree + j]);
            let denominator = Interval::point(knots[span + j - level + 1]) - left;
            let alpha = t
                .sub(Jet::constant(left))
                .scale(divide(Interval::point(1.0), denominator)?);
            work[j] = std::array::from_fn(|k| {
                one.sub(alpha)
                    .mul(work[j - 1][k])
                    .add(alpha.mul(work[j][k]))
            });
        }
    }
    work.get(degree)
        .copied()
        .ok_or_else(|| Error::internal("empty de Boor net"))
}
fn control(p: Point3, reference: Point3, w: f64, scale: f64) -> Result<[Jet; 4]> {
    let weight = divide(Interval::point(w), Interval::point(scale))?;
    let p = [p.x, p.y, p.z];
    let o = [reference.x, reference.y, reference.z];
    Ok(std::array::from_fn(|k| {
        Jet::constant(if k == 3 {
            weight
        } else {
            (Interval::point(p[k]) - Interval::point(o[k])) * weight
        })
    }))
}
fn cartesian(h: [Jet; 4], reference: Point3, weight_bound: Interval) -> Result<NurbsBounds> {
    let w = Interval::new(
        h[3].v.lo().max(weight_bound.lo()),
        h[3].v.hi().min(weight_bound.hi()),
    )
    .ok_or_else(|| Error::internal("homogeneous weight enclosure is inconsistent"))?;
    let mut point = [Interval::point(0.0); 3];
    let mut derivatives = [[Interval::point(0.0); 3]; 2];
    let origin = [reference.x, reference.y, reference.z];
    for k in 0..3 {
        let q = divide(h[k].v, w)?;
        point[k] = q + Interval::point(origin[k]);
        for (axis, d) in derivatives.iter_mut().enumerate() {
            d[k] = divide(h[k].d[axis] - q * h[3].d[axis], w)?;
        }
    }
    Ok(NurbsBounds {
        point,
        du: derivatives[0],
        dv: derivatives[1],
    })
}
fn weights(values: impl Iterator<Item = f64>) -> Result<(f64, Interval)> {
    let values: Vec<_> = values.collect();
    let max = values.iter().copied().fold(0.0, f64::max);
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let lo = divide(Interval::point(min), Interval::point(max))?.lo();
    if lo <= 0.0 {
        return Err(Error::unimplemented(
            "weight dynamic range exceeds interval arithmetic",
        ));
    }
    Ok((
        max,
        Interval::new(lo, 1.0).ok_or_else(|| Error::degenerate("invalid weight range"))?,
    ))
}
/// Nonempty knot spans, in parameter order.
pub(crate) fn knot_spans(knots: &[f64], degree: usize, count: usize) -> Vec<ParamRange> {
    (degree..count)
        .filter(|&i| knots[i] < knots[i + 1])
        .map(|i| ParamRange::bounded(knots[i], knots[i + 1]))
        .collect()
}
/// Bound an arbitrary-degree positive-weight NURBS curve on one knot span.
pub fn bound_curve(curve: &BCurve, r: ParamRange) -> Result<NurbsBounds> {
    let p = curve.degree();
    let s = span(curve.knots(), p, curve.points().len(), r)?;
    let weight = |i: usize| curve.weights().map_or(1.0, |w| w[i]);
    let (scale, wb) = weights((s - p..=s).map(weight))?;
    let reference = curve.points()[s - p];
    let net = (s - p..=s)
        .map(|i| control(curve.points()[i], reference, weight(i), scale))
        .collect::<Result<Vec<_>>>()?;
    cartesian(deboor(net, curve.knots(), p, s, r, 0)?, reference, wb)
}
/// Bound an arbitrary-degree positive-weight NURBS surface on one knot rectangle.
pub fn bound_surface(surface: &BSurface, u: ParamRange, v: ParamRange) -> Result<NurbsBounds> {
    let (p, q) = surface.degrees();
    let (nu, nv) = surface.counts();
    let (uk, vk) = surface.knots();
    let us = span(uk, p, nu, u)?;
    let vs = span(vk, q, nv, v)?;
    let ids = (us - p..=us)
        .flat_map(|i| (vs - q..=vs).map(move |j| i * nv + j))
        .collect::<Vec<_>>();
    let weight = |i: usize| surface.weights().map_or(1.0, |w| w[i]);
    let (scale, wb) = weights(ids.iter().map(|&i| weight(i)))?;
    let reference = surface.points()[(us - p) * nv + vs - q];
    let mut net = Vec::new();
    for i in us - p..=us {
        let row = (vs - q..=vs)
            .map(|j| {
                let k = i * nv + j;
                control(surface.points()[k], reference, weight(k), scale)
            })
            .collect::<Result<Vec<_>>>()?;
        net.push(deboor(row, vk, q, vs, v, 1)?);
    }
    cartesian(deboor(net, uk, p, us, u, 0)?, reference, wb)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    proptest! {
        #[test]
        fn rational_quartic_interval_contains_position_and_derivative(t in 0.1f64..0.8, width in 0.001f64..0.1, offset in -100.0f64..100.0) {
            let c=BCurve::new(4,[vec![0.0;5],vec![1.0;5]].concat(),(0..5).map(|i|Point3::new(offset+i as f64,(i*i) as f64, (i%2) as f64)).collect(),Some(vec![1.0,2.0,0.5,3.0,1.0])).unwrap();
            let bounds=bound_curve(&c,ParamRange::bounded(t,t+width)).unwrap();
            for a in [t,t+width*0.3,t+width] {let e=c.eval(a);for(k,(x,d))in [(e.point.x,e.d1.x),(e.point.y,e.d1.y),(e.point.z,e.d1.z)].into_iter().enumerate(){prop_assert!(bounds.point[k].contains(x));prop_assert!(bounds.du[k].contains(d));}}
        }
    }
}
