//! Curves.

use std::f64::consts::TAU;

use gensolid_core::{Error, Precision, Result};
use gensolid_math::{Frame, Point3, Vec3};

use crate::nurbs::BCurve;
use crate::param::ParamRange;

/// Point and first derivative of a curve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurveDerivs {
    /// Position.
    pub point: Point3,
    /// First derivative with respect to the parameter.
    pub d1: Vec3,
}

/// Infinite line `origin + t * direction`, with unit `direction`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line {
    origin: Point3,
    direction: Vec3,
}

impl Line {
    /// Creates a line; `direction` is normalised.
    pub fn new(origin: Point3, direction: Vec3, precision: &Precision) -> Result<Self> {
        if !origin.is_finite() {
            return Err(Error::invalid_argument("line origin is not finite"));
        }
        let direction = direction
            .normalized(precision.degenerate_length())
            .ok_or_else(|| Error::degenerate("line direction has zero length"))?;
        Ok(Self { origin, direction })
    }

    /// Line through two distinct points, parameterised by arc length from `a`.
    pub fn through(a: Point3, b: Point3, precision: &Precision) -> Result<Self> {
        Self::new(a, b - a, precision)
    }

    /// Point at `t = 0`.
    pub const fn origin(&self) -> Point3 {
        self.origin
    }

    /// Unit direction.
    pub const fn direction(&self) -> Vec3 {
        self.direction
    }

    /// Point and derivative at `t`.
    pub fn eval(&self, t: f64) -> CurveDerivs {
        CurveDerivs {
            point: self.origin + self.direction * t,
            d1: self.direction,
        }
    }

    /// Parameter of the closest point to `p`.
    pub fn closest_param(&self, p: Point3) -> f64 {
        (p - self.origin).dot(self.direction)
    }
}

/// Circle in the `x`-`y` plane of `frame`, centred on its origin; `t` is the angle from `x`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle {
    frame: Frame,
    radius: f64,
}

impl Circle {
    /// Creates a circle; the radius must exceed linear precision.
    pub fn new(frame: Frame, radius: f64, precision: &Precision) -> Result<Self> {
        if !(radius.is_finite() && radius > precision.linear()) {
            return Err(Error::degenerate(
                "circle radius must exceed linear precision",
            ));
        }
        Ok(Self { frame, radius })
    }

    /// Frame: origin is the centre, `z` the axis.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Radius.
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// Point and derivative at angle `t`.
    pub fn eval(&self, t: f64) -> CurveDerivs {
        let (s, c) = t.sin_cos();
        let (x, y) = (self.frame.x(), self.frame.y());
        CurveDerivs {
            point: self.frame.origin() + (x * c + y * s) * self.radius,
            d1: (y * c - x * s) * self.radius,
        }
    }

    /// Angle in `[0, 2π)` of the closest point to `p` (0 for points on the axis).
    pub fn closest_param(&self, p: Point3) -> f64 {
        let local = self.frame.to_local(p);
        if local.x == 0.0 && local.y == 0.0 {
            return 0.0;
        }
        local.y.atan2(local.x).rem_euclid(TAU)
    }
}

/// Ellipse in the `x`-`y` plane of `frame`; the major axis lies along `x`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ellipse {
    frame: Frame,
    major: f64,
    minor: f64,
}

impl Ellipse {
    /// Creates an ellipse with `major >= minor > linear precision`.
    pub fn new(frame: Frame, major: f64, minor: f64, precision: &Precision) -> Result<Self> {
        if !(minor.is_finite() && minor > precision.linear()) {
            return Err(Error::degenerate(
                "ellipse minor radius must exceed linear precision",
            ));
        }
        if !(major.is_finite() && major >= minor) {
            return Err(Error::invalid_argument(
                "ellipse major radius must be at least the minor",
            ));
        }
        Ok(Self {
            frame,
            major,
            minor,
        })
    }

    /// Frame: origin is the centre, `x` the major axis.
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// Major and minor radii.
    pub const fn radii(&self) -> (f64, f64) {
        (self.major, self.minor)
    }

    /// Point and derivative at `t`.
    pub fn eval(&self, t: f64) -> CurveDerivs {
        let (s, c) = t.sin_cos();
        let (x, y) = (self.frame.x(), self.frame.y());
        CurveDerivs {
            point: self.frame.origin() + x * (self.major * c) + y * (self.minor * s),
            d1: y * (self.minor * c) - x * (self.major * s),
        }
    }
}

/// Every curve class of the kernel.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Curve {
    /// Infinite line.
    Line(Line),
    /// Full circle.
    Circle(Circle),
    /// Full ellipse.
    Ellipse(Ellipse),
    /// NURBS curve.
    BCurve(BCurve),
}

impl Curve {
    /// Parameter range.
    pub fn range(&self) -> ParamRange {
        match self {
            Self::Line(_) => ParamRange::UNBOUNDED,
            Self::Circle(_) | Self::Ellipse(_) => ParamRange::FULL_TURN,
            Self::BCurve(c) => c.range(),
        }
    }

    /// Point and first derivative at `t`.
    pub fn eval(&self, t: f64) -> CurveDerivs {
        match self {
            Self::Line(c) => c.eval(t),
            Self::Circle(c) => c.eval(t),
            Self::Ellipse(c) => c.eval(t),
            Self::BCurve(c) => c.eval(t),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use super::*;

    #[test]
    fn circle_eval_and_inverse() {
        let p = Precision::DEFAULT;
        let frame = Frame::new(Point3::new(1.0, 0.0, 0.0), Vec3::Z, Vec3::X, &p).unwrap();
        let c = Circle::new(frame, 2.0, &p).unwrap();
        let d = c.eval(FRAC_PI_2);
        assert!(d.point.distance(Point3::new(1.0, 2.0, 0.0)) <= p.linear());
        assert!((d.d1 - Vec3::new(-2.0, 0.0, 0.0)).length() <= p.linear());
        assert!((c.closest_param(Point3::new(1.0, 5.0, 3.0)) - FRAC_PI_2).abs() <= p.angular());
    }

    #[test]
    fn rejects_degenerate_line() {
        let p = Precision::DEFAULT;
        assert!(Line::through(Point3::ORIGIN, Point3::ORIGIN, &p).is_err());
    }
}
