//! Arithmetic expressions retain outward-rounded enclosures until a branch is chosen.

use std::ops::{Add, Div, Mul, Neg, Sub};

use gensolid_core::{Error, Result};
use gensolid_math::{Interval, Point3, Sign, Vec3};

#[derive(Clone, Copy)]
pub(super) struct Scalar {
    pub value: f64,
    bounds: Interval,
}

impl Scalar {
    pub fn new(value: f64) -> Self {
        Self {
            value,
            bounds: Interval::point(value),
        }
    }

    pub fn sqr(self) -> Self {
        Self {
            value: self.value * self.value,
            bounds: self.bounds.sqr(),
        }
    }

    pub fn sqrt(self) -> Self {
        let bounds = match Interval::new(
            self.bounds.lo().max(0.0).sqrt().next_down().max(0.0),
            self.bounds.hi().max(0.0).sqrt().next_up(),
        ) {
            Some(bounds) => bounds,
            None => Interval::ENTIRE,
        };
        Self {
            value: self.value.max(0.0).sqrt(),
            bounds,
        }
    }

    /// Resolve exact/interval signs first. Only an enclosure containing zero may
    /// use the explicit tolerance fallback; callers supply the expression's units.
    pub fn sign(self, tolerance: f64) -> Result<Sign> {
        if !self.value.is_finite() {
            return Err(Error::unsupported("intersection arithmetic is not finite"));
        }
        if let Some(sign) = self.bounds.sign() {
            return Ok(sign);
        }
        // Explicit tolerance fallback for an uncertified sign, never for a
        // certified nonzero discriminant or distance.
        if self.value.abs() <= tolerance {
            Ok(Sign::Zero)
        } else {
            Err(Error::unsupported(
                "intersection sign is numerically unresolved",
            ))
        }
    }

    /// A domain boundary must be certified unless it was separately established
    /// as a polynomial endpoint root by the endpoint residual classification. Never snap an interior root into a
    /// span merely because its parameter lies near a boundary.
    pub fn strict_sign(self) -> Result<Sign> {
        if !self.value.is_finite() {
            return Err(Error::unsupported("intersection arithmetic is not finite"));
        }
        self.bounds
            .sign()
            .ok_or_else(|| Error::unsupported("intersection root domain is numerically unresolved"))
    }
}

macro_rules! arithmetic {
    ($trait:ident, $method:ident, $operator:tt) => {
        impl $trait for Scalar {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                Self { value: self.value $operator rhs.value, bounds: self.bounds $operator rhs.bounds }
            }
        }
    };
}
arithmetic!(Add, add, +);
arithmetic!(Sub, sub, -);
arithmetic!(Mul, mul, *);

impl Div for Scalar {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let reciprocal = if rhs.bounds.contains(0.0) {
            Interval::ENTIRE
        } else {
            match Interval::new(
                (1.0 / rhs.bounds.hi()).next_down(),
                (1.0 / rhs.bounds.lo()).next_up(),
            ) {
                Some(bounds) => bounds,
                None => Interval::ENTIRE,
            }
        };
        Self {
            value: self.value / rhs.value,
            bounds: self.bounds * reciprocal,
        }
    }
}

impl Neg for Scalar {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            bounds: -self.bounds,
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct Vector {
    pub x: Scalar,
    pub y: Scalar,
    pub z: Scalar,
}

impl Vector {
    pub fn new(value: Vec3) -> Self {
        Self {
            x: Scalar::new(value.x),
            y: Scalar::new(value.y),
            z: Scalar::new(value.z),
        }
    }

    pub fn between(a: Point3, b: Point3) -> Self {
        Self::new(a.to_vec()) - Self::new(b.to_vec())
    }

    pub fn dot(self, other: Self) -> Scalar {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn norm2(self) -> Scalar {
        self.x.sqr() + self.y.sqr() + self.z.sqr()
    }

    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn value(self) -> Vec3 {
        Vec3::new(self.x.value, self.y.value, self.z.value)
    }

    pub fn zero(self, tolerance: f64) -> Result<bool> {
        Ok(self.x.sign(tolerance)? == Sign::Zero
            && self.y.sign(tolerance)? == Sign::Zero
            && self.z.sign(tolerance)? == Sign::Zero)
    }
}

impl Sub for Vector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<Scalar> for Vector {
    type Output = Self;
    fn mul(self, rhs: Scalar) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}
