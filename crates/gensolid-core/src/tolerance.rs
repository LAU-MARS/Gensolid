//! Session precision and local tolerances.
//!
//! This is the only kernel module allowed to contain literal tolerance values; the rule is
//! enforced by `verifier/policy`. Everything else derives thresholds from [`Precision`] or
//! from the [`Tolerance`] stored on an edge or vertex.

use crate::error::{Error, Result};

/// Session-wide precision.
///
/// - `linear`: two points closer than this are coincident (model units).
/// - `angular`: two directions closer than this are parallel (radians).
/// - `size`: all geometry must lie inside a cube of this edge length centred on the origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Precision {
    linear: f64,
    angular: f64,
    size: f64,
}

impl Precision {
    /// Linear precision 1e-8, angular precision 1e-11 rad, size box 1000 model units.
    pub const DEFAULT: Self = Self {
        linear: 1.0e-8,
        angular: 1.0e-11,
        size: 1.0e3,
    };

    /// Creates a precision, rejecting non-finite, non-positive or inconsistent values.
    pub fn new(linear: f64, angular: f64, size: f64) -> Result<Self> {
        let positive = |v: f64| v.is_finite() && v > 0.0;
        if !positive(linear) || !positive(angular) || !positive(size) {
            return Err(Error::invalid_argument(
                "precision values must be finite and positive",
            ));
        }
        if linear >= size {
            return Err(Error::invalid_argument(
                "linear precision must be smaller than the size box",
            ));
        }
        Ok(Self {
            linear,
            angular,
            size,
        })
    }

    /// Linear precision.
    pub const fn linear(&self) -> f64 {
        self.linear
    }

    /// Angular precision in radians.
    pub const fn angular(&self) -> f64 {
        self.angular
    }

    /// Edge length of the size box.
    pub const fn size(&self) -> f64 {
        self.size
    }

    /// Tolerance of an edge or vertex that carries no local tolerance: half the linear
    /// precision, so two exact entities within tolerance of each other are coincident.
    pub const fn default_local(&self) -> f64 {
        self.linear * 0.5
    }

    /// Lengths at or below this are degenerate.
    pub const fn degenerate_length(&self) -> f64 {
        self.linear
    }

    /// Step size at which iterative solvers stop: two orders below linear precision.
    pub const fn solver(&self) -> f64 {
        self.linear * 1.0e-2
    }
}

impl Default for Precision {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Local tolerance attached to an edge or vertex (tolerant modelling).
///
/// An entity with a local tolerance is a tube (edge) or ball (vertex) of that radius rather
/// than an exact curve or point. Tolerances only grow during modelling.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Tolerance(f64);

impl Tolerance {
    /// Creates a tolerance. It must be finite, at least the default local tolerance, and
    /// smaller than the size box.
    pub fn new(value: f64, precision: &Precision) -> Result<Self> {
        if !value.is_finite() || value < precision.default_local() || value >= precision.size() {
            return Err(Error::invalid_argument(format!(
                "tolerance {value} outside [{}, {})",
                precision.default_local(),
                precision.size()
            )));
        }
        Ok(Self(value))
    }

    /// Radius of the tolerance tube or ball.
    pub const fn value(self) -> f64 {
        self.0
    }

    /// The larger of two tolerances.
    pub const fn max(self, other: Self) -> Self {
        if other.0 > self.0 { other } else { self }
    }

    /// Grows the tolerance by `amount` (accumulated error from an operation).
    pub fn widened(self, amount: f64) -> Result<Self> {
        if !amount.is_finite() || amount < 0.0 {
            return Err(Error::invalid_argument(
                "tolerance growth must be finite and non-negative",
            ));
        }
        Ok(Self(self.0 + amount))
    }
}

/// Tolerance in effect for an entity: its local tolerance, or the session default.
pub fn effective_tolerance(local: Option<Tolerance>, precision: &Precision) -> f64 {
    local.map_or(precision.default_local(), Tolerance::value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_inconsistent_precision() {
        assert!(Precision::new(1.0, 1.0e-3, 0.5).is_err());
        assert!(Precision::new(f64::NAN, 1.0e-3, 10.0).is_err());
        assert!(Precision::new(1.0e-6, 1.0e-9, 100.0).is_ok());
    }

    #[test]
    #[allow(clippy::float_cmp)] // The values are passed through unchanged, so exactly equal.
    fn local_tolerance_never_below_default() {
        let p = Precision::DEFAULT;
        assert!(Tolerance::new(p.default_local() * 0.5, &p).is_err());
        let t = Tolerance::new(1.0e-5, &p).unwrap();
        assert_eq!(effective_tolerance(Some(t), &p), 1.0e-5);
        assert_eq!(effective_tolerance(None, &p), p.default_local());
    }
}
