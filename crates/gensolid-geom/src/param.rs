//! Parameter ranges of curves and surfaces.

use std::f64::consts::TAU;

/// Parameter range in one direction. Bounds may be infinite (lines, planes).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParamRange {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
    /// Whether the geometry repeats with period `hi - lo`.
    pub periodic: bool,
}

impl ParamRange {
    /// `(-∞, ∞)`.
    pub const UNBOUNDED: Self = Self {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
        periodic: false,
    };

    /// `[0, 2π)`, periodic.
    pub const FULL_TURN: Self = Self {
        lo: 0.0,
        hi: TAU,
        periodic: true,
    };

    /// Bounded, non-periodic range.
    pub const fn bounded(lo: f64, hi: f64) -> Self {
        Self {
            lo,
            hi,
            periodic: false,
        }
    }

    /// Whether both bounds are finite.
    pub fn is_bounded(&self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }

    /// Period length, if periodic.
    pub fn period(&self) -> Option<f64> {
        self.periodic.then_some(self.hi - self.lo)
    }

    /// Maps `t` into the range: wraps periodic ranges, clamps the others.
    pub fn normalize(&self, t: f64) -> f64 {
        match self.period() {
            Some(period) if period > 0.0 => self.lo + (t - self.lo).rem_euclid(period),
            _ => t.clamp(self.lo, self.hi),
        }
    }
}
