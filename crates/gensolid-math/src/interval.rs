//! Interval arithmetic with outward rounding.
//!
//! Every IEEE-754 operation in round-to-nearest mode is within half an ulp of the exact
//! result, so widening each computed bound by one ulp (`next_down` / `next_up`) encloses
//! the exact result. The enclosures are slightly looser than with directed rounding, but
//! portable and free of global FPU state.

use std::ops::{Add, Mul, Neg, Sub};

use crate::Sign;

/// Closed interval `[lo, hi]` guaranteed to contain the exact real value it represents.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    /// The whole real line; the result of any operation that produced NaN.
    pub const ENTIRE: Self = Self {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    /// Creates `[lo, hi]`, or `None` if `lo > hi` or either bound is NaN.
    pub fn new(lo: f64, hi: f64) -> Option<Self> {
        (lo <= hi).then_some(Self { lo, hi })
    }

    /// Degenerate interval containing exactly `value`. NaN maps to [`Interval::ENTIRE`].
    pub fn point(value: f64) -> Self {
        if value.is_nan() {
            Self::ENTIRE
        } else {
            Self {
                lo: value,
                hi: value,
            }
        }
    }

    /// Lower bound.
    pub const fn lo(self) -> f64 {
        self.lo
    }

    /// Upper bound.
    pub const fn hi(self) -> f64 {
        self.hi
    }

    /// Width, rounded up.
    pub fn width(self) -> f64 {
        (self.hi - self.lo).next_up()
    }

    /// Whether `value` lies inside the interval.
    pub fn contains(self, value: f64) -> bool {
        self.lo <= value && value <= self.hi
    }

    /// Smallest interval containing both operands.
    pub fn hull(self, other: Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.max(other.hi),
        }
    }

    /// Certified sign, or `None` when the interval straddles zero.
    pub fn sign(self) -> Option<Sign> {
        if self.lo > 0.0 {
            Some(Sign::Positive)
        } else if self.hi < 0.0 {
            Some(Sign::Negative)
        } else if self.lo == 0.0 && self.hi == 0.0 {
            Some(Sign::Zero)
        } else {
            None
        }
    }

    /// Square, which unlike `x * x` is never negative.
    pub fn sqr(self) -> Self {
        let a = self.lo * self.lo;
        let b = self.hi * self.hi;
        let hi = a.max(b).next_up();
        let lo = if self.contains(0.0) {
            0.0
        } else {
            a.min(b).next_down().max(0.0)
        };
        Self::checked(lo, hi)
    }

    fn checked(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() {
            Self::ENTIRE
        } else {
            Self { lo, hi }
        }
    }
}

impl Add for Interval {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::checked((self.lo + o.lo).next_down(), (self.hi + o.hi).next_up())
    }
}

impl Sub for Interval {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::checked((self.lo - o.hi).next_down(), (self.hi - o.lo).next_up())
    }
}

impl Mul for Interval {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let products = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        if products.iter().any(|p| p.is_nan()) {
            return Self::ENTIRE;
        }
        let lo = products.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::checked(lo.next_down(), hi.next_up())
    }
}

impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    #[test]
    fn one_tenth_sum_is_enclosed() {
        let tenth = Interval::point(0.1);
        let sum = (0..10).fold(Interval::point(0.0), |acc, _| acc + tenth);
        assert!(sum.contains(1.0));
        assert!(sum.lo() < sum.hi());
    }

    #[test]
    fn sign_is_certified_or_absent() {
        assert_eq!(Interval::point(2.0).sign(), Some(Sign::Positive));
        assert_eq!(Interval::new(-1.0, 1.0).and_then(Interval::sign), None);
        assert_eq!(Interval::point(0.0).sign(), Some(Sign::Zero));
    }

    proptest! {
        #[test]
        fn operations_enclose_point_results(a in -1.0e6f64..1.0e6, b in -1.0e6f64..1.0e6) {
            let (ia, ib) = (Interval::point(a), Interval::point(b));
            prop_assert!((ia + ib).contains(a + b));
            prop_assert!((ia - ib).contains(a - b));
            prop_assert!((ia * ib).contains(a * b));
            prop_assert!(ia.sqr().contains(a * a));
            prop_assert!(ia.sqr().lo() >= 0.0);
        }
    }
}
