//! Exact geometric predicates (Shewchuk's adaptive-precision arithmetic, via `robust`).
//!
//! The returned sign is the sign of the exact determinant of the input coordinates, as if
//! computed with infinite precision. Inputs must be finite.

use robust::{Coord, Coord3D};

use crate::Sign;
use crate::vector::{Point2, Point3};

fn sign_of(value: f64) -> Option<Sign> {
    if value.is_nan() {
        None
    } else if value > 0.0 {
        Some(Sign::Positive)
    } else if value < 0.0 {
        Some(Sign::Negative)
    } else {
        Some(Sign::Zero)
    }
}

/// Orientation of `c` relative to the directed line `a → b`: `Positive` if `a, b, c` turn
/// counter-clockwise, `Zero` if collinear. `None` for non-finite input.
pub fn orient2d(a: Point2, b: Point2, c: Point2) -> Option<Sign> {
    if !(a.x.is_finite()
        && a.y.is_finite()
        && b.x.is_finite()
        && b.y.is_finite()
        && c.x.is_finite()
        && c.y.is_finite())
    {
        return None;
    }
    let coord = |p: Point2| Coord { x: p.x, y: p.y };
    sign_of(robust::orient2d(coord(a), coord(b), coord(c)))
}

/// Orientation of `d` relative to the plane through `a, b, c`: `Positive` if `d` lies on the
/// side from which `a, b, c` appear clockwise (that is, below the plane when `a, b, c` are
/// counter-clockwise seen from above), `Zero` if coplanar. `None` for non-finite input.
pub fn orient3d(a: Point3, b: Point3, c: Point3, d: Point3) -> Option<Sign> {
    if !(a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite()) {
        return None;
    }
    let coord = |p: Point3| Coord3D {
        x: p.x,
        y: p.y,
        z: p.z,
    };
    sign_of(robust::orient3d(coord(a), coord(b), coord(c), coord(d)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient2d_basic() {
        let (a, b) = (Point2::new(0.0, 0.0), Point2::new(1.0, 0.0));
        assert_eq!(orient2d(a, b, Point2::new(0.5, 1.0)), Some(Sign::Positive));
        assert_eq!(orient2d(a, b, Point2::new(0.5, -1.0)), Some(Sign::Negative));
        assert_eq!(orient2d(a, b, Point2::new(7.0, 0.0)), Some(Sign::Zero));
    }

    #[test]
    fn orient2d_resolves_one_ulp_offsets() {
        let a = Point2::new(0.5, 0.5);
        let b = Point2::new(12.0, 12.0);
        let c = Point2::new(24.0, 24.0);
        assert_eq!(orient2d(a, b, c), Some(Sign::Zero));
        let nudged = Point2::new(24.0, 24.0_f64.next_up());
        assert_eq!(orient2d(a, b, nudged), Some(Sign::Positive));
    }

    #[test]
    fn orient3d_sign_convention() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        assert_eq!(
            orient3d(a, b, c, Point3::new(0.0, 0.0, -1.0)),
            Some(Sign::Positive)
        );
        assert_eq!(
            orient3d(a, b, c, Point3::new(0.0, 0.0, 1.0)),
            Some(Sign::Negative)
        );
        assert_eq!(
            orient3d(a, b, c, Point3::new(3.0, 3.0, 0.0)),
            Some(Sign::Zero)
        );
    }

    #[test]
    fn non_finite_input_has_no_sign() {
        let a = Point2::new(f64::NAN, 0.0);
        assert_eq!(
            orient2d(a, Point2::new(1.0, 0.0), Point2::new(0.0, 1.0)),
            None
        );
        for invalid in [f64::NEG_INFINITY, f64::INFINITY, f64::NAN] {
            assert_eq!(
                orient2d(
                    Point2::new(invalid, 0.0),
                    Point2::new(1.0, 0.0),
                    Point2::new(0.0, 1.0)
                ),
                None
            );
            assert_eq!(
                orient3d(
                    Point3::ORIGIN,
                    Point3::new(invalid, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0)
                ),
                None
            );
        }
    }
}
