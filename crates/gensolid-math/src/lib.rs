//! Numerical foundations of the Gensolid kernel.
//!
//! Topological decisions (which side, coincident or not) must never be taken on a raw
//! floating-point comparison. They go through [`predicates`] (exact) or [`interval`]
//! (certified enclosure) and fall back to tolerances only when both are inconclusive.

pub mod bbox;
pub mod frame;
pub mod interval;
pub mod predicates;
pub mod vector;

pub use bbox::Aabb;
pub use frame::Frame;
pub use interval::Interval;
pub use predicates::{orient2d, orient3d};
pub use vector::{Point2, Point3, Vec2, Vec3};

/// Sign of an exactly or certifiably evaluated quantity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Sign {
    /// Strictly negative.
    Negative,
    /// Exactly zero.
    Zero,
    /// Strictly positive.
    Positive,
}
