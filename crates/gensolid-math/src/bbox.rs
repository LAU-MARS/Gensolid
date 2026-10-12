//! Axis-aligned bounding boxes.

use crate::vector::{Point3, Vec3};

/// Axis-aligned box with `min <= max` componentwise.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb {
    min: Point3,
    max: Point3,
}

impl Aabb {
    /// Smallest box containing all points, or `None` for an empty iterator.
    pub fn from_points(points: impl IntoIterator<Item = Point3>) -> Option<Self> {
        let mut points = points.into_iter();
        let first = points.next()?;
        Some(points.fold(
            Self {
                min: first,
                max: first,
            },
            |b, p| b.including(p),
        ))
    }

    /// Lower corner.
    pub const fn min(&self) -> Point3 {
        self.min
    }

    /// Upper corner.
    pub const fn max(&self) -> Point3 {
        self.max
    }

    /// The box grown to contain `p`.
    pub fn including(self, p: Point3) -> Self {
        Self {
            min: Point3::new(
                self.min.x.min(p.x),
                self.min.y.min(p.y),
                self.min.z.min(p.z),
            ),
            max: Point3::new(
                self.max.x.max(p.x),
                self.max.y.max(p.y),
                self.max.z.max(p.z),
            ),
        }
    }

    /// Smallest box containing both boxes.
    pub fn union(self, other: Self) -> Self {
        self.including(other.min).including(other.max)
    }

    /// The box grown by `amount` on every side.
    pub fn expanded(self, amount: f64) -> Self {
        let d = Vec3::new(amount, amount, amount);
        Self {
            min: self.min - d,
            max: self.max + d,
        }
    }

    /// Whether `p` lies inside or on the box.
    pub fn contains(&self, p: Point3) -> bool {
        (self.min.x..=self.max.x).contains(&p.x)
            && (self.min.y..=self.max.y).contains(&p.y)
            && (self.min.z..=self.max.z).contains(&p.z)
    }

    /// Diagonal vector from `min` to `max`.
    pub fn diagonal(&self) -> Vec3 {
        self.max - self.min
    }
}
