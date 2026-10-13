//! Orthonormal frames (the "axis2" / basis set of commercial kernels).

use gensolid_core::{Error, Precision, Result};

use crate::vector::{Point3, Vec3};

/// Right-handed orthonormal frame: an origin, a principal axis `z` and a reference
/// direction `x`; `y = z × x`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    origin: Point3,
    x: Vec3,
    y: Vec3,
    z: Vec3,
}

impl Frame {
    /// The world frame.
    pub const WORLD: Self = Self {
        origin: Point3::ORIGIN,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };

    /// Builds a frame whose `z` is `axis` and whose `x` is `ref_direction` projected onto
    /// the plane perpendicular to `axis`.
    pub fn new(
        origin: Point3,
        axis: Vec3,
        ref_direction: Vec3,
        precision: &Precision,
    ) -> Result<Self> {
        if !origin.is_finite() {
            return Err(Error::invalid_argument("frame origin is not finite"));
        }
        let z = axis
            .normalized(precision.degenerate_length())
            .ok_or_else(|| Error::degenerate("frame axis has zero length"))?;
        let ref_length = ref_direction.length();
        if !(ref_length.is_finite() && ref_length > precision.degenerate_length()) {
            return Err(Error::degenerate(
                "frame reference direction has zero length",
            ));
        }
        // Cross products avoid cancellation in `ref - z * dot(ref, z)`, whose
        // rounded residual can retain a measurable component along the axis.
        let y_raw = z.cross(ref_direction);
        let y_length = y_raw.length();
        if y_length <= precision.angular() * ref_length {
            return Err(Error::degenerate(
                "frame reference direction is parallel to the axis",
            ));
        }
        let y = y_raw / y_length;
        let x = y.cross(z).normalized(precision.angular()).ok_or_else(|| {
            Error::degenerate("frame reference direction cannot be orthogonalized")
        })?;
        Ok(Self { origin, x, y, z })
    }

    /// Builds a frame with the given axis and an arbitrary perpendicular reference direction.
    pub fn from_axis(origin: Point3, axis: Vec3, precision: &Precision) -> Result<Self> {
        let candidates = [Vec3::X, Vec3::Y, Vec3::Z];
        let least_aligned = candidates
            .into_iter()
            .min_by(|a, b| a.dot(axis).abs().total_cmp(&b.dot(axis).abs()))
            .unwrap_or(Vec3::X);
        Self::new(origin, axis, least_aligned, precision)
    }

    /// Origin.
    pub const fn origin(&self) -> Point3 {
        self.origin
    }

    /// Reference direction.
    pub const fn x(&self) -> Vec3 {
        self.x
    }

    /// Second in-plane direction, `z × x`.
    pub const fn y(&self) -> Vec3 {
        self.y
    }

    /// Principal axis.
    pub const fn z(&self) -> Vec3 {
        self.z
    }

    /// World position of local coordinates `local`.
    pub fn to_world(&self, local: Vec3) -> Point3 {
        self.origin + self.vector_to_world(local)
    }

    /// World direction of a local vector.
    pub fn vector_to_world(&self, local: Vec3) -> Vec3 {
        self.x * local.x + self.y * local.y + self.z * local.z
    }

    /// Local coordinates of a world point.
    pub fn to_local(&self, p: Point3) -> Vec3 {
        let d = p - self.origin;
        Vec3::new(d.dot(self.x), d.dot(self.y), d.dot(self.z))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projects_reference_direction() {
        let p = Precision::DEFAULT;
        let f = Frame::new(Point3::ORIGIN, Vec3::Z, Vec3::new(2.0, 0.0, 5.0), &p).unwrap();
        assert_eq!(f.x(), Vec3::X);
        assert_eq!(f.y(), Vec3::Y);
    }

    #[test]
    fn rejects_parallel_reference_direction() {
        let p = Precision::DEFAULT;
        assert!(Frame::new(Point3::ORIGIN, Vec3::Z, Vec3::new(0.0, 0.0, 3.0), &p).is_err());
    }

    #[test]
    fn local_world_round_trip() {
        let p = Precision::DEFAULT;
        let f = Frame::from_axis(Point3::new(1.0, 2.0, 3.0), Vec3::new(1.0, 1.0, 0.0), &p).unwrap();
        let q = Point3::new(-4.0, 0.5, 7.0);
        let back = f.to_world(f.to_local(q));
        assert!(back.distance(q) <= p.linear());
    }
}
