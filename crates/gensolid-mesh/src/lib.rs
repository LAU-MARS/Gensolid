//! Output-only tessellation. B-rep operations and mass integration do not depend on
//! this crate. Unsupported curvature/trim domains return explicit errors.

mod patch;
mod planar;

pub use patch::surface_tessellate_patch;
pub use planar::body_tessellate_planar;

use gensolid_math::Point3;

/// Triangle output with a certified map from an analytic parameter patch or a
/// planar face's exact projected vertex polygon, including its trim enclosures.
#[derive(Clone, Debug)]
pub struct Tessellation {
    /// Mesh vertices in model coordinates.
    pub vertices: Vec<Point3>,
    /// Counter-clockwise triangles, referring to `vertices`.
    pub triangles: Vec<[usize; 3]>,
    /// Bound for patch interpolation, or planar shadows with bounded line-carrier
    /// endpoint errors, edge tubes and vertex balls. Arbitrary curved trims are
    /// outside the body tessellator's supported domain.
    pub interpolation_bound: f64,
    /// Outward-rounded bound on the vertex evaluation arithmetic.
    pub vertex_roundoff_bound: f64,
}

impl Tessellation {
    /// Guaranteed map deviation including floating-point vertex evaluation enclosures.
    pub fn deviation_bound(&self) -> f64 {
        (self.interpolation_bound + self.vertex_roundoff_bound).next_up()
    }
}
