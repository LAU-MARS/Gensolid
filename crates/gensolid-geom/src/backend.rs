//! Owned geometry boundary for optional numerical algorithm providers.
//! The native Rust provider lives in this crate; foreign providers live above it.
//! Topology never owns foreign pointers.
use crate::{BCurve, BSurface};
use gensolid_core::{Precision, Result};
use gensolid_math::{Point3, Vec3};

/// Surface construction request; all input geometry is owned by the kernel.
#[derive(Clone, Debug)]
pub enum SurfaceConstruction {
    /// Skin through section curves in the supplied order.
    Loft(Vec<BCurve>),
    /// Translate a profile by a displacement.
    Extrude(BCurve, Vec3),
    /// Sweep a profile along a path using the provider's framing convention.
    Sweep(BCurve, BCurve),
    /// Keep profile orientation fixed: `S(u,v) = profile(u) + path(v) - path(v_min)`.
    /// This explicit translation convention also applies to curved spatial paths.
    TranslationSweep {
        /// Profile, kept in its original orientation at the path start.
        profile: BCurve,
        /// Path supplying displacement, without rotating the profile.
        path: BCurve,
    },
    /// Rotate a profile about an axis through an origin, by radians.
    Revolve {
        /// Profile to revolve.
        profile: BCurve,
        /// Point on the axis.
        origin: Point3,
        /// Axis direction.
        axis: Vec3,
        /// Rotation in radians.
        angle: f64,
    },
}

/// A numerical intersection candidate, not a completeness certificate.
#[derive(Clone, Debug)]
pub struct SurfaceIntersectionCandidate {
    /// World-space curve.
    pub curve: BCurve,
    /// UV curve on the first surface.
    pub first_uv: BCurve,
    /// UV curve on the second surface.
    pub second_uv: BCurve,
    /// Provider diagnostic only; never used as a certified bound.
    pub reported_deviation: f64,
}

/// Optional geometry algorithms. Results must be validated by the modelling caller.
/// No Send/Sync requirement: providers may own thread-confined native sessions.
pub trait GeometryBackend {
    /// Precision used by the provider.
    fn precision(&self) -> Precision;
    /// Construct a surface without constructing topology.
    fn construct_surface(&self, request: &SurfaceConstruction) -> Result<BSurface>;
    /// Candidate curves with their paired UV representations.
    fn intersect_surfaces(
        &self,
        first: &BSurface,
        second: &BSurface,
    ) -> Result<Vec<SurfaceIntersectionCandidate>>;
    /// Approximate normal offset; callers must verify the intended distance.
    fn offset_surface(&self, surface: &BSurface, distance: f64) -> Result<BSurface>;
}
