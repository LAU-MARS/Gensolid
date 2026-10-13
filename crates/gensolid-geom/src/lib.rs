//! Geometry of the Gensolid kernel.
//!
//! The geometry classes mirror those of the Parasolid XT format so models round-trip
//! without approximation:
//!
//! | Kind       | Curves                          | Surfaces                                   |
//! |------------|---------------------------------|--------------------------------------------|
//! | Analytic   | line, circle, ellipse           | plane, cylinder, cone, sphere, torus       |
//! | Free-form  | B-curve (NURBS)                 | B-surface (NURBS)                          |
//! | Procedural | *planned:* SP, trimmed, intersection | swept, spun; *planned:* offset, blend |
//!
//! Procedural geometry refers to other geometry by id, so all geometry lives in a
//! [`GeomStore`] and is evaluated through it.

pub mod curve;
pub mod intersect;
pub mod nurbs;
pub mod param;
mod project;
pub mod store;
pub mod surface;

pub use curve::{Circle, Curve, CurveDerivs, Ellipse, Line};
pub use intersect::{
    IntersectionCase, IntersectionScore, IntersectionSet, analytic_intersection_score,
    intersect_curve_surface, intersect_surfaces, intersection_score,
};
pub use nurbs::{BCurve, BSurface};
pub use param::ParamRange;
pub use store::{
    CurveId, CurveProjection, GeomStore, ProjectionKind, SurfaceId, SurfaceProjection,
};
pub use surface::{Cone, Cylinder, Plane, Sphere, Spun, Surface, SurfaceDerivs, Swept, Torus};

#[cfg(test)]
mod robustness_tests;
