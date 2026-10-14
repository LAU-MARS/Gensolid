//! Boundary-representation topology.
//!
//! The entity hierarchy follows Parasolid's so that XT files map one-to-one:
//!
//! ```text
//! Body ─┬─ Region ── Shell ─┬─ Face ── Loop ── Fin ── Edge ── Vertex
//!       │                   ├─ wire Edge
//!       │                   └─ acorn Vertex
//!       └─ regions[0] is always the infinite void region
//! ```
//!
//! - Non-manifold from the start: an edge has any number of fins, a face separates a
//!   front and a back shell, a body may hold several solid and void regions.
//! - Tolerant modelling: edges and vertices may carry a local [`Tolerance`]; without one
//!   they are exact at session precision.
//!
//! [`Tolerance`]: gensolid_core::Tolerance

mod block;
pub mod entity;
mod euler;
mod extrusion;
mod inquire;
pub mod model;
mod parametric;
mod polyhedron;

pub use block::make_block;
pub use entity::*;
pub use euler::{EdgeSplit, edge_parameter_bounds, split_line_edge};
pub use extrusion::{make_homothetic_loft, make_linear_sweep, make_prism};
pub use inquire::{body_bounding_box, body_volume};
pub use model::{FinGeometry, Model, RawEditor};
pub use parametric::{LoftSection, ParametricBody, make_loft_sections, make_revolution};
pub use polyhedron::make_convex_polyhedron;
