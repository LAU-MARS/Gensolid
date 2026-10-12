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
pub mod model;

pub use block::make_block;
pub use entity::*;
pub use model::{Model, RawEditor};
