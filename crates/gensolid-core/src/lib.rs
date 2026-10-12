//! Core types shared by every Gensolid crate.
//!
//! - [`arena`]: generational arenas with typed identifiers, used for all model entities.
//! - [`error`]: the kernel-wide error type with stable numeric codes for the C API.
//! - [`tolerance`]: session precision and local (per edge / per vertex) tolerances.

pub mod arena;
pub mod error;
pub mod tolerance;

pub use arena::{Arena, Id};
pub use error::{Error, ErrorCode, Result};
pub use tolerance::{Precision, Tolerance, effective_tolerance};
