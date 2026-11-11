//! Pure-Rust SMLib-compatible curve/surface backend powered by `solidcore`.
//!
//! Replaces the legacy C++ `usd-brep` FFI bridge and companion WebAssembly transport
//! with safe Rust (`gensolid_geom::solidcore_backend`) while preserving the full
//! `gensolid::smlib` API contract.

pub use gensolid_geom::solidcore_backend::{
    Backend, Curve, CurveEvaluation, Direction, IntersectionBranch, Parameterization,
    ProjectedCurves, SOLIDCORE_VERSION, SolidcoreBackend, Solution, Surface, SurfaceEvaluation,
    UPSTREAM_REVISION,
};

/// Backwards-compatible browser initialization hook.
///
/// Because `solidcore` compiles directly into the Rust WebAssembly module without a
/// separate C++ companion binary, initialization succeeds immediately as a no-op.
#[cfg(all(target_arch = "wasm32", not(target_os = "emscripten")))]
pub fn initialize(_module: wasm_bindgen::JsValue) -> gensolid_core::Result<()> {
    Ok(())
}
