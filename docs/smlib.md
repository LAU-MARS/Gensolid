# Pure-Rust `solidcore` curve and surface backend

> Phase 1 update: the external C++ `usd-brep` FFI bridge (`native/smlib.cpp`), Python C++
> builder (`tools/build-smlib.py`), process-global session mutex, and companion Emscripten
> WebAssembly transport have been completely replaced by pure safe Rust via `solidcore`
> (`gensolid_geom::solidcore_backend::SolidcoreBackend`, re-exported as `gensolid::smlib::Backend`).

Default Gensolid builds remain lightweight and self-contained. Enabling the `solidcore`
feature on `gensolid-geom` / `gensolid-api` (or the backwards-compatible `smlib` feature
on `gensolid-ffi` / `gensolid-wasm`) compiles the pure-Rust `solidcore` 1.0.1 crate
(Apache-2.0, ported from NVIDIA's `usd-brep` SMLib/SM_API) directly into the Rust workspace
with `#![forbid(unsafe_code)]`.

## Build and run (Windows, Linux, macOS, and `wasm32-unknown-unknown`)

No C++17 compiler, `SMLIB_SOURCE_DIR` checkout, `ar`, or Emscripten SDK is needed on any
platform:

```sh
cargo test -p gensolid-geom --features solidcore
cargo test -p gensolid-api --test solidcore_backend
cargo test -p gensolid-ffi --features smlib
cargo run -p gensolid-ffi --features smlib --example smlib_curves_surfaces
cargo run -p gensolid-ffi --features smlib --example smlib_modeling
cargo build -p gensolid-wasm --target wasm32-unknown-unknown --release --features smlib
```

## Architecture and certification boundary

1. **Zero `unsafe` and isolated sessions**: Every `SolidcoreBackend` owns an isolated
   `solidcore::Kernel` context in safe Rust without any process-wide mutex or raw pointer.
2. **Candidate + Certificate contract**: Curves, surfaces, intersections, and offsets
   returned by `SolidcoreBackend` enter Gensolid through `gensolid_geom::backend::GeometryBackend`.
   Every resulting B-rep sheet, prism, cylindrical/planar offset, or oblique cut is still
   independently proved by Gensolid's interval/Bernstein certificates and checked by the
   immutable `gensolid-check` verifier before any entity or tag is published.
3. **Provenance**: `solidcore` is licensed under Apache-2.0 (derived from NVIDIA's
   Apache-2.0 `usd-brep` SMLib/SM_API without OCCT or USD/Hydra modules). Copyright and
   attribution notices are recorded in [`NOTICE`](../NOTICE). The full reference workspace
   for subsequent M1–M5 algorithm ports lives locally under `baseline/solidcore-revo`.
