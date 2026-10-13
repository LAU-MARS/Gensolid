# SMLib curve and surface backend

This branch wraps NVIDIA's [USD BRep / SMLib](https://github.com/NVIDIA-Omniverse/usd-brep)
through a private C ABI and a safe Rust API in `gensolid::smlib`. The pinned revision is
`c9c5979f04cea8c38213004d4f04df0dcfb40d05` (Apache-2.0). This is an optional native
geometry dependency, not a Rust translation of the upstream kernel.

## Build and run

Default Gensolid builds remain Rust-only. Enabling `smlib` builds the unmodified
SMLib C++ core, three SM_API implementation files, and our adapter. It does not
build or link OpenUSD, Python, OCCT, the USD importers, or the OCCT bridges.
The builder uses the standard allocator and serial execution; TBB is not needed.

Requirements: Rust/Cargo, Python 3, Git, a C++17 compiler, and `ar`. The builder
supports native macOS and Linux targets. Windows, cross-compilation and WASM
with this feature are not supported. The existing Rust-only WASM build remains
available.

```sh
git clone https://github.com/NVIDIA-Omniverse/usd-brep.git target/usd-brep
git -C target/usd-brep checkout c9c5979f04cea8c38213004d4f04df0dcfb40d05
export SMLIB_SOURCE_DIR="$PWD/target/usd-brep"
# Optional shared build cache to avoid compiling the kernel again per Cargo profile:
export SMLIB_BUILD_DIR="$PWD/target/smlib-native-cache"
cargo run -p gensolid-ffi --features smlib --example smlib_curves_surfaces
cargo test -p gensolid-ffi --features smlib
cargo clippy -p gensolid-ffi --all-targets --features smlib -- -D warnings
```

Cargo never downloads native source automatically. The builder verifies the
revision and rejects modifications to the two upstream source directories.
It rebuilds native files when source/compiler/flags change and copies the upstream
license and third-party notices beside the native archive. `CXX` can select a
compiler executable. The source and build directories are machine-local and are
not committed. Use an absolute path for `SMLIB_BUILD_DIR`.

```rust,no_run
use gensolid::smlib::Backend;
use gensolid_core::Precision;
use gensolid_math::Point3;

# fn main() -> gensolid_core::Result<()> {
let backend = Backend::new(Precision::DEFAULT)?;
let curve = backend.circle(Point3::new(0.0, 0.0, 0.0), 2.0)?;
let [start, end] = curve.domain()?;
let sample = curve.evaluate((start + end) / 2.0)?;
let native_nurbs = curve.nurbs()?;
# let _ = (sample, native_nurbs);
# Ok(())
# }
```

## Coverage

The agreed scope is **curve/surface geometry**, not all public classes in the
roughly 534,000-line SMLib core. No solid boolean, fillet, face-loop editing,
healing, file import/export, tessellation or body API is exposed by this adapter.

| Area | Rust interface |
|---|---|
| NURBS interchange | `import_curve`, `import_surface`, `Curve::nurbs`, `Surface::nurbs`; degrees, expanded knots, Euclidean points and rational weights |
| Evaluation | natural domains, curve position/first/second derivatives, surface position and all first/second partials, curve continuity |
| Point projection | closest-point parameter solutions on curves and surfaces |
| Intersections | curve/curve and curve/surface solutions, including overlap intervals; surface/surface curves, both UV curves, upstream type and reported deviation |
| Curve construction | line segment, circle, ellipse, arc, helix, control-point curve, interpolation, approximation, rectangle and regular polygon |
| Surface construction | corner patch, control grid, interpolated grid, unordered-point fit, extrusion, sweep, ruled surface, revolution, loft/skin |
| Projection and offsets | curve drop, trimmed projection, UV lift, explicit-plane or detected-plane curve offset, surface offset |
| Curve editing | compatibility, ordering, knot simplification, line/arc decomposition, trim, reverse, degree elevation, knot insertion |
| Surface editing | rectangular domain trim, U/V reversal, isoparametric curve, degree elevation, knot insertion, U/V split |

Curves constructed in XY can also be supplied in arbitrary frames through NURBS
interchange. Angles in the Rust API are radians. Native NURBS control nets map to
Gensolid's `u_index * v_count + v_index` layout. UV curves carry `u,v` in their
`x,y` coordinates. Constructors requiring NURBS reject unsupported native classes
with an error rather than reinterpret pointers. Rational weights must be positive,
as required by the existing Gensolid NURBS model. Imports reject upstream knot
normalization that changes the supplied representation.

`Solution<D>` preserves upstream solution kinds and start/end parameter tuples;
single-valued solutions repeat their start at the end. Unknown kinds remain
observable. A surface intersection branch contains a world-space curve and paired
parameter-space curves. These objects are **numerical results**, not an
`IntersectionSet` completeness certificate. Empty success is not a certified empty
intersection, and reported deviation is not a rigorous bound. Coincidence or
singular cases may return an upstream error. Consumers must not infer a valid
B-rep from a successful SMLib status alone.

Rectangular trimming changes the surface parameter domain. It is not arbitrary
loop trimming and does not create a Gensolid face. No new modelling operation is
registered in `Session`, so this feature does not bypass `check_body`. Integrating
these curves into Gensolid topology remains a separate modelling task, subject
to the existing verifier and property-test requirements.

## Ownership, errors and concurrency

`Backend`, `Curve`, and `Surface` share a native session arena using `Rc`.
A curve/surface remains usable after its `Backend` variable is dropped. Native
objects are freed when the final handle to that session disappears; long-running
batch applications should use one session per bounded job. Individual handle
clones do not allocate native geometry. Mixing objects from different sessions
is rejected before crossing the ABI. Handles are neither `Send` nor `Sync`.

Every native entry point is serialized by one mutex, including creation and
teardown. This is necessary because SM_API also uses a process-global default
context. Independent parallel jobs should use separate processes. That upstream
default context is initialized once and retained until process exit.

The private ABI exchanges only pointers, fixed-width status/op fields, `size_t`
counts, and primitive arrays. Native exceptions are caught at the boundary;
errors carry a Gensolid `ErrorCode` and the original upstream numeric status.
Error results do not become empty success. Native result buffers are copied and
freed by their allocating side. Editing operations work on copies so upstream
partial modification on failure cannot change existing Rust handles.

All Rust `unsafe` remains in `gensolid-ffi`. The dependency direction stays
`ffi -> api -> modeling -> geom`: geometry crates do not depend on FFI. Native
backend calls are currently an explicit API for callers; there is no hidden
fallback in Gensolid's exact geometry dispatcher.

## Tolerances and upstream limits

The adapter supplies `Precision::linear()` / `angular()` to low-level surface
intersection, curve intersections, point projection, sweep, offsets, projection,
lift, knot simplification, line/arc decomposition and iso-curve construction.
It adds no literal tolerance constants. Upstream internal numerical thresholds
and legacy context behavior still exist: passing a tolerance does not prove
Gensolid's predicate or certification requirements.

Control-point and interpolation construction use the corresponding upstream
policies. Approximation from point clouds retains the upstream fitter's policy;
it does not promise a session-precision fit. Helix fitting has upstream input and
accuracy limits; unsupported requests fail. `order_curves` does not reverse curves
or certify a correctly oriented closed loop. The upstream [known issues](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/KNOWN_ISSUES.md)
include seam/pole projection defects. Those constraints are not removed by FFI.

## Provenance and distribution

The adapter is original glue against upstream public headers. No OCCT source was
read, translated, adapted or linked. The external checkout includes OCCT-related
folders, but they are excluded from this build. Preserve NVIDIA's Apache-2.0
license and applicable copyright/notice information when distributing a binary
with the feature enabled. See the project `NOTICE` and the license copies staged
beside `libgensolid_smlib.a`.

## Validation on this branch

Validated on Apple Silicon macOS with Rust 1.99.0:

- `cargo fmt --all` and the full workspace Clippy check with warnings denied.
- `cargo test` for the existing default workspace, including repository policy checks.
- The full workspace Clippy check with `gensolid-ffi/smlib` enabled.
- Eleven native integration/property tests: rational data round-trips, unequal
  surface dimensions, derivatives, random coordinate frames and weights, near
  degeneracy, projection/lift, overlap intervals, paired UV intersection curves,
  constructors, edits, lifetime management and invalid inputs. The two property
  tests each run 24 cases. All call the real C++ library, not a mock.
- C++ adapter syntax/warnings checked with `-Wall -Wextra -Werror` and upstream
  headers treated as system headers.
- Regenerated `gensolid.h` with cbindgen; the public C header is unchanged.
- The runnable example returns one intersection branch and zero residual at its
  midpoint. This is a smoke test, not a completeness proof.

Linux is supported by the build driver but has not been executed in this session.
