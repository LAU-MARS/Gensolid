# SMLib curve and surface backend

> Migration update (2026-10-23): the supported core modelling paths also have a native
> `RustGeometryBackend`, without the build/runtime requirements below. This document
> still describes the optional original C++ backend. See [the staged Rust migration](rust-migration.md)
> for the exact capability boundary; the full SMLib facade has not been replaced.


This branch wraps NVIDIA's [USD BRep / SMLib](https://github.com/NVIDIA-Omniverse/usd-brep)
through a private C ABI and a safe Rust API in `gensolid::smlib`. The pinned revision is
`c9c5979f04cea8c38213004d4f04df0dcfb40d05` (Apache-2.0). This is an optional C++
geometry dependency, not a Rust translation of the upstream kernel.

## Build and run

Default Gensolid builds remain Rust-only. Enabling `smlib` builds the unmodified
SMLib C++ core, three SM_API implementation files, and our adapter. It does not
build or link OpenUSD, Python, OCCT, the USD importers, or the OCCT bridges.
The builder uses the standard allocator and serial execution; TBB is not needed.

Requirements: Rust/Cargo, Python 3, Git, a C++17 compiler, and `ar`. The builder
supports native macOS and Linux targets. Native Windows and general native
cross-compilation are not yet implemented. Browser Wasm uses the separate package
build below; the existing Rust-only Wasm build also remains available.

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

## Browser Wasm release

For new SMLib browser releases, use the [single-Wasm build](smlib-single-wasm.md).
It links Rust and C++ into one module with a shared memory. The dual-module route
below remains available for existing wasm-bindgen consumers and comparison.

### Dual-module build

The **complete Rust curve/surface facade** now works on `wasm32-unknown-unknown`
with `smlib` enabled. The deployable package has **two Wasm modules**, accompanying
JavaScript, TypeScript declarations, notices, and a SHA-256 manifest. It is not a
single statically linked Wasm file. The same package runs in supported browsers
on macOS, Linux and Windows; end users do not install or compile C++.

This preserves the existing `wasm-bindgen` target. Rust invokes the same SMLib
operation dispatcher through a synchronous JS transport; the entire selected C++
core is compiled by Emscripten. No algorithm is replaced by a JS approximation.
Inputs and results are copied between separate memories, which adds call/copy
overhead. Batch geometry work in Rust when practical. C++ exceptions are caught
at the C ABI, with Emscripten's JS exception runtime enabled. No WASI host,
filesystem, pthreads, SharedArrayBuffer, or cross-origin isolation is required.
Separate workers must initialize their own package and own sessions.

Validated tool versions: Emscripten **6.0.12**, Rust **1.99.0**, and wasm-bindgen CLI
**0.2.129**, matching Cargo.lock. Install/activate that Emscripten SDK and source
`emsdk_env.sh` (Windows: `emsdk_env.bat`). The SDK build host must have Python 3,
Git and Rust; build-host coverage is distinct from browser runtime portability.

```sh
rustup target add wasm32-unknown-unknown
cargo install wasm-bindgen-cli --version 0.2.129 --locked
export SMLIB_SOURCE_DIR="$PWD/target/usd-brep" # pinned checkout described above
python3 tools/package-smlib-wasm.py --output target/smlib-web
```

`cargo build -p gensolid-wasm --target wasm32-unknown-unknown --release --features smlib`
builds the Rust module only. Use the packaging command for a complete release;
it compiles the companion, runs wasm-bindgen, and stages all runtime files and
licenses. Serve the **whole output directory** over HTTP(S), including `snippets/`.
No network/CDN request is required by the package beyond fetching its own assets.

```js
import load, { SmlibBackend } from './smlib-web/index.mjs';
await load();
const backend = new SmlibBackend();
const curve = backend.circle([0, 0, 0], 2);
backend.free(); // curve retains the session
console.log(curve.evaluate(curve.domain()[0])); // XYZ, d1, d2
curve.free();
```

The JS convenience classes currently expose circle/line/corner-patch construction,
evaluation, domains and curve reversal. **All previously wrapped curve/surface
operations are available to Rust code** through `gensolid::smlib` (also re-exported
by `gensolid-wasm`), and can be used in application-specific wasm-bindgen exports.
The browser regression suite exercises that full Rust facade, including
intersections, projection, offsets and editing. This does not add solid modelling
operations or change the existing body API.

For custom Rust Wasm applications, enable the `gensolid-ffi/smlib` feature and call
`gensolid::smlib::initialize` with the initialized companion module before making a
`Backend`. The provided JS loader handles initialization once per JS realm.
Replacing a companion while handles exist is rejected. Each cross-module result
is copied before the C++ allocation is freed; JS views are reacquired after memory
growth. wasm32 result offsets are compile-time asserted in the C++ adapter.

The SMLib module reserves an **8 MiB stack** with overflow checks. Emscripten's
small default stack was insufficient for the tested intersection paths and caused
memory corruption; this setting is required. Linear memory grows on demand and
remains bounded by Emscripten's wasm32 limits. Call `.free()` on exported JS handles
when done. Releasing arenas returns allocations for reuse but does not shrink a
Wasm memory. Large or adversarial geometry can still exceed finite memory/stack
limits and trap; the regression cases are not an arbitrary-input memory bound.

### Repeatable browser validation

```sh
python3 tools/package-smlib-wasm.py --validation --output target/smlib-web-validation
npm install --prefix tools/wasm
npx --prefix tools/wasm playwright install chromium firefox webkit
node tools/wasm/test-browser.mjs target/smlib-web-validation
BROWSER=firefox node tools/wasm/test-browser.mjs target/smlib-web-validation
BROWSER=webkit node tools/wasm/test-browser.mjs target/smlib-web-validation
```

`BROWSER_EXECUTABLE` selects an existing compatible browser executable. Alternatively,
serve that output directory and open `validation.html`. Production packages exclude
the validation export/page. For C++ AddressSanitizer use `--validation --sanitize`
and a different output directory; this diagnostic build uses substantially more
memory and should not be shipped. The validation feature reuses nine existing
native tests without weakening their assertions, plus 48 seeded random samples
covering rational surfaces, frames, derivatives and near-degenerate inputs. It
repeats the suite 12 times (684 checks), forces memory growth, tests initialization
and live-handle replacement errors, retains handles after backend release, and
checks the existing Rust body API. An independent module Worker runs another 57
checks through the production loader. Native proptest tests remain enabled unchanged.

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
loop trimming and does not itself create a Gensolid face. The separate checked
modelling layer now consumes this backend through borrowed `GeometryBackend`
calls in `Session`, constructing graph sheets and directional prisms, validating
restricted bounded face cuts, and certifying planar sheet offsets. See the
[modelling pipeline support matrix](modeling-pipeline.md) and run the
`smlib_modeling` example. General curved Booleans and shelling remain unsupported.

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

Browser validation on this branch (same artifacts, macOS host):

- Edge / Chromium 155.0.4283.45, Firefox 153.0 and WebKit 26.5: 684 main-thread
  checks plus 57 Worker checks each, with forced memory growth, lifetime checks and explicit C++ exception recovery.
- Edge: production-package loader and geometry smoke test; validation exports absent.
- Edge: C++ AddressSanitizer build passes the main-thread and Worker suites without an invalid
  memory access report. This is not a proof of absence of all upstream defects.
- Native SMLib tests and the existing workspace tests remain passing; default and
  feature-enabled Clippy checks pass, including the browser target.

Linux/Windows operating-system browser runs and Emscripten build hosts were not
executed here. The shipped browser files contain Wasm and JavaScript, not a
macOS binary dependency. Native Windows library support is a separate outstanding
item, as described above.
