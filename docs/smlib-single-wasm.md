# Single-Wasm SMLib package

The selected route is a unified `wasm32-unknown-emscripten` build. Rust, SMLib and
the C++ runtime are statically linked into **one `gensolid.wasm`**, with an ES-module
loader. The JS loader remains necessary; one Wasm module does not mean one file
containing all application JavaScript.

## Why this route

| Route | Benefit | Maintenance cost |
|---|---|---|
| Rust + C++ on Emscripten (selected) | Direct Rust/C++ calls and one memory; uses the existing C++ runtime | Pin the SDK and Rust toolchain, rebuild std, maintain a small browser ABI |
| Keep wasm-bindgen/unknown-unknown and integrate a C++ runtime | Keeps the existing binding generator | Own the runtime, allocation, exception and import compatibility integration |

The current browser-facing API is small, while SMLib relies heavily on C++ runtime
features. Moving the binding layer is the smaller change. The second route is not
impossible, but it has not been implemented or experimentally compared here.
[The Rust target documentation](https://doc.rust-lang.org/stable/rustc/platform-support/wasm32-unknown-emscripten.html)
describes the C/C++ interoperability and recommends rebuilding std with matching
Emscripten ABI settings.

The geometry algorithms and the safe Rust curve/surface facade are shared with
the native build. This target selects direct C ABI calls instead of the dual-module
JS transport. Results still use the existing ownership/copy rules; this is not a
claim that every API becomes zero-copy. The single memory removes the *additional*
copy between separate Rust and C++ Wasm memories.

## Build

Activate Emscripten **6.0.12**, then install the fixed Rust compiler and source:

```sh
rustup toolchain install nightly-2026-10-01 --profile minimal --component rust-src \
  --target wasm32-unknown-emscripten
export SMLIB_SOURCE_DIR="$PWD/target/usd-brep" # same pinned checkout as docs/smlib.md
python3 tools/package-single-wasm.py --output target/smlib-single-web
```

The script rebuilds Rust std using `-Zbuild-std=std,panic_abort`, builds the SMLib
archive and links everything with Emscripten. Rust panics abort; C++ exceptions
are caught inside the existing C ABI using Emscripten's JS exception runtime.
The native public C header and the Rust-only wasm-bindgen build remain separate
supported entry points. The older dual-Wasm package command remains available
for comparison and existing consumers.

The final directory contains exactly one `.wasm`, Emscripten JS, `index.mjs`,
TypeScript declarations, notices and a SHA-256 manifest. No SMLib companion module,
wasm-bindgen snippets, or browser-to-SMLib JS transport is required. Deploy the
whole directory together.

```js
import load, {Session, SmlibBackend} from './smlib-single-web/index.mjs';
await load();
const backend = new SmlibBackend();
const curve = backend.circle([0, 0, 0], 2);
backend.free();
console.log(curve.evaluate(curve.domain()[0]));
curve.free();
const session = new Session();
const body = session.bodyCreateSolidBlock(2, 3, 4);
console.log(session.bodyCheck(body), session.bodyAskVolume(body));
session.free();
```

The existing `Session`, `SmlibBackend`, `SmlibCurve` and `SmlibSurface` method names
and typed-array return shapes are preserved. Initialization is now just `await
load()`; there is no separate `initializeSmlib` call. Freed/stale handles are
rejected. Geometry retains its backend through Rust reference counting. JS error
objects returned by the Rust dispatcher include a numeric `code`.
All Rust curve/surface methods remain available for application-specific Rust
exports. The convenience JS classes have the same scope as the dual-module API;
this change does not add previously unexposed modelling operations.

## Validation

```sh
python3 tools/package-single-wasm.py --validation --output target/smlib-single-validation
npm install --prefix tools/wasm
node tools/wasm/test-browser.mjs target/smlib-single-validation
BROWSER=firefox node tools/wasm/test-browser.mjs target/smlib-single-validation
BROWSER=webkit node tools/wasm/test-browser.mjs target/smlib-single-validation
node tools/wasm/test-browser.mjs target/smlib-single-web
```

The same shared Rust geometry assertions run through direct C++ calls, including
rational NURBS, derivatives, intersections, construction, editing and projection.
The browser harness checks forced memory growth, retained geometry after backend
release, C++ exception recovery, the Rust body checker, Worker isolation, and the
absence of imported SMLib JS bridge functions. Production packages omit diagnostic
exports. The eight-MiB stack and stack-overflow checks are preserved from the
validated SMLib browser configuration.

C++ AddressSanitizer is available with `--validation --sanitize` and a separate
output directory. It instruments the C++ side and allocator, not Rust instructions.
This diagnostic package is not the production artifact.

## Results on this branch

Validated on a macOS Apple Silicon build/test host:

- Edge / Chromium 155.0.4283.45, Firefox 153.0 and WebKit 26.5 each passed the
  single-module suite: 684 main-thread cases plus 57 cases in a module Worker.
- The harness also passed the existing Session query/placed-block API, numeric
  error codes, released-handle rejection, forced memory growth and C++ exceptions.
- The production loader passed its separate smoke test, with diagnostic exports
  absent and exactly one Wasm binary in the release directory (about 5.2 MiB).
- The C++ AddressSanitizer variant passed the main-thread and Worker suite in Edge.
- Default workspace tests, native SMLib tests, workspace Clippy and the Emscripten
  target Clippy check passed. The public C header was regenerated and is unchanged.
- The earlier dual-module browser regression still passes after sharing the Rust
  validation runner; its assertions have not been relaxed.

Linux/Windows browser OS runs and build hosts were not executed here. The released
Wasm files are independent of the macOS build host. No speedup ratio is claimed:
removing the extra cross-module transport is structural; workload benchmarking is
separate from correctness validation.
