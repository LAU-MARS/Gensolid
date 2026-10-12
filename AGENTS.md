# Agent Development Rules

This file is written for the coding agents (and humans) working on Gensolid. The rules are maintained by humans; agents must not modify this file.

## What this project is

A B-rep kernel benchmarked against Parasolid. The data model, naming, and behaviour take Parasolid as the reference; correctness is decided by `gensolid-check` and `verifier/`, not by "it seems to run".

## Hard rules (CI enforces most of them)

1. **Don't touch the ruler.** Do not modify `verifier/`, `crates/gensolid-check/`, `crates/gensolid-core/src/tolerance.rs`, `AGENTS.md`, or `.github/`, and never change test assertions to make tests pass. If you believe a test is wrong, say so in the PR description and hand it to a human.
2. **No literal tolerances.** Every threshold is derived from `Precision` or the entity's `Tolerance`. Only `tolerance.rs` may contain constants such as `1e-8` (`verifier/policy` checks this).
3. **No tessellation shortcuts.** Exact B-rep operations — Booleans, fillets, intersections — must never tessellate first and compute on the mesh. B-rep crates must not depend on the mesh crate (`verifier/policy` checks this).
4. **No bare floating-point comparisons in topological decisions.** Use the exact predicates from `gensolid_math::predicates`, or the deterministic signs of `Interval`; fall back to tolerance only when neither can decide, and write that fallback explicitly in the code.
5. **Library code never panics.** `unwrap`, `expect`, `panic!`, and `todo!` are forbidden; return errors carrying an `ErrorCode`. When an algorithm is not implemented yet, return `ErrorCode::Unimplemented` — never return a plausible-looking wrong result.
6. **`unsafe` only in `gensolid-ffi`.** Every `unsafe` block carries a `// SAFETY:` justification.
7. **The output of every new modelling operation must pass `check_body`**, backed by proptest property tests covering random sizes, random coordinate frames, and inputs near degenerate boundaries.
8. **Code provenance:**
   - Do not translate or adapt OCCT source code (LGPL). Reading its documentation, papers, and the inputs and outputs of its test cases is fine.
   - Porting code from truck (Apache-2.0) is allowed: record the source repository and commit in the file header, and add the copyright notice to `NOTICE`.
   - Do not paste text from the Parasolid manuals. The XT format specification is used only to implement format compatibility.

## Conventions

- Face normal = surface normal × `Face::sense`; the `front` shell lies on the side the normal points to.
- The fins of a loop run counter-clockwise around the face normal, with material on the left.
- `Body::regions[0]` is always the infinite empty region.
- Surface parameterizations make the natural outward normal `du × dv` (away from the axis; for planes, the frame's `z`).
- API naming: `<class>_<verb>_<object>`, e.g. `body_create_solid_block`, `body_ask_faces`; the C interface adds the `GS_` prefix and upper-cases the class name.

## Before committing

```sh
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test
```

After changing the C interface, regenerate the header:

```sh
cbindgen --config crates/gensolid-ffi/cbindgen.toml --crate gensolid-ffi \
  --output crates/gensolid-ffi/include/gensolid.h crates/gensolid-ffi
```
