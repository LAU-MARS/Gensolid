# Planar void STEP export: trusted capability contract v1

This contract was prepared separately from candidate implementation on 2026-10-22, following the project owner's explicit request to resolve frozen “unimplemented” tests before a new RSI run.

## Maintenance boundary

The historical final assertion in `independent_occt_void_export_preserves_regions_and_volume` stated that `write_planar_solid` must return `Unimplemented`. That described the old capability state, rather than an invariant that should survive implementation.

Its replacement has two strict branches: an unimplemented result must carry exactly `ErrorCode::Unimplemented`; a successful result must preserve the complete solid and its finite void. All preceding assertions and all existing invalid, malformed, touching, reversed, outside and overlapping input tests remain in force. No checker, tolerance policy, corpus, CI policy or existing public benchmark was changed.

The two `#[cfg(test)] mod ...` declarations in STEP `lib.rs` were moved to the end without changing the compiled tests. RSI's conservative inline-test hash previously treated every production function following those early declarations as frozen.

This is a **pre-run maintenance migration**, not a worker exemption. Once this seed commit is selected, candidates may not edit the migrated tests, the probe, its helper, fixtures, or metric definitions. Future capability-state migrations require another separately reviewed seed contract before a run; candidates cannot rewrite expectations to admit themselves.

## Fixed workload

The executable probe is `cargo run -p gensolid-step --example rsi_void_export_contract`. Its independently defined analytic boxes, shell structure, coordinates, orientation, volume and checker results are checked by `examples/support/void_export_assertions.rs`.

The contract contains 16 fixed cases:

- Eight positive export cases: independent OCCT hollow-block fixture; simple/complex one-cavity and two-cavity inputs; transformed frames; small/large scales; and a near-precision case.
- One ordinary single-shell export control.
- Seven rejection cases: an off-carrier vertex, broken fin owner, empty shell, duplicate shell, broken region backpointer, a geometrically valid closed cavity outside the outer shell, and geometrically valid closed cavities overlapping each other.

The last two cases retain valid individual carriers/manifolds, demonstrated before calling the writer, so calling `check_body` alone is insufficient. Containment and separation must also be established.

Every positive case first validates its input independently. `Unimplemented` records an unsupported outcome without failing the seed's ordinary regression checks. A successful export must reimport as exactly one valid body with no unchecked checks, preserve all expected regions and shell membership, and match independent coordinates and analytic material volume. Returning an outer box, dropping a cavity, memorizing one fixture, or regrouping shells does not satisfy the contract.

## Baseline and advancement

Before implementation, the measured baseline is **8 passed / 16**, with **all eight positive void-export cases unsupported**. The passed cases are controls and explicit rejections; this is not “50% void-export capability.”

The RSI v3 adapter binds the probe, assertion helper and independent fixture hashes, validates all case IDs and counts, and records `step_planar_void_export_certified`. A changed contract is incomparable; failed checks reject the candidate; loss of any previously passing case is a regression even if the aggregate count rises.

This run's acceptance target is **16/16**, along with the unchanged public benchmark, formatting, Clippy, the full test suite, exact-commit independent review, and a fully re-evaluated integration merge. Passing the contract proves only this bounded planar export behavior, not general curved void handling, arbitrary B-rep correctness or Parasolid equivalence.

The probe uses Gensolid's reader backed by independent analytic invariants. External STEP interoperability is a separate check, rather than a claim derived from self round-tripping alone.
