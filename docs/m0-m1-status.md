# M0/M1 public development evidence and M2 foundation

Recorded on 2026-10-13. Acceptance scope is the public development benchmark, with an aggregate interface reserved for future private scores. The human-owned verifier, checker, tolerance module and policy assertions are unchanged by this work. Existing workspace changes are preserved and protected-file hashes are compared with the start-of-task snapshot.

Final validation: `cargo test` passed 198 tests, zero failed or ignored; `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings` and `git diff --check` passed. The C header was regenerated. All 17 protected-file hashes and the protected-file inventory match the start-of-task snapshot. The complete public benchmark exited successfully and verified all 22 frozen files, including 21 CAD inputs.

Generated evidence: [complete report](../target/m0-m1-m2-report.json), [compact metrics](../target/m0-m1-m2-summary.json) and [test log](../target/m1-m2-tests.log). These generated target files are local build artifacts and are reproduced by the commands below.

## Reproduce

```sh
python tools/fetch-abc-smoke.py
cargo run -p gensolid-cli -- check corpus/abc-smoke
cargo run -p gensolid-cli -- benchmark --output target/m0-m1-m2-report.json
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test
```

`check` accepts files and recursively scanned directories. Exit 0 means every selected file is checked valid, 1 means a supported body is invalid, and 2 means missing, malformed, unsupported or unchecked input, or an empty workload. Every selected input gets a JSON diagnostic. `benchmark` returns 0 when report generation succeeds and frozen input integrity is intact; manifest/hash failures return 2 while retaining the report. Successful report generation is separate from geometry coverage and recorded benchmark failures.

Use `--corpus <directory>` for a broader local STEP/XT subset. Keep the frozen input list when comparing runs. The included ABC sample is the single model linked by its authors: one STEP document and three native XT parts. URLs, archive names and SHA-256 checksums are in [the manifest](../corpus/abc-smoke/manifest.json); the fetcher verifies hashes and refuses to overwrite changed data.

The default corpus contains 21 CAD documents: four ABC documents, twelve original randomly sized/rotated OCCT solids, two deliberately damaged STEP files and three OCCT void solids. Each directory has a frozen manifest. Expected CAD entries are unioned with discovered files; missing files remain in the denominator. Hash failures are separate from body validity. Empty or malformed manifests cannot establish integrity. The CLI does not scan private cases.

## M0

The real STEP document imports three bodies. All three native text XT files import one body each using exact schema `SCH_3001168_30100_13006`. All six bodies pass `check_body`, with zero faults and zero skipped checks. Four-file diagnostic coverage is 4/4.

OCCT runs out of process. The immutable verifier measures four public blocks, agreeing 4/4. A development helper reuses its binding selector and `BRepCheck_Analyzer`, transfers every source solid including void bodies, and refuses a definitive verdict on partial transfer. Complete corpus results are 15 valid files, one invalid file, four explicitly unsupported files and one malformed file: definitive Gensolid coverage 16/21. STEP comparisons agree 12/13 (92.31%), with comparison coverage 13/18. Positive supported STEP inputs agree 12/12. One perturbed-point file differs: OCCT retains the point but increases an imported vertex tolerance to about one model unit, while Gensolid preserves source coordinates and selected precision. A dangling-reference file makes this installed OCCT process fail with exit code `0xc0000005`. All disagreements, unsupported inputs and oracle failures remain visible.

These are small-sample results. OCCT has no native XT reader here. Paired STEP/XT evidence compares per-body face/edge counts and original XT vertex positions. STEP retains extra vertices on ring edges; native XT can represent them without vertices. Two bodies have measured vertex residuals below session linear precision; one has counts-only evidence. Pairing does not compare curved interiors.

XT decoding uses pinned `parasolid-core` 0.3.7 with an original adapter auditing embedded-field provenance and source/target topology reachability. Unknown schemas, unsupported wire/acorn bodies, convergent carriers and nested-body extensions cannot yield a partial body reported as checked. Neutral-binary evidence currently includes independently encoded synthetic blocks; real binary producer evidence remains to be added.

The STEP reader validates the Part21 envelope and all references across DATA sections, preserves source geometry and topology directions, and distinguishes malformed, invalid and unsupported input. Tests include real OCCT exports, the ABC sample, random frames, rational B-spline faces and near-degenerate dimensions. `BREP_WITH_VOIDS` supports certifiable convex planar shells with straight edges, strict containment and independently separated cavities. The hollow block has volume 168 and passes all checks. Randomized single/double-cavity tests cover size, frames and near-degenerate gaps. Curved/nonconvex or uncertain void relationships, assembly placements, unsupported body classes, spherical vertex loops, unsupported procedural/trimmed geometry and external references return explicit errors. Planar writing refuses finite voids rather than losing them. Numeric source coordinates are preserved; callers choose precision in source units.

## M1

Geometry includes analytic/NURBS evaluation and first derivatives, exact orientation predicates, interval sign classification, analytic projections and explicitly labelled numerical projection candidates. Candidate projection does not certify the global minimum on arbitrary NURBS geometry.

Supported intersection families include line with plane/sphere/cylinder/cone; circle with plane/sphere; ellipse with plane; plane with plane/sphere/cylinder and supported cone/torus sections; sphere with sphere; parallel cylinders; and coaxial sphere/cylinder. Additional complete families include:

- Circle/cylinder with circle plane perpendicular to the axis, arbitrary circle centres, tangencies and overlap.
- Ellipse/sphere with sphere centre on the ellipse axis, including four isolated points.
- Positive-weight rational BCurve/plane of degree one, two or three, with homogeneous span extraction and root isolation. Cubics include nonuniform, repeated and unclamped knots and multiple roots; unresolved root clusters return `Unimplemented`.
- Positive-weight rational bilinear BSurface/plane, retaining disconnected branches as rational quadratic curves where partition signs can be certified. Unresolved common-zero crossings and area overlap remain unsupported.
- Oblique plane/cone ellipses, apex points and apex generator lines on either nappe.
- Meridional plane/torus and coaxial sphere/torus, with complete circle branches.
- Off-axis sphere/torus certified empty sets, extreme-point tangencies and a special complete one/two-meridian-circle family. General quartic curves remain unsupported.

The fixed analytic catalogue passes 24/24 in two frames. Original deterministic workload: seed `20261013`, fourteen known relative configurations and 64 random scales/rigid frames, 896/896 matched, zero unimplemented and zero failed. Maximum sampled point residual is about `9.99e-11` source units. The unchanged extension of cubic curve/plane, bilinear NURBS surface/plane and oblique cone/plane configurations now passes 192/192, zero unimplemented and zero failed, maximum sampled residual about `9.81e-11`. Combined randomized coverage is 1088/1088 (100%). Original inputs, expected branch counts and denominators are unchanged. Residuals are within `Precision::DEFAULT.linear()`.

Scores compare known reference branch counts and residuals on isolated points, curve samples and original-geometry overlap witnesses. They do not prove branch completeness for arbitrary random pairs without an independent answer. Higher-degree NURBS curve/surface intersections now have a separate interval-certified path; general NURBS surface/surface intersections have an explicitly non-certified numerical tracing report (see [advanced operations](advanced-modeling.md)). Exact nonempty general NURBS surface/surface results, skew cylinders, non-apex hyperbolic/parabolic cone sections, general oblique torus sections and several off-axis round-surface pairs remain `Unimplemented`. Mixtures of isolated points and partially overlapping spans require a broader result representation. The requested remaining 192 development cases are complete; M1 is not a solver for all curve/surface classes.

Earlier randomized regressions fixed frame-orthogonalization cancellation and projection line-search rejection on a distance-rounding plateau. The actual failing proptest seed remains recorded.

## M2 initial implementation

`gensolid-modeling` is a checked layer above topology and the immutable checker. It provides blocks, convex polyhedra, convex polygon extrusion, linear vector sweeps and two-section positive homothetic lofts with planar side faces. A line-edge split preserves its carrier and updates incident fins; the checked entry point validates every affected body before and after mutation. The Rust session API exposes these operations. Randomized tests cover sizes, coordinate frames, orientations, oblique sweeps, shifted/scaled lofts, near-degenerate dimensions and edge splitting. Unsupported or degenerate input returns an error.

Unit-density mass properties integrate original planar B-rep loops directly: volume, material boundary area, centre of mass and inertia at the centre in world axes. Oriented holes and finite cavities are included; material/material interfaces do not contribute to boundary area. Signed-volume and oriented-area decisions use outward-rounded intervals. Coordinate subtraction is enclosed before volume certification; compensated sums and direct transverse second-moment sums retain small axial inertia in thin rods. Curved carriers remain unsupported for mass integration.

`gensolid-mesh` is output-only. No B-rep operation or mass calculation depends on it. It triangulates supported straight-edged planar faces and analytic rectangular parameter patches on planes, cylinders, spheres, cones and tori. Analytic patches use derivative bounds for interpolation error and interval trigonometric enclosures for vertex roundoff. Returned results must satisfy the requested chord tolerance. The patch API does not claim watertight meshes of arbitrary trimmed curved bodies. Unsupported trims, planar holes, NURBS curvature bounds and resource exhaustion return explicit errors.

Planar ear classification uses interval enclosures of the exact plane projection, including coordinate subtraction and UV arithmetic. Uncertain orientations, including unresolved collinear ears, return `Unimplemented`. Its deviation budget includes each line carrier's endpoint projection error plus edge tube radius, vertex tolerance balls and the plane-shadow displacement. A regression moves a shared edge toward the exterior corner of two planes: the body checker accepts the tolerant representation, but a chord budget below the combined displacement is correctly rejected.

M2 is started, not complete. Convex off-axis revolution, corresponding convex multi-section solid lofts and compatible NURBS surface lofts are now available in the [advanced operations](advanced-modeling.md). General Euler operators, other analytic curved solid primitives, curved-path sweeps, unrestricted solid lofts, curved-body mass integration and complete trimmed-body tessellation remain future work. The public report records native OCCT mass comparisons and mesh deviation evidence for the supported foundation.

The M2 run uses seed `20261013`, 16 frames/scales and three operations per frame: 48/48 checked bodies, 48/48 independent input-volume answers and 48/48 native OCCT mass comparisons within precision-derived budgets. Maximum relative differences:

| Quantity | Maximum difference |
|---|---:|
| Volume | 2.9105508855529445e-12 |
| Surface area | 1.2077874848834528e-13 |
| Centroid distance / body diameter | 1.964771689522121e-12 |
| Inertia maximum-component norm | 8.241001499376445e-7 |

The inertia statistic divides the largest absolute tensor-component difference by the largest absolute OCCT tensor component. Its maximum occurs on a small translated loft; the original difference is retained rather than normalized away by moving the benchmark inputs. All 48 planar meshes pass a chord budget of `6.4e-7`, with maximum returned bound `5.0000948130463056e-9`. All three analytic patch meshes pass a chord budget of `0.00065536`; the maximum returned bound is `0.00026630901864713245`. These are model-unit deviations for the explicitly supported domains above.

## Reserved private score contract

An independent scorer can supply an aggregate receipt outside the repository:

```json
{"suite_id":"held-out-v1","total":200,"passed":187}
```

```sh
cargo run -p gensolid-cli -- benchmark --private-score /path/to/receipt.json --output target/m0-m1-m2-report.json
```

Receipts require a nonempty suite ID, positive total and `passed <= total`. Unknown fields are rejected, excluding raw cases. The report records receipt supply; it does not independently certify the service or declare L0 acceptance. No private test set is created or required for this run.

## Public authorities

XT compatibility uses Siemens' [JT v10 Rev-B, Appendix E](https://www.plm.automation.siemens.com/en_us/Images/JT-v10-file-format-reference-rev-B_tcm1023-233786.pdf), describing XT B-rep and embedded schemas. The historic 2008 XT Format Reference Siemens URL was unavailable. No Parasolid manual prose or OCCT implementation code is included.

STEP grammar/entity semantics were checked against [Part21 third edition](https://www.steptools.com/stds/step/IS_final_p21e3.html), [EDGE_CURVE](https://www.steptools.com/stds/stp_aim/html/t_edge_curve.html) and [SURFACE_CURVE](https://www.steptools.com/stds/stp_aim/html/t_surface_curve.html). Corpus provenance: [ABC authors' project page](https://deep-geometry.github.io/abc-dataset/). Mass integration uses original boundary formulas; mathematical context: [Mirtich, JGT 1(2), 1996](https://people.eecs.berkeley.edu/~jfc/mirtich/massProps.html), with no reference code copied.

Native mass comparisons use the documented [BRepGProp](https://dev.opencascade.org/doc/refman/html/class_b_rep_g_prop.html) and [GProp_GProps](https://dev.opencascade.org/doc/refman/html/class_g_prop___g_props.html) APIs. `UseTriangulation` is false; solids are checked closed and valid before mass evaluation. OCCT independently reconstructs the same planar boundary vertices and loops. This checks mass integration on those carriers; closed-form profile-volume checks separately test the intended extrusion/sweep/loft volume.
