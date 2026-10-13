# Original OCCT development inputs

Twelve original STEP solids generated through OCCT's public Python API, plus two intentionally mutated documents. These files are independent black-box inputs, separate from the ABC sample and the human-owned verifier.

The positive solids cover randomly sized and rotated blocks, cylinders, cones, spheres, ring tori and blocks with through holes. The sphere cases currently expose unsupported spherical vertex loops. Negative inputs move or delete one referenced point without silently repairing the file.

`manifest.json` records the seed, producer version, construction parameters, file sizes and SHA-256 digests. Reproduce geometry into a fresh directory with `python tools/build-development-corpus.py --output <empty-directory>`; STEP producer metadata may vary between runs, so freeze and reuse the resulting manifest when comparing scores. The generator refuses to overwrite an existing corpus.

The perturbed-point file intentionally reveals different validity policies: Gensolid's original-coordinate import reports a checker fault, while OCCT's STEP import retains the point and increases an imported vertex tolerance to about one model unit, after which `BRepCheck` reports valid. The benchmark records the disagreement. The deleted-point document causes this installed OCCT helper to fail; that failure is recorded separately from valid/invalid verdicts.
