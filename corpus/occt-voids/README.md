# Original OCCT void solids

Three unmodified black-box OCCT exports: a hollow block, a hollow solid with a nonconvex outer boundary, and a hollow cylinder. Source fixtures and hashes are recorded in `manifest.json`.

The hollow block's outer dimensions are 4 × 6 × 8 and its cavity dimensions are 2 × 3 × 4; expected material volume is 168. Its STEP import proves containment and passes the immutable body checker. The other two files are valid in OCCT and currently exceed the certified Gensolid void-import scope, so their reports retain explicit unsupported outcomes.
