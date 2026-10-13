# ABC development sample

The unmodified example linked on the [ABC dataset project page](https://deep-geometry.github.io/abc-dataset/): one STEP document containing three solids, and the associated Parasolid archive with three native text XT documents. This is a smoke sample from one published model, not a representative ABC benchmark or a private L0 set.

`manifest.json` records source URLs, archive members, byte counts and SHA-256 checksums. Reproduce or verify the files with `python tools/fetch-abc-smoke.py`. The fetcher refuses to replace modified files. No model is repaired, normalized, simplified or converted before import.

The CAD models retain their creators' copyright. See the dataset authors' licensing statement and linked Onshape terms on the project page. The repository's Apache-2.0 licence applies to Gensolid source code, not these third-party CAD models.

Citation: Sebastian Koch et al., *ABC: A Big CAD Model Dataset for Geometric Deep Learning*, CVPR 2019.
