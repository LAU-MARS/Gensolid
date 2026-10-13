# ABC subset

Put redistributable STEP files from the ABC dataset here (`.step` or `.stp`). `gensolid-oracle` reads every file in this directory, checks it with Gensolid, and when pythonocc is installed compares that verdict with OCCT's `BRepCheck_Analyzer`.

The agreement rate is `agreed / compared`. Files either kernel cannot read are counted separately and are not part of the rate. This directory is empty until those files are added; the parametric blocks in `../blocks.json` are the cases that run without an external download.
