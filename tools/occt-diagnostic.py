"""Public OCCT STEP diagnostics, isolated from the Rust kernel and L0 verifier.

Reuses the immutable verifier's binding loader and BRepCheck_Analyzer. Product
roots and geometry-only MANIFOLD_SOLID_BREP/BREP_WITH_VOIDS entities are transferred
through OCCT's documented API. No OCCT implementation code is incorporated.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True, type=Path)
    args = parser.parse_args()
    helper = Path(__file__).resolve().parents[1] / "verifier/oracle/occt_check.py"
    specification = importlib.util.spec_from_file_location("gensolid_immutable_occt_loader", helper)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load the existing OCCT binding selector")
    module = importlib.util.module_from_spec(specification)
    # Loading a read-only verifier helper must not create bytecode in its directory.
    sys.dont_write_bytecode = True
    specification.loader.exec_module(module)
    reader_type, done, analyzer_type = module.load()
    reader = reader_type()
    if reader.ReadFile(str(args.step)) != done:
        print(json.dumps({"read": False, "valid": False, "expected_solids": 0,
                          "transferred_solids": 0, "partial": False}))
        return
    reader.TransferRoots()
    model = reader.StepModel()
    ranks = [rank for rank in range(1, model.NbEntities() + 1)
             if model.Entity(rank).DynamicType().Name() in {
                 "StepShape_ManifoldSolidBrep", "StepShape_BrepWithVoids"}]
    # Transfer every source solid even if another root already produced a shape.
    # A skipped source solid prevents a definitive whole-file validity verdict.
    transferred = sum(bool(reader.TransferOne(rank)) for rank in ranks)
    partial = transferred != len(ranks)
    count = reader.NbShapes()
    read = count > 0 and not partial and bool(ranks)
    valid = read and all(bool(analyzer_type(reader.Shape(rank)).IsValid())
                              for rank in range(1, count + 1))
    print(json.dumps({"read": read, "valid": valid, "expected_solids": len(ranks),
                      "transferred_solids": transferred, "partial": partial}))


if __name__ == "__main__":
    main()
