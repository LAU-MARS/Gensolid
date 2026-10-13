"""Ask Open CASCADE whether a STEP file is a valid solid.

Gensolid does not link OCCT. `--probe` exits 0 when pythonocc imports.
`--step path` prints {"read": bool, "valid": bool}.
"""

from __future__ import annotations

import argparse
import json
import sys


def load():
    try:
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.IFSelect import IFSelect_RetDone
        from OCP.STEPControl import STEPControl_Reader
    except ImportError:
        from OCC.Core.BRepCheck import BRepCheck_Analyzer
        from OCC.Core.IFSelect import IFSelect_RetDone
        from OCC.Core.STEPControl import STEPControl_Reader
    return STEPControl_Reader, IFSelect_RetDone, BRepCheck_Analyzer


def check(path):
    reader_type, done, analyzer_type = load()
    reader = reader_type()
    status = reader.ReadFile(str(path))
    if status != done:
        return {"read": False, "valid": False}
    # A geometry-only file has no product root, so transfer each solid by rank.
    model = reader.StepModel()
    for index in range(1, model.NbEntities() + 1):
        entity = model.Entity(index)
        if entity.DynamicType().Name() == "StepShape_ManifoldSolidBrep":
            reader.TransferOne(index)
    if reader.NbShapes() < 1:
        return {"read": False, "valid": False}
    valid = True
    for index in range(1, reader.NbShapes() + 1):
        valid = valid and bool(analyzer_type(reader.Shape(index)).IsValid())
    return {"read": True, "valid": valid}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--step")
    args = parser.parse_args()
    if args.probe:
        load()
        return 0
    if not args.step:
        parser.error("missing --step")
    json.dump(check(args.step), sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ImportError as exc:
        print(f"pythonocc is not importable: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
