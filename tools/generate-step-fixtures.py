"""Generate independent STEP test inputs through OCCT's public Python API.

Only black-box inputs/outputs are used; no OCCT source code is incorporated into
Gensolid. Generated files describe original primitive combinations, not ABC data.
Requires the optional OCP package; it is never linked into the Rust kernel.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepPrimAPI import (
        BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeCone,
        BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus,
    )
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCP.gp import gp_Pnt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    block = BRepPrimAPI_MakeBox(2.0, 3.0, 4.0).Shape()
    shapes = {
        "block": block,
        "cylinder": BRepPrimAPI_MakeCylinder(2.0, 3.0).Shape(),
        "cone": BRepPrimAPI_MakeCone(2.0, 1.0, 3.0).Shape(),
        "sphere": BRepPrimAPI_MakeSphere(2.0).Shape(),
        "torus": BRepPrimAPI_MakeTorus(3.0, 1.0).Shape(),
        "block_with_hole": BRepAlgoAPI_Cut(
            block, BRepPrimAPI_MakeBox(gp_Pnt(0.5, 0.5, -1.0), 1.0, 1.0, 6.0).Shape()
        ).Shape(),
    }
    for name, shape in shapes.items():
        if not BRepCheck_Analyzer(shape).IsValid():
            raise RuntimeError(f"OCCT generated invalid {name}")
        writer = STEPControl_Writer()
        if writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone:
            raise RuntimeError(f"OCCT transfer failed: {name}")
        path = args.output / f"{name}.step"
        if writer.Write(str(path)) != IFSelect_RetDone:
            raise RuntimeError(f"OCCT write failed: {path}")


if __name__ == "__main__":
    main()
