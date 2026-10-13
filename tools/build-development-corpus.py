"""Generate independent OCCT STEP development inputs and intentional negative cases.

Uses OCCT's public Python API as a black box; no OCCT implementation is copied or
linked into Gensolid. These original solids are not ABC data or private L0 cases.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
from pathlib import Path


def main() -> None:
    import OCP
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepPrimAPI import (
        BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeCone,
        BRepPrimAPI_MakeSphere, BRepPrimAPI_MakeTorus,
    )
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from OCP.gp import gp_Ax1, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "corpus/occt-development")
    parser.add_argument("--seed", type=int, default=20261013)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("output must be absent or empty; existing frozen inputs are preserved")
    args.output.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    files = []
    first_block = None
    for index in range(2):
        for kind in ["block", "cylinder", "cone", "sphere", "torus", "block_with_hole"]:
            dimensions = [rng.uniform(0.2, 8.0) for _ in range(3)]
            x, y, z = dimensions
            if kind == "block":
                shape = BRepPrimAPI_MakeBox(x, y, z).Shape()
            elif kind == "cylinder":
                shape = BRepPrimAPI_MakeCylinder(x, z).Shape()
            elif kind == "cone":
                shape = BRepPrimAPI_MakeCone(x, x * 0.5, z).Shape()
            elif kind == "sphere":
                shape = BRepPrimAPI_MakeSphere(x).Shape()
            elif kind == "torus":
                shape = BRepPrimAPI_MakeTorus(x * 3.0, x).Shape()
            else:
                shape = BRepAlgoAPI_Cut(
                    BRepPrimAPI_MakeBox(x, y, z).Shape(),
                    BRepPrimAPI_MakeBox(gp_Pnt(x * 0.25, y * 0.25, -z),
                                       x * 0.5, y * 0.5, z * 3.0).Shape()).Shape()
            axis = [rng.uniform(-1.0, 1.0) for _ in range(3)]
            angle = rng.uniform(-3.0, 3.0)
            translation = [rng.uniform(-20.0, 20.0) for _ in range(3)]
            transform = gp_Trsf()
            transform.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(*axis)), angle)
            transform.SetTranslationPart(gp_Vec(*translation))
            shape = BRepBuilderAPI_Transform(shape, transform, True).Shape()
            if not BRepCheck_Analyzer(shape).IsValid():
                raise RuntimeError(f"OCCT generated invalid {kind}/{index}")
            writer = STEPControl_Writer()
            if writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone:
                raise RuntimeError("OCCT STEP transfer failed")
            filename = f"{kind}-{index}.step"
            path = args.output / filename
            if writer.Write(str(path)) != IFSelect_RetDone:
                raise RuntimeError("OCCT STEP write failed")
            data = path.read_bytes()
            files.append({"file": filename, "sha256": hashlib.sha256(data).hexdigest(),
                          "bytes": len(data), "construction": kind, "dimensions": dimensions,
                          "rotation_axis": axis, "rotation_angle": angle,
                          "translation": translation, "producer_valid": True})
            if first_block is None and kind == "block":
                first_block = data.decode("ascii")
    if first_block is None:
        raise RuntimeError("no block generated")
    # Perturb one original vertex carrier without repairing adjacent geometry.
    vertex = re.search(r"#\d+\s*=\s*VERTEX_POINT\s*\([^,]*,\s*#(\d+)\s*\)", first_block)
    if vertex is None:
        raise RuntimeError("STEP producer emitted no expected vertex entity")
    point_id = vertex.group(1)
    point_pattern = re.compile(rf"(#{point_id}\s*=\s*CARTESIAN_POINT\s*\([^,]*,\s*\()([^)]*)(\)\s*\)\s*;)")
    point = point_pattern.search(first_block)
    if point is None:
        raise RuntimeError("STEP producer emitted no expected point entity")
    coordinates = [float(value) for value in point.group(2).split(",")]
    coordinates[0] += 1.0
    replacement = point.group(1) + ",".join(format(value, ".17g") for value in coordinates) + point.group(3)
    negatives = {
        "invalid-vertex.step": point_pattern.sub(lambda _: replacement, first_block, count=1),
        "missing-reference.step": point_pattern.sub("", first_block, count=1),
    }
    for filename, text in negatives.items():
        data = text.encode("ascii")
        (args.output / filename).write_bytes(data)
        files.append({"file": filename, "sha256": hashlib.sha256(data).hexdigest(),
                      "bytes": len(data), "construction": "intentional mutation of block-0",
                      "mutation": "move a referenced point" if filename.startswith("invalid") else "delete a referenced point"})
    (args.output / "manifest.json").write_text(json.dumps({
        "dataset": "original OCCT development solids and negative inputs",
        "producer": f"OCP {OCP.__version__}", "seed": args.seed, "files": files,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Generated and froze {len(files)} STEP inputs in {args.output}")


if __name__ == "__main__":
    main()
