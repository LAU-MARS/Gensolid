"""Measure one solid box with Open CASCADE, out of process.

Gensolid does not link OCCT. This helper is the oracle: same sizes and frame,
JSON on stdout. `--probe` exits 0 when pythonocc can be imported.

Install either `cadquery` (provides the `OCP` bindings) or `pythonocc-core`
(`OCC.Core`). The box is centred on `location`, matching `body_create_solid_block`.
OCCT's `MakeBox` places its corner at the frame origin, so the corner is shifted
by half the size along each frame axis.
"""

from __future__ import annotations

import argparse
import json
import sys


def _load():
    try:
        from OCP.BRepBndLib import BRepBndLib
        from OCP.BRepGProp import BRepGProp
        from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
        from OCP.Bnd import Bnd_Box
        from OCP.GProp import GProp_GProps
        from OCP.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
        from OCP.TopExp import TopExp
        from OCP.TopTools import TopTools_IndexedMapOfShape
        from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
    except ImportError:
        from OCC.Core.BRepBndLib import BRepBndLib
        from OCC.Core.BRepGProp import BRepGProp
        from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
        from OCC.Core.Bnd import Bnd_Box
        from OCC.Core.GProp import GProp_GProps
        from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX
        from OCC.Core.TopExp import TopExp
        from OCC.Core.TopTools import TopTools_IndexedMapOfShape
        from OCC.Core.gp import gp_Ax2, gp_Dir, gp_Pnt
    return {
        "MakeBox": BRepPrimAPI_MakeBox,
        "Ax2": gp_Ax2,
        "Dir": gp_Dir,
        "Pnt": gp_Pnt,
        "GProps": GProp_GProps,
        "BRepGProp": BRepGProp,
        "Bnd_Box": Bnd_Box,
        "BRepBndLib": BRepBndLib,
        "TopExp": TopExp,
        "Map": TopTools_IndexedMapOfShape,
        "FACE": TopAbs_FACE,
        "EDGE": TopAbs_EDGE,
        "VERTEX": TopAbs_VERTEX,
    }


def _call(namespace, names, *args):
    for name in names:
        fn = getattr(namespace, name, None)
        if fn is None:
            continue
        try:
            return fn(*args)
        except TypeError:
            continue
    raise RuntimeError(f"none of {names} accepted the arguments")


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _scale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _norm(a):
    length = (_dot(a, a)) ** 0.5
    if length == 0.0:
        raise ValueError("zero-length direction")
    return _scale(a, 1.0 / length)


def _frame(axis, ref):
    z = _norm(axis)
    x = _sub(ref, _scale(z, _dot(ref, z)))
    x = _norm(x)
    y = _cross(z, x)
    return x, y, z


def _count(api, shape, kind):
    indexed = api["Map"]()
    _call(api["TopExp"], ("MapShapes_s", "MapShapes"), shape, kind, indexed)
    extent = getattr(indexed, "Extent", None) or getattr(indexed, "Size")
    return int(extent())


def _xyz(text):
    parts = [float(item) for item in text.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected x,y,z")
    return parts


def measure(args):
    api = _load()
    x_axis, y_axis, z_axis = _frame(args.axis, args.ref)
    corner = _sub(tuple(args.location), _scale(x_axis, args.x * 0.5))
    corner = _sub(corner, _scale(y_axis, args.y * 0.5))
    corner = _sub(corner, _scale(z_axis, args.z * 0.5))
    frame = api["Ax2"](
        api["Pnt"](*corner),
        api["Dir"](*z_axis),
        api["Dir"](*x_axis),
    )
    shape = api["MakeBox"](frame, args.x, args.y, args.z).Shape()

    props = api["GProps"]()
    _call(api["BRepGProp"], ("VolumeProperties_s", "VolumeProperties"), shape, props)
    box = api["Bnd_Box"]()
    _call(api["BRepBndLib"], ("Add_s", "Add"), shape, box)
    bounds = box.Get()
    # Add() enlarges the box by the shape tolerance (OCCT confusion, 1e-7).
    # Report the geometric box so it can be compared with Gensolid's vertex box.
    gap = float(box.GetGap()) if hasattr(box, "GetGap") else 0.0
    report = {
        "faces": _count(api, shape, api["FACE"]),
        "edges": _count(api, shape, api["EDGE"]),
        "vertices": _count(api, shape, api["VERTEX"]),
        "volume": float(props.Mass()),
        "box_min": [float(bounds[0]) + gap, float(bounds[1]) + gap, float(bounds[2]) + gap],
        "box_max": [float(bounds[3]) - gap, float(bounds[4]) - gap, float(bounds[5]) - gap],
    }
    json.dump(report, sys.stdout)
    sys.stdout.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--x", type=float)
    parser.add_argument("--y", type=float)
    parser.add_argument("--z", type=float)
    parser.add_argument("--location", type=_xyz)
    parser.add_argument("--axis", type=_xyz)
    parser.add_argument("--ref", type=_xyz)
    args = parser.parse_args()
    if args.probe:
        _load()
        return 0
    missing = [
        name
        for name in ("x", "y", "z", "location", "axis", "ref")
        if getattr(args, name) is None
    ]
    if missing:
        parser.error("missing " + ", ".join("--" + name for name in missing))
    measure(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ImportError as exc:
        print(f"pythonocc is not importable: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
