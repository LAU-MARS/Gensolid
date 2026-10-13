"""Native OCCT polyhedral B-rep construction and geometry-based mass diagnostics.

Input is one UTF-8 JSON batch on stdin. World vertices and outward polygon loops
are reconstructed with MakePolygon, MakeFace, Sewing and MakeSolid. This compares
the same planar carriers; it is not an independent extrusion/loft implementation.
Volume and surface properties explicitly disable triangulation. No OCCT source
code is incorporated, and the immutable verifier is only read for its binding probe.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path


def load():
    sys.dont_write_bytecode = True
    helper = Path(__file__).resolve().parents[1] / "verifier/oracle/occt_check.py"
    spec = importlib.util.spec_from_file_location("gensolid_immutable_occt_probe", helper)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load the immutable OCCT binding probe")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.load()
    try:
        from OCP.BRepBuilderAPI import (
            BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace,
            BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
        )
        from OCP.BRepCheck import BRepCheck_Analyzer
        from OCP.BRepGProp import BRepGProp
        from OCP.BRepLib import BRepLib
        from OCP.GProp import GProp_GProps
        from OCP.TopoDS import TopoDS
        from OCP.gp import gp_Pnt
        return ("OCP", BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace,
                BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
                BRepCheck_Analyzer, GProp_GProps, gp_Pnt, TopoDS.Shell_s,
                BRepGProp.VolumeProperties_s, BRepGProp.SurfaceProperties_s,
                BRepLib.OrientClosedSolid_s)
    except ImportError:
        from OCC.Core.BRepBuilderAPI import (
            BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace,
            BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
        )
        from OCC.Core.BRepCheck import BRepCheck_Analyzer
        from OCC.Core.BRepGProp import brepgprop
        from OCC.Core.BRepLib import breplib
        from OCC.Core.GProp import GProp_GProps
        from OCC.Core.TopoDS import topods
        from OCC.Core.gp import gp_Pnt
        return ("pythonocc", BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace,
                BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid,
                BRepCheck_Analyzer, GProp_GProps, gp_Pnt, topods.Shell,
                brepgprop.VolumeProperties, brepgprop.SurfaceProperties,
                breplib.OrientClosedSolid)


def check_case(case, tolerance, api):
    (_, polygon_type, face_type, sewing_type, solid_type, analyzer_type,
     props_type, point_type, shell_cast, volume_properties, surface_properties,
     orient_closed_solid) = api
    points = case["points"]
    faces = case["faces"]
    if len(points) < 4 or len(faces) < 4:
        raise ValueError("a closed polyhedron needs at least four points and faces")
    if any(len(point) != 3 or not all(math.isfinite(x) for x in point) for point in points):
        raise ValueError("non-finite or malformed world coordinates")
    sewing = sewing_type(tolerance)
    for corners in faces:
        if len(corners) < 3 or len(set(corners)) != len(corners):
            raise ValueError("malformed polygon loop")
        polygon = polygon_type()
        for index in corners:
            if not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < len(points):
                raise ValueError("polygon index outside the vertex array")
            polygon.Add(point_type(*points[index]))
        polygon.Close()
        if not polygon.IsDone():
            raise RuntimeError("OCCT MakePolygon failed")
        face = face_type(polygon.Wire(), True)
        if not face.IsDone():
            raise RuntimeError("OCCT MakeFace could not construct a planar carrier")
        sewing.Add(face.Face())
    sewing.Perform()
    shell = shell_cast(sewing.SewedShape())
    if not shell.Closed():
        return {"id": case["id"], "status": "open_shell", "valid": False,
                "closed": False, "properties": None}
    builder = solid_type(shell)
    if not builder.IsDone():
        raise RuntimeError("OCCT MakeSolid failed")
    solid = builder.Solid()
    if not orient_closed_solid(solid):
        return {"id": case["id"], "status": "unorientable_solid", "valid": False,
                "closed": True, "properties": None}
    valid = bool(analyzer_type(solid).IsValid())
    if not valid:
        return {"id": case["id"], "status": "invalid_solid", "valid": False,
                "closed": True, "properties": None}
    # A nearby accumulation reference avoids cancellation on small translated models.
    reference = point_type(*points[0])
    volume = props_type(reference)
    area = props_type(reference)
    volume_properties(solid, volume, True, False, False)
    surface_properties(solid, area, False, False)
    center = volume.CentreOfMass()
    inertia = volume.MatrixOfInertia()
    properties = {"volume": volume.Mass(), "surface_area": area.Mass(),
                  "centroid": [center.X(), center.Y(), center.Z()],
                  "inertia": [[inertia.Value(i, j) for j in range(1, 4)]
                              for i in range(1, 4)]}
    values = [properties["volume"], properties["surface_area"], *properties["centroid"],
              *(x for row in properties["inertia"] for x in row)]
    if not all(math.isfinite(x) for x in values) or properties["volume"] <= 0:
        raise RuntimeError("OCCT returned non-finite or non-positive mass properties")
    return {"id": case["id"], "status": "ok", "valid": True, "closed": True,
            "properties": properties}


def main():
    sys.stdin.reconfigure(encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    request = json.load(sys.stdin)
    tolerance = float(request["precision"]["linear"])
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("sewing tolerance must be a finite positive session precision")
    try:
        api = load()
    except ImportError as error:
        print(json.dumps({"available": False, "status": "bindings_unavailable",
                          "message": str(error), "cases": []}))
        return
    rows = []
    for case in request["cases"]:
        try:
            rows.append(check_case(case, tolerance, api))
        except Exception as error:
            rows.append({"id": case.get("id"), "status": "construction_error",
                         "valid": False, "closed": False, "properties": None,
                         "message": str(error)})
    print(json.dumps({"available": True, "status": "ok", "backend": api[0],
                      "use_triangulation": False, "cases": rows}, allow_nan=False))


if __name__ == "__main__":
    main()
