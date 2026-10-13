#!/usr/bin/env python3
"""Reproducible public advanced-operation probes; not the immutable verifier.

Tests real Rust scenes and analytic answers. Curved mass/OCCT comparison and
SSI branch completeness are deliberately NOT claimed by this development suite.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def assess(request, scene):
    mode = request["scene"]
    metrics = scene["metrics"]
    linear = scene["precision"]["linear"]
    if mode in ("revolve", "loft"):
        return scene["status"] == "checked" and metrics["checker_faults"] == 0 and metrics["unchecked"] == 0
    if mode == "smooth_loft":
        return metrics["section_residual"] <= linear
    if mode == "curve_surface":
        return (scene["status"] == "certified" and metrics["roots"] == 3
                and metrics["unresolved_cells"] == 0
                and metrics["root_position_bound"] <= linear
                and all(abs(params[0] - expected) <= linear
                        for params, expected in zip(scene["root_parameters"], (0.2, 0.5, 0.8))))
    # Independent implicit circle answer. This checks all returned sample points,
    # without treating a numerical trace as proof that every branch was found.
    radius = request["radius"]
    points = [point for line in scene["lines"] for point in line["points"]]
    maximum = max((max(abs(math.hypot(p[0], p[1]) - radius), abs(p[2])) for p in points), default=math.inf)
    metrics["analytic_circle_sample_error"] = maximum
    return bool(points) and maximum <= linear and metrics["possible_cells"] > 0 and scene["status"] == "candidate"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "target" / "advanced-development-report.json")
    args = parser.parse_args()
    subprocess.run(["cargo", "build", "-p", "gensolid-viewer"], cwd=ROOT, check=True)
    binary = ROOT / "target" / "debug" / ("gensolid-viewer.exe" if os.name == "nt" else "gensolid-viewer")
    cases = ([{"scene": "revolve", "angle": a} for a in (10, 180, 360)]
             + [{"scene": "loft", "sections": n} for n in (2, 4, 6)]
             + [{"scene": "smooth_loft", "bend": b} for b in (-1, 0, 1)]
             + [{"scene": "curve_surface"}]
             + [{"scene": "surface_surface", "radius": r} for r in (0.3, 0.65)])
    results = []
    for request in cases:
        start = time.monotonic()
        try:
            run = subprocess.run([str(binary)], input=json.dumps(request), text=True,
                                 encoding="utf-8", capture_output=True, cwd=ROOT, timeout=90)
            scene = json.loads(run.stdout)
            passed = run.returncode == 0 and assess(request, scene)
            results.append({"request": request, "passed": passed, "status": scene.get("status"),
                            "metrics": scene.get("metrics"), "message": scene.get("message"),
                            "seconds": time.monotonic() - start})
        except (ValueError, OSError, subprocess.TimeoutExpired) as error:
            results.append({"request": request, "passed": False, "error": str(error)})
        print(json.dumps(results[-1], ensure_ascii=True), flush=True)
    report = {"suite": "advanced-development-v1", "total": len(cases),
              "passed": sum(item["passed"] for item in results), "cases": results,
              "ssi_branch_completeness": "not certified",
              "curved_mass_occt": "not run; curved B-rep mass integration is unsupported",
              "display_chord_error": "not certified; preview only"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    raise SystemExit(0 if report["passed"] == report["total"] else 1)


if __name__ == "__main__":
    main()
