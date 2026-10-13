"""Fetch the unmodified STEP/XT example linked by the ABC dataset authors.

This development smoke sample is outside the human-owned verifier. Existing files
with different content are refused. The manifest pins every downloaded/extracted byte.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from io import BytesIO
from pathlib import Path
from urllib.request import urlopen
from zipfile import ZipFile


BASE = "https://deep-geometry.github.io/abc-dataset/data/"
SOURCES = {
    "00000050_80d90bfdd2e74e709956122a_step_000.step":
        "579e4b3e18709f6888a4179f189926e40c35a35c156392ec593105cf1650a60c",
    "00000050_80d90bfdd2e74e709956122a_parasolid_000.zip":
        "57ce8354ada505bfad2c2b4b7439761490f19feb85ed88b0041ae365812a888e",
}
MEMBERS = {
    "Part Studio 1 - Part 1.x_t": "ea041f7bd5c3a5afb055b2163e981a7ea1e00cdf994b65ebe7e4e58ea349f888",
    "Part Studio 1 - Part 2.x_t": "6ccbf08a97e5769d2b458b2156872a11b4065021ac364df385344ef738202347",
    "Part Studio 1 - Part 3.x_t": "b5169644405fd3209249dd08f5ff0dc7b44fd5d3bec8a575395f33ed876bac91",
}


def store(path: Path, data: bytes, expected: str) -> None:
    actual = hashlib.sha256(data).hexdigest()
    if actual != expected:
        raise ValueError(f"checksum mismatch for {path.name}: {actual}")
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError(f"refusing to overwrite modified fixture: {path}")
    else:
        path.write_bytes(data)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "corpus/abc-smoke")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    files = []
    for name, checksum in SOURCES.items():
        path = args.output / name
        if path.exists():
            data = path.read_bytes()
        else:
            with urlopen(BASE + name, timeout=45) as response:
                data = response.read()
        store(path, data, checksum)
        files.append({"file": name, "url": BASE + name,
                      "sha256": checksum, "bytes": len(data)})
        if name.endswith(".zip"):
            with ZipFile(BytesIO(data)) as archive:
                for member, member_checksum in sorted(MEMBERS.items()):
                    extracted = archive.read(member)
                    store(args.output / member, extracted, member_checksum)
                    files.append({"file": member, "archive": name, "member": member,
                                  "sha256": member_checksum, "bytes": len(extracted)})
    manifest = {
        "dataset": "ABC v00 official example",
        "source": "https://deep-geometry.github.io/abc-dataset/",
        "selection": "The single example linked by the dataset authors; development smoke sample, not a private acceptance corpus.",
        "files": files,
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {len(files)} files in {args.output}")


if __name__ == "__main__":
    main()
