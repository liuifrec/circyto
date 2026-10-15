"""Inspect wheel/sdist/.conda contents before installation or publication.

For .conda use a conda-build environment (conda_package_streaming is required).
Print archive SHA-256 and retain a complete member inventory in JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import zipfile


REQUIRED = [
    "circyto/pipeline/import_ciri3.py",
    "circyto/resources/tools/CIRI3/README.md",
    "circyto/resources/tools/CIRI3/LICENSE.txt",
    "circyto/resources/tools/CIRI-full_v2.0/bin/ciri2_adapter.sh",
    "circyto/resources/tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh",
    *[f"circyto/resources/tools/CIRI-full_v2.0/bin/{name}/{file}" for name, files in [
        ("CIRI_v2.0.6", ["CIRI2.pl", "CIRI2_manual.txt", "LICENSE.txt"]),
        ("CIRI_AS_v1.2", ["CIRI_AS_v1.2.pl", "CIRI_AS_manual.txt", "LICENSE.txt"]),
    ] for file in files],
]


def members(path):
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            assert archive.testzip() is None
            return archive.namelist()
    if path.suffix == ".conda":
        from conda_package_streaming.package_streaming import stream_conda_component

        names = []
        for component in ("pkg", "info"):
            for _, member in stream_conda_component(str(path), component=component):
                names.append(member.name)
        return names
    with tarfile.open(path) as archive:
        return archive.getnames()


def check(path):
    names = members(path)
    forbidden = [n for n in names if Path(n).suffix.lower() in {".jar", ".class"}
                 or Path(n).name == "CIRI-full_manual.txt"]
    assert not forbidden, f"Excluded assets: {forbidden}"
    for suffix in REQUIRED:
        assert any(n.endswith(suffix) for n in names), f"Missing {suffix}"
    assert any(n.endswith("THIRD_PARTY_NOTICES.md") for n in names)
    if path.name.endswith(".tar.gz"):
        for suffix in ("examples/installed_release_check.py", "examples/native_ciri3/result.txt",
                       "examples/native_ciri3/LICENSE.upstream.txt"):
            assert any(n.endswith(suffix) for n in names), suffix
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size, "excluded_assets": [], "members": sorted(names)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    results = [check(path) for path in args.archives]
    args.out.write_text(json.dumps(results, indent=2) + "\n")
    for record in results:
        print(record["sha256"], record["path"], len(record["members"]), "members: PASS")


if __name__ == "__main__":
    main()
