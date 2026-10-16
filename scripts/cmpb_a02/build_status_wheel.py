#!/usr/bin/env python3
"""Build the local A02 repair artifact and a separate pinned installed environment.

No uploads or version edits. Run under run_stage.py's authorized repair ledger.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile

from run_stage import digest, dump


def inspect_members(names):
    if any(Path(n).suffix.lower() in {".jar", ".class"} or n.endswith("CIRI-full_manual.txt") for n in names):
        raise AssertionError("Excluded external-detector asset in distribution")
    for suffix in (
        "circyto/resources/tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh",
        "circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl",
        "circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/LICENSE.txt",
        "circyto/resources/tools/CIRI-full_v2.0/bin/CIRI_AS_v1.2/LICENSE.txt",
        "circyto/pipeline/import_ciri3.py",
    ):
        assert any(n.endswith(suffix) for n in names), suffix


def build(root, repo, commit):
    commands = []

    def run(argv, cwd=None):
        commands.append({"argv": [str(x) for x in argv], "cwd": str(cwd or root)})
        dump(root / "records/build_commands.json", commands)
        subprocess.run(argv, cwd=cwd or root, check=True)

    source = root / "build-source"
    source.mkdir()
    archive_path = root / "inputs/repair-source.tar"
    # Select the same build inputs as the existing distribution regression;
    # keep unrelated Git history, manuscripts and external tool trees out.
    inputs = ["circyto", "docs", "examples", "pyproject.toml", "setup.py", "MANIFEST.in",
              "README.md", "LICENSE", "THIRD_PARTY_NOTICES.md", "CITATION.cff"]
    run(["git", "-C", str(repo), "archive", "--format=tar", "--output=" + str(archive_path), commit, *inputs])
    with tarfile.open(archive_path) as archive:
        archive.extractall(source, filter="data")
    run([sys.executable, "-m", "build", "--sdist", "--wheel", "--outdir", str(root / "artifacts"), str(source)])
    wheel = next((root / "artifacts").glob("*.whl"))
    sdist = next((root / "artifacts").glob("*.tar.gz"))
    run([sys.executable, "-m", "twine", "check", str(wheel), str(sdist)])
    with tarfile.open(sdist) as archive:
        inspect_members(archive.getnames())
        archive.extractall(root / "sdist-source", filter="data")
    expanded = next((root / "sdist-source").iterdir())
    run([sys.executable, "-m", "build", "--wheel", "--outdir", str(root / "rebuilt"), str(expanded)])
    rebuilt = next((root / "rebuilt").glob("*.whl"))
    run([sys.executable, "-m", "twine", "check", str(rebuilt)])
    members = {}
    for artifact in [wheel, rebuilt]:
        with zipfile.ZipFile(artifact) as archive:
            inspect_members(archive.namelist())
            assert archive.testzip() is None
            python_members = {}
            for name in archive.namelist():
                if name.startswith("circyto/") and name.endswith(".py"):
                    data = archive.read(name)
                    assert data == (source / name).read_bytes(), name
                    python_members[name] = hashlib.sha256(data).hexdigest()
            members[str(artifact.relative_to(root))] = python_members
    assert list(members.values())[0] == list(members.values())[1]
    env = root / "installed-env"
    run([sys.executable, "-m", "venv", str(env)])
    py = str(env / "bin/python")
    run([py, "-m", "pip", "install", "-r", str(root / "inputs/constraints.txt")])
    run([py, "-m", "pip", "install", "--no-deps", str(wheel)])
    run([py, "-m", "pip", "check"])
    shutil.copytree(repo / "examples/native_ciri3", root / "native_ciri3")
    for path in [repo / "examples/validate_native_ciri3.py", repo / ".github/scripts/check_installed_distribution.py"]:
        shutil.copyfile(path, root / path.name)
    run([py, str(root / "check_installed_distribution.py")])
    for args in [["--help"], ["import-circ", "--help"], ["report", "--help"]]:
        run([str(env / "bin/circyto"), *args])
    run([py, str(root / "validate_native_ciri3.py"), "--fixtures", str(root / "native_ciri3"), "--out", str(root / "installed-native-acceptance")])
    dump(root / "records/build_acceptance.json", {
        "source_commit": commit, "source_archive_sha256": digest(archive_path),
        "artifacts": {str(p.relative_to(root)): {"sha256": digest(p), "bytes": p.stat().st_size}
                      for p in [wheel, sdist, rebuilt]},
        "wheel_python_members": members, "excluded_assets_absent": True,
        "sdist_derived_python_members_equal": True, "clean_installed_acceptance": "PASS",
        "version_policy": "Local unpublished repair retaining metadata 0.10.0; identify by source commit and SHA-256, not the frozen manuscript artifact.",
        "build_environment": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()},
    })


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--commit", required=True)
    args = p.parse_args()
    build(args.root.resolve(), args.repo.resolve(), args.commit)
