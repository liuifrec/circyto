#!/usr/bin/env python3
"""Record installed A02 origins and external tools; does not run detection."""

import argparse
import importlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import zipfile

from run_stage import digest, dump, stamp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args()
    root = args.root.resolve()
    if (
        Path(sys.prefix).resolve() != root / "env"
        or platform.python_version() != "3.10.22"
    ):
        raise RuntimeError("Wrong Python prefix/version")
    if os.environ.get("PYTHONPATH") or os.environ.get("PYTHONNOUSERSITE") != "1":
        raise RuntimeError("Python isolation not applied")
    expected = dict(
        anndata="0.11.4", mudata="0.3.10", numpy="2.2.6", scipy="1.15.3", pandas="2.3.3"
    )
    versions = {
        dist.metadata["Name"]: dist.version for dist in metadata.distributions()
    }
    for name, version in expected.items():
        if metadata.version(name) != version:
            raise RuntimeError("Dependency constraint mismatch: " + name)
    origins = {
        name: str(Path(importlib.import_module(name).__file__).resolve())
        for name in ["circyto", "anndata", "mudata", "numpy", "scipy", "pandas", "h5py"]
    }
    if not all(Path(value).is_relative_to(root / "env") for value in origins.values()):
        raise RuntimeError("A module was imported outside the isolated environment")
    wheel = root / "inputs/artifacts/circyto-0.10.0-py3-none-any.whl"
    installed_parent = Path(origins["circyto"]).parent.parent
    matches = []
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            if name.startswith("circyto/") and name.endswith(".py"):
                if archive.read(name) != (installed_parent / name).read_bytes():
                    raise RuntimeError("Installed source differs from wheel: " + name)
                matches.append(name)
    helpers = root / "helpers/runtime"
    helpers.mkdir(parents=True, exist_ok=True)
    for name in ["run_stage.py", "tool_audit.py"]:
        target = helpers / name
        if target.exists():
            raise FileExistsError(target)
        shutil.copyfile(Path(__file__).with_name(name), target)
    bindir = root / "audit-bin"
    bindir.mkdir(exist_ok=True)
    binaries = {"bwa": root / "env-bwa/bin/bwa", "java": root / "env/bin/java"}
    for name, actual in binaries.items():
        argv = [
            "/usr/bin/python3",
            str(helpers / "tool_audit.py"),
            "--root",
            str(root),
            "--tool",
            name,
            "--real",
            str(actual),
            "--",
        ]
        wrapper = bindir / name
        wrapper.write_text("#!/bin/sh\nexec " + shlex.join(argv) + ' "$@"\n')
        wrapper.chmod(0o755)
    env = dict(
        PATH=str(bindir)
        + ":"
        + str(root / "env/bin")
        + ":"
        + str(root / "env-bwa/bin")
        + ":/usr/bin:/bin",
        PYTHONNOUSERSITE="1",
        CIRCYTO_CIRI3_JAR=str(root / "external/CIRI3_Java_18.0.1.jar"),
        CIRCYTO_CIRI3_JAVA=str(bindir / "java"),
        CIRCYTO_CIRI3_EXTRA_ARGS="-S 0 -Ma 0",
        JAVA_TOOL_OPTIONS="-Xmx16g -XX:ActiveProcessorCount=8",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        NUMEXPR_NUM_THREADS="1",
    )
    os.environ.update(env)
    dump(root / "records/execution_environment.json", env)
    commands = []

    def capture(name, argv, acceptable=(0,)):
        result = subprocess.run(argv, capture_output=True, text=True, check=False)
        target = root / "records" / (name + ".txt")
        target.write_text(result.stdout + result.stderr)
        commands.append(
            {
                "argv": argv,
                "exit_code": result.returncode,
                "output": str(target.relative_to(root)),
                "sha256": digest(target),
            }
        )
        if result.returncode not in acceptable:
            raise RuntimeError("Preflight command failed: " + name)
        return result.stdout + result.stderr

    capture("pip-check", [sys.executable, "-m", "pip", "check"])
    capture("pip-freeze", [sys.executable, "-m", "pip", "freeze", "--all"])
    bwa = capture("bwa-version", [str(bindir / "bwa")], (1,))
    java = capture("java-version", [str(bindir / "java"), "-version"])
    if "0.7.17-r1188" not in bwa or 'version "17.' not in java:
        raise RuntimeError("External executable version mismatch")
    capture("cli-help", [str(root / "env/bin/circyto"), "--help"])
    capture("doctor", [str(root / "env/bin/circyto"), "doctor"])
    # Runtime configuration inspection is deliberately outside the value oracle.
    from circyto.detectors.ciri3 import inspect_ciri3_runtime

    runtime = inspect_ciri3_runtime()
    if not runtime["direct_ready"] or runtime["jar"] != env["CIRCYTO_CIRI3_JAR"]:
        raise RuntimeError("Configured CIRI3 is not ready at the reviewed JAR path")
    for name, prefix in [("python-java", root / "env"), ("bwa", root / "env-bwa")]:
        capture(
            "conda-explicit-" + name,
            [
                "/home/liuyuchen/.local/bin/micromamba",
                "--no-rc",
                "list",
                "--prefix",
                str(prefix),
                "--explicit",
            ],
        )
    dump(
        root / "records/environment_manifest.json",
        {
            "recorded_utc": stamp(),
            "python": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
            "module_origins": origins,
            "resolved_distributions": versions,
            "wheel_sha256": digest(wheel),
            "installed_python_members_equal_wheel": matches,
            "external_executables": {
                name: {"path": str(path.resolve()), "sha256": digest(path.resolve())}
                for name, path in binaries.items()
            },
            "external_jar_sha256": digest(root / "external/CIRI3_Java_18.0.1.jar"),
            "ciri3_configuration": runtime,
            "environment": env,
            "commands": commands,
            "helper_sha256": {
                name: digest(helpers / name)
                for name in ["run_stage.py", "tool_audit.py"]
            },
        },
    )
    print(
        json.dumps(
            {
                "environment_validated": True,
                "python": platform.python_version(),
                "wheel_python_members_matched": len(matches),
                "bwa": "0.7.17-r1188",
                "java_major": 17,
            }
        )
    )


if __name__ == "__main__":
    main()
