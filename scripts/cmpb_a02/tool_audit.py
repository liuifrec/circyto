#!/usr/bin/env python3
"""Transparent A02 subprocess audit wrapper. No production normalization/filtering.

The actual executable receives exactly the original argument vector and inherited
stdin/stdout/stderr. Nonzero/partial outputs remain evidence, not completion.
"""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

from run_stage import digest, dump, stamp


def option(argv, key):
    return Path(argv[argv.index(key) + 1]) if key in argv else None


def run(root, tool, executable, args):
    executable = Path(executable).resolve(strict=True)
    directory = root / "records/tool_commands"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (
        str(time.time_ns()) + "-" + str(os.getpid()) + "-" + tool + ".json"
    )
    inputs, output, kind = [], None, "runtime_inspection"
    if tool == "bwa" and args and args[0] in ["index", "mem"]:
        kind = "index" if args[0] == "index" else "alignment"
        inputs = (
            [Path(args[-1])] if kind == "index" else [Path(args[-2]), Path(args[-1])]
        )
        if kind == "alignment":
            target = Path(os.readlink("/proc/self/fd/1"))
            if target.is_absolute() and target.is_file():
                output = target
    if tool == "java" and "-jar" in args:
        inputs = [option(args, "-jar")]
        if "-I" in args:
            kind = "detector"
            inputs.extend(
                [option(args, flag) for flag in ["-I", "-F", "-A"] if flag in args]
            )
            output = option(args, "-O")
    inputs = [p.resolve(strict=True) for p in inputs]
    record = {
        "schema": "circyto.cmpb.a02.tool.v1",
        "kind": kind,
        "started_utc": stamp(),
        "tool": tool,
        "executable": str(executable),
        "executable_sha256": digest(executable),
        "audit_helper_sha256": digest(__file__),
        "audit_library_sha256": digest(Path(__file__).with_name("run_stage.py")),
        "argv": [str(executable), *args],
        "pid": os.getpid(),
        "input_sha256": {str(p): digest(p) for p in inputs},
        "output_path": str(output) if output is not None else None,
        "status": "running",
        "completion_bound": False,
        "environment": {
            k: os.environ[k]
            for k in [
                "JAVA_TOOL_OPTIONS",
                "CIRCYTO_CIRI3_EXTRA_ARGS",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
            ]
            if k in os.environ
        },
    }
    dump(path, record)
    started = time.monotonic()
    result = subprocess.run([str(executable), *args], check=False)
    unchanged = all(digest(p) == record["input_sha256"][str(p)] for p in inputs)
    output_sha = digest(output) if output is not None and output.is_file() else None
    record.update(
        {
            "finished_utc": stamp(),
            "elapsed_seconds": time.monotonic() - started,
            "exit_code": result.returncode,
            "inputs_unchanged": unchanged,
            "output_sha256": output_sha,
            "status": "exited",
            "completion_bound": result.returncode == 0
            and unchanged
            and output_sha is not None,
        }
    )
    dump(path, record)
    if not unchanged:
        return 125
    return result.returncode if result.returncode >= 0 else 128 - result.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--tool", required=True, choices=["bwa", "java"])
    parser.add_argument("--real", required=True)
    parser.add_argument("argv", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    argv = args.argv[1:] if args.argv[:1] == ["--"] else args.argv
    return run(args.root, args.tool, args.real, argv)


if __name__ == "__main__":
    sys.exit(main())
