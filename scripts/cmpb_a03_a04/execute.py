#!/usr/bin/env python3
"""Sequential analysis orchestration, never alignment/detection or package installation."""

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

from oracle import verify, truth, sha
from run_stage import dump, stamp, footprint


def command(config_path, arm, out):
    config = json.loads(config_path.read_text())
    if arm == "circyto":
        argv = [
            str(Path(sys.executable).with_name("circyto")),
            "import-circ",
            {
                "native": "--ciri3-manifest",
                "calls": "--calls",
                "matrix": "--matrix-dir",
            }[config["mode"]],
            config["input"],
            "--cells",
            config["cells"],
        ]
        for key in [
            "genome_build",
            "coordinate_convention",
            "support_definition",
            "detector",
        ]:
            argv += ["--" + key.replace("_", "-"), config[key]]
        if config.get("rna"):
            argv += ["--rna", config["rna"]]
        return argv + ["--outdir", str(out)]
    return [
        sys.executable,
        "-B",
        str(
            Path(__file__).with_name(
                "reference.py" if arm == "reference" else "unsafe.py"
            )
        ),
        "--config",
        str(config_path),
        "--outdir",
        str(out),
    ]


def run_arm(config_path, arm, directory):
    directory.mkdir(parents=True, exist_ok=False)
    out = directory / "output"
    argv = command(config_path, arm, out)
    timing = directory / "time.tsv"
    launched = ["/usr/bin/time", "-f", "%e\t%M\t%U\t%S", "-o", str(timing), *argv]
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    start = time.perf_counter()
    with (
        (directory / "stdout.log").open("x") as stdout,
        (directory / "stderr.log").open("x") as stderr,
    ):
        p = subprocess.run(launched, stdout=stdout, stderr=stderr, env=env, check=False)
    elapsed = time.perf_counter() - start
    rows = timing.read_text().splitlines()
    t = rows[-1].split("\t")
    result = dict(
        arm=arm,
        argv=argv,
        started_environment={
            k: env.get(k)
            for k in [
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            ]
        },
        affinity=sorted(os.sched_getaffinity(0)),
        elapsed_seconds=elapsed,
        time_wall_seconds=float(t[0]),
        peak_rss_kib=int(t[1]),
        user_seconds=float(t[2]),
        system_seconds=float(t[3]),
        exit_code=p.returncode,
        output_bytes=footprint(out) if out.exists() else 0,
        stdout_sha256=sha(directory / "stdout.log"),
        stderr_sha256=sha(directory / "stderr.log"),
    )
    dump(directory / "execution.json", result)
    return result, out


def freeze(root, label):
    sources = {str(p): sha(p) for p in Path(__file__).parent.glob("*.py")}
    sources.update(
        {
            str(p): sha(p)
            for p in (Path(__file__).parents[2] / "manuscript/cmpb/a03_a04").glob(
                "*.md"
            )
        }
    )
    sources[str(Path(__file__).parents[1] / "cmpb_a02/independent_checks.py")] = sha(
        Path(__file__).parents[1] / "cmpb_a02/independent_checks.py"
    )
    path = root / "records" / f"{label}.json"
    if path.exists():
        old = json.loads(path.read_text())
        if any(not Path(p).exists() or sha(p) != h for p, h in old["sources"].items()):
            raise ValueError("Frozen source changed")
        return old
    record = dict(recorded_utc=stamp(), sources=sources)
    dump(path, record)
    return record


def a03(root, action):
    C = root.parent / "continuation_20261016"
    folder = root / "a03"
    config_path = folder / "config.json"
    if not config_path.exists():
        dump(
            config_path,
            dict(
                mode="native",
                input=str(C / "inputs/native_four.tsv"),
                cells=str(C / "inputs/cells_four.tsv"),
                rna=str(C / "inputs/rna/imr90_original_overlap_profile.four.h5ad"),
                genome_build="GRCh38",
                coordinate_convention="1-based-closed",
                support_definition="CIRI3_junction_reads",
                detector="CIRI3",
            ),
        )
    config = json.loads(config_path.read_text())
    if action == "control":
        freeze(root, "reference_precontrol_sources")
        order = [("control", "circyto"), ("control", "reference")]
    else:
        gate = json.loads((root / "records/a03_control.json").read_text())
        if not gate["passed"]:
            raise ValueError("Real-input control failed")
        freeze(root, "reference_precontrol_sources")
        freeze(root, "frozen_a03_sources")
        order = [
            ("warmup", "circyto"),
            ("warmup", "reference"),
            ("measurement1", "circyto"),
            ("measurement1", "reference"),
            ("measurement2", "reference"),
            ("measurement2", "circyto"),
            ("measurement3", "circyto"),
            ("measurement3", "reference"),
        ]
    results = []
    for phase, arm in order:
        r, out = run_arm(config_path, arm, folder / (phase + "-" + arm))
        r["phase"] = phase
        r["measured"] = phase.startswith("measurement")
        if r["exit_code"] == 0:
            r["acceptance"], prevalence = verify(config, out, arm)
            r["input_bytes"] = sum(p.stat().st_size for p in truth(config)[4].values())
        else:
            r["acceptance"] = dict(passed=False, error="process exited nonzero")
        results.append(r)
        dump(
            root / "records" / f"a03_{action}.json",
            dict(passed=all(x["acceptance"]["passed"] for x in results), runs=results),
        )
        if not r["acceptance"]["passed"]:
            raise RuntimeError(
                "STOP: valid real-input arm failed; retain output and inspect " + arm
            )
    print(json.dumps(dict(action=action, runs=len(results), passed=True)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["control", "measure"])
    ap.add_argument("--root", type=Path, required=True)
    a = ap.parse_args()
    a03(a.root.resolve(), a.action)
