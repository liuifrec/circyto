from __future__ import annotations

import sys
import shutil
from pathlib import Path
from typing import List
from circyto.cli.doctor_meta import DETECTOR_SPECS

import typer

DETECTORS = {
    "find-circ3": {
        "required_cmds": ["bowtie2", "samtools"],
        "required_assets": [],
    },
    "ciri-full": {
        "required_cmds": ["bwa", "java"],
        "required_assets": ["CIRI-full-jar"],
    },
}

doctor_app = typer.Typer(
    help="Diagnose circyto environment and external dependencies.",
    invoke_without_command=True,
)


def _print(status: str, name: str, msg: str):
    typer.echo(f"[{status:<4}] {name:<10} {msg}")


@doctor_app.callback()
def doctor():
    """
    Run environment diagnostics.
    """
    typer.echo("circyto doctor\n")

    fatal_missing: List[str] = []

    # ---- Python -------------------------------------------------------------
    ver = sys.version_info
    py_ver = f"{ver.major}.{ver.minor}.{ver.micro}"
    if (ver.major, ver.minor) >= (3, 10):
        _print("OK", "python", py_ver)
    else:
        _print("MISS", "python", f"{py_ver} (>=3.10 required)")
        fatal_missing.append("python>=3.10")
    found_cmds = {}
    found_assets = {}

    # ---- Required executables ----------------------------------------------
    def check(cmd: str, required_for: str | None = None):
        path = shutil.which(cmd)
        found_cmds[cmd] = path is not None
        if path:
            _print("OK", cmd, path)
        else:
            if required_for:
                _print("MISS", cmd, f"required for {required_for}")
                fatal_missing.append(cmd)
            else:
                _print("WARN", cmd, "not found (optional)")


    check("bowtie2", "find-circ3")
    check("samtools", "find-circ3")
    check("bwa", "ciri-full")
    check("java", "ciri-full")
    check("STAR")  # optional

    # ---- Detector assets ----------------------------------------------------
    tools_dir = Path(__file__).resolve().parents[2] / "tools"
    jars = list(tools_dir.glob("CIRI-full*.jar")) if tools_dir.exists() else []
    found_assets["CIRI-full-jar"] = bool(jars)

    if jars:
        _print("OK", "CIRI-full", "jar found in tools/")
    else:
        _print("WARN", "CIRI-full", "jar not found in tools/")

    # ---- Summary ------------------------------------------------------------
    typer.echo("\nDetector readiness:")
    for spec in DETECTOR_SPECS:
        missing_cmds = [c for c in spec.required_cmds if not found_cmds.get(c, False)]
        missing_assets = [a for a in spec.required_assets if not found_assets.get(a, False)]

        if not missing_cmds and not missing_assets:
            status = "READY"
            msg = ""
        elif missing_cmds:
            status = "NOT READY"
            msg = f"(missing: {', '.join(missing_cmds)})"
        else:
            status = "PARTIAL"
            msg = f"(missing: {', '.join(missing_assets)})"

        typer.echo(f"  {spec.name:<10}: {status}" + (f" {msg}" if msg else ""))


    typer.echo("\nSummary:")
    if fatal_missing:
        typer.echo("  ❌ Environment NOT ready")
        typer.echo("  Missing:")
        for m in sorted(set(fatal_missing)):
            typer.echo(f"   - {m}")
        raise typer.Exit(code=1)
    else:
        typer.echo("  ✅ Environment looks good")
        raise typer.Exit(code=0)

