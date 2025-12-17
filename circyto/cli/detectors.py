from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Dict, Any

import typer

from circyto.cli.doctor_meta import DETECTOR_SPECS

detectors_app = typer.Typer(
    help="List available detectors and their dependency requirements.",
    invoke_without_command=True,
)


def _have_cmd(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _have_asset(asset: str) -> bool:
    # Keep this consistent with doctor.py
    if asset == "CIRI-full-jar":
        tools_dir = Path(__file__).resolve().parents[2] / "tools"
        if not tools_dir.exists():
            return False
        return len(list(tools_dir.glob("CIRI-full*.jar"))) > 0
    return False


def _status_for(spec) -> str:
    missing_cmds = [c for c in spec.required_cmds if not _have_cmd(c)]
    missing_assets = [a for a in spec.required_assets if not _have_asset(a)]

    if not missing_cmds and not missing_assets:
        return "READY"
    if missing_cmds:
        return "NOT READY"
    return "PARTIAL"


@detectors_app.callback()
def detectors(
    json_out: bool = typer.Option(False, "--json", help="Output machine-readable JSON."),
    hints: bool = typer.Option(False, "--hints", help="Print install hints."),
):
    """
    Print detectors and dependency readiness.
    """
    rows = []
    for spec in DETECTOR_SPECS:
        status = _status_for(spec)
        needs = spec.required_cmds + spec.required_assets
        rows.append(
            {
                "name": spec.name,
                "status": status,
                "type": spec.det_type,
                "needs": needs,
                "hints": spec.hint_lines,
            }
        )

    if json_out:
        typer.echo(json.dumps({"detectors": rows}, indent=2))
        raise typer.Exit(code=0)

    # Pretty table (simple, no extra deps)
    typer.echo("NAME         STATUS     TYPE   NEEDS")
    for r in rows:
        needs_str = ", ".join(r["needs"]) if r["needs"] else "-"
        typer.echo(f"{r['name']:<12} {r['status']:<10} {r['type']:<5} {needs_str}")

    if hints:
        typer.echo("\nHints:")
        for r in rows:
            typer.echo(f"  {r['name']}:")
            for line in r["hints"]:
                typer.echo(f"    - {line}")

    raise typer.Exit(code=0)
