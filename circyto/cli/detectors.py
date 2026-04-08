from __future__ import annotations

import json
import shutil
from typing import Dict, Any

import typer

from circyto.cli.doctor_meta import DETECTOR_SPECS, resolve_asset

detectors_app = typer.Typer(
    help="List available detectors and their dependency requirements.",
    invoke_without_command=True,
)


def _have_cmd(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _have_asset(asset: str) -> bool:
    return resolve_asset(asset).found


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
        asset_details = {}
        for asset in spec.required_assets:
            resolution = resolve_asset(asset)
            asset_details[asset] = {
                "found": resolution.found,
                "path": str(resolution.resolved_path) if resolution.resolved_path else None,
                "source": resolution.source,
                "checked_paths": [str(path) for path in resolution.checked_paths],
            }
        rows.append(
            {
                "name": spec.name,
                "status": status,
                "type": spec.det_type,
                "needs": needs,
                "assets": asset_details,
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
        for asset, details in r["assets"].items():
            if details["path"]:
                typer.echo(f"  asset: {asset} -> {details['path']}")
            elif details["checked_paths"]:
                typer.echo(f"  asset: {asset} -> missing")
                for checked in details["checked_paths"]:
                    typer.echo(f"    checked: {checked}")

    if hints:
        typer.echo("\nHints:")
        for r in rows:
            typer.echo(f"  {r['name']}:")
            for line in r["hints"]:
                typer.echo(f"    - {line}")

    raise typer.Exit(code=0)
