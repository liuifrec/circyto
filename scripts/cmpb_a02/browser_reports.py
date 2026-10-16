#!/usr/bin/env python3
"""Open the two actual A02 reports offline and compare browser DOM values."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from independent_checks import Tables
from run_stage import digest, dump, stamp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--workflow-directory", default="workflow_firstcell")
    ap.add_argument("--import-directory", default="import_firstcell")
    ap.add_argument("--audit-directory", default="records/acceptance")
    args = ap.parse_args()
    root = args.root.resolve()
    chrome = shutil.which("google-chrome")
    if chrome is None:
        raise RuntimeError("No local headless Chrome available for report opening")
    results = []
    for label, directory in [
        ("new_workflow", args.workflow_directory),
        ("native_import", args.import_directory),
    ]:
        source = root / directory / "qc/report.html"
        if not source.is_file():
            raise FileNotFoundError(source)
        dest = root / args.audit_directory / ("browser_" + label)
        dest.mkdir(parents=True, exist_ok=False)
        screenshot = dest / "report.png"
        dom_path = dest / "dom.html"
        argv = [
            chrome,
            "--headless",
            "--disable-gpu",
            "--disable-dev-shm-usage",
            "--disable-background-networking",
            "--disable-sync",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-extensions",
            "--disable-component-update",
            "--host-resolver-rules=MAP * ~NOTFOUND",
            "--user-data-dir=" + str(root / "tmp" / ("chrome-" + label)),
            "--window-size=1440,1200",
            "--screenshot=" + str(screenshot),
            "--dump-dom",
            source.as_uri(),
        ]
        with dom_path.open("x") as out, (dest / "stderr.log").open("x") as err:
            result = subprocess.run(argv, stdout=out, stderr=err, check=False)
        if result.returncode or not screenshot.is_file():
            raise RuntimeError("Offline browser failed: " + label)
        original = Tables()
        original.feed(source.read_text())
        rendered = Tables()
        rendered.feed(dom_path.read_text())
        if (
            original.tables != rendered.tables
            or original.cards != rendered.cards
            or original.svgs != rendered.svgs
        ):
            raise AssertionError(
                "Browser DOM changed displayed report values: " + label
            )
        results.append(
            {
                "report": label,
                "argv": argv,
                "exit_code": result.returncode,
                "source_sha256": digest(source),
                "DOM_sha256": digest(dom_path),
                "screenshot_sha256": digest(screenshot),
                "tables_equal": len(original.tables),
                "cards_equal": len(original.cards),
                "SVG_values_equal": len(original.svgs),
                "network": "DNS disabled; local file and inline assets only",
                "screenshot": str(screenshot.relative_to(root)),
            }
        )
        dump(
            root / args.audit_directory / "browser_acceptance.json",
            {
                "recorded_utc": stamp(),
                "results": results,
                "complete": len(results) == 2,
            },
        )
    print(
        json.dumps(
            {"offline_reports_opened": len(results), "DOM_values_preserved": True}
        )
    )


if __name__ == "__main__":
    main()
