#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from audit_manuscript_objects import collect_summary
from _mudata_utils import json_dumps, write_table
from manuscript_object_specs import EXPECTED_HOST_GENE_FIELDS, MANUSCRIPT_OBJECTS


def load_expectations(path: Path | None) -> dict[str, dict[str, object]]:
    if path is None:
        return {key: dict(value) for key, value in MANUSCRIPT_OBJECTS.items()}
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit("Expectations JSON must be an object keyed by dataset name.")
    return data


def check_equal(rows: list[dict[str, object]], dataset: str, metric: str, observed: Any, expected: Any) -> None:
    ok = observed == expected
    rows.append(
        {
            "dataset": dataset,
            "metric": metric,
            "observed": json_dumps(observed) if isinstance(observed, (dict, list)) else observed,
            "expected": json_dumps(expected) if isinstance(expected, (dict, list)) else expected,
            "status": "pass" if ok else "fail",
            "message": "" if ok else f"{metric} mismatch for {dataset}",
        }
    )


def check_float(
    rows: list[dict[str, object]],
    dataset: str,
    metric: str,
    observed: float,
    expected: float,
    *,
    tolerance: float = 1e-9,
) -> None:
    ok = abs(float(observed) - float(expected)) <= tolerance
    rows.append(
        {
            "dataset": dataset,
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "status": "pass" if ok else "fail",
            "message": "" if ok else f"{metric} mismatch for {dataset}",
        }
    )


def check_required_fields(rows: list[dict[str, object]], dataset: str, observed_missing: list[str]) -> None:
    ok = len(observed_missing) == 0
    rows.append(
        {
            "dataset": dataset,
            "metric": "required_host_gene_fields",
            "observed": ",".join(observed_missing) if observed_missing else "all_present",
            "expected": ",".join(EXPECTED_HOST_GENE_FIELDS),
            "status": "pass" if ok else "fail",
            "message": "" if ok else f"Missing host-gene provenance fields: {', '.join(observed_missing)}",
        }
    )


def check_local_private_metadata(rows: list[dict[str, object]], dataset: str, summary: dict[str, object]) -> None:
    count = int(summary["local_paths_detected"])
    ok = count == 0
    rows.append(
        {
            "dataset": dataset,
            "metric": "local_private_metadata_findings",
            "observed": count,
            "expected": 0,
            "status": "pass" if ok else "warn",
            "message": ""
            if ok
            else (
                "Source-object metadata contains local/private-looking strings. "
                "Public audit outputs redact values; inspect manuscript_object_audit.md before release."
            ),
        }
    )


def validate_one(dataset: str, path: Path, expectation: dict[str, object]) -> tuple[dict[str, object], list[dict[str, object]]]:
    summary = collect_summary(dataset, path)
    rows: list[dict[str, object]] = []

    expected_shapes = expectation.get("expected_shapes", {})
    if not isinstance(expected_shapes, dict):
        raise SystemExit(f"{dataset} expected_shapes must be an object.")
    for modality, shape in expected_shapes.items():
        observed_shape = summary["modality_shapes"].get(modality)
        check_equal(rows, dataset, f"{modality}_shape", observed_shape, shape)

    check_equal(rows, dataset, "shared_cells", summary["shared_cells"], expectation.get("expected_shared_cells"))
    check_equal(
        rows,
        dataset,
        "host_gene_annotated",
        summary["host_gene_annotated"],
        expectation.get("expected_host_gene_annotated"),
    )
    check_equal(
        rows,
        dataset,
        "host_gene_total",
        summary["host_gene_total"],
        expectation.get("expected_host_gene_total"),
    )
    expected_total = expectation.get("expected_host_gene_total")
    expected_annotated = expectation.get("expected_host_gene_annotated")
    if expected_total:
        expected_rate = float(expected_annotated) / float(expected_total)
        check_float(rows, dataset, "host_gene_rate", float(summary["host_gene_rate"]), expected_rate, tolerance=5e-4)
    if "expected_median_circRNA_count" in expectation:
        check_float(
            rows,
            dataset,
            "median_circRNA_count",
            float(summary["median_circRNA_count"]),
            float(expectation["expected_median_circRNA_count"]),
        )
    if "expected_median_circRNA_total_support" in expectation:
        check_float(
            rows,
            dataset,
            "median_circRNA_total_support",
            float(summary["median_circRNA_total_support"]),
            float(expectation["expected_median_circRNA_total_support"]),
        )

    check_required_fields(rows, dataset, list(summary["missing_host_gene_fields"]))
    check_local_private_metadata(rows, dataset, summary)
    public_ok = str(summary["public_data_only"]).startswith("yes_")
    rows.append(
        {
            "dataset": dataset,
            "metric": "public_cell_identifier_assessment",
            "observed": summary["public_data_only"],
            "expected": "yes_public_accession_or_canonical_cell_ids",
            "status": "pass" if public_ok else "fail",
            "message": "" if public_ok else "Cell identifiers require manual review before public release.",
        }
    )

    return summary, rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate manuscript-scale MuData objects against current manuscript claims without rerunning detectors."
    )
    parser.add_argument("--smartseq3", type=Path, default=MANUSCRIPT_OBJECTS["Smart-seq3"]["default_path"])
    parser.add_argument("--hap1", type=Path, default=MANUSCRIPT_OBJECTS["HAP1"]["default_path"])
    parser.add_argument("--imr90", type=Path, default=MANUSCRIPT_OBJECTS["IMR90"]["default_path"])
    parser.add_argument("--outdir", type=Path, default=Path("manuscript/results/validation_check"))
    parser.add_argument(
        "--expectations-json",
        type=Path,
        default=None,
        help="Optional JSON expectations for tests or future manuscript revisions.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    expectations = load_expectations(args.expectations_json)
    inputs = {
        "Smart-seq3": args.smartseq3,
        "HAP1": args.hap1,
        "IMR90": args.imr90,
    }
    all_summaries: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    for dataset, path in inputs.items():
        if dataset not in expectations:
            continue
        summary, rows = validate_one(dataset, path, expectations[dataset])
        all_summaries.append(summary)
        all_rows.extend(rows)

    args.outdir.mkdir(parents=True, exist_ok=True)
    validation = pd.DataFrame(all_rows)
    write_table(validation, args.outdir / "validation_summary.tsv")
    failures = validation[validation["status"] == "fail"].copy()
    warnings = validation[validation["status"] == "warn"].copy()
    write_table(failures, args.outdir / "validation_failures.tsv")
    (args.outdir / "validation_summary.json").write_text(
        json.dumps(
            {
                "status": "pass" if failures.empty else "fail",
                "n_checks": int(len(validation)),
                "n_failures": int(len(failures)),
                "n_warnings": int(len(warnings)),
                "failures": failures.to_dict(orient="records"),
                "warnings": warnings.to_dict(orient="records"),
                "objects": [
                    {
                        "dataset": s["dataset"],
                        "relative_path": s["relative_path"],
                        "sha256": s["sha256"],
                        "file_size_bytes": s["file_size_bytes"],
                        "modalities": s["modalities"],
                        "modality_shapes": s["modality_shapes"],
                        "shared_cells": s["shared_cells"],
                    }
                    for s in all_summaries
                ],
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    if not failures.empty:
        print(f"[ERROR] Manuscript object validation failed with {len(failures)} failing checks.")
        for row in failures.to_dict(orient="records"):
            observed = row["observed"]
            if isinstance(observed, float) and math.isnan(observed):
                observed = "nan"
            print(f"  - {row['dataset']} {row['metric']}: observed={observed} expected={row['expected']}")
        return 1
    if not warnings.empty:
        print(f"[OK] Manuscript object validation passed ({len(validation)} checks, {len(warnings)} warnings).")
        return 0
    print(f"[OK] Manuscript object validation passed ({len(validation)} checks).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
