#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from _mudata_utils import (
    CIRC_MODALITY_CANDIDATES,
    common_obs_names,
    find_modality,
    get_or_compute_circ_metric,
    host_gene_recovery,
    host_gene_source_counts,
    json_dumps,
    load_mudata,
    write_table,
)
from manuscript_object_specs import EXPECTED_HOST_GENE_FIELDS, MANUSCRIPT_OBJECTS


PATH_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("unix_home_path", re.compile(r"(^|[\s=:;,])/(home|Users)/[^ \t\n\r;,\"]+")),
    ("wsl_or_mount_path", re.compile(r"(^|[\s=:;,])/(mnt|media)/[^ \t\n\r;,\"]+")),
    ("server_path", re.compile(r"(^|[\s=:;,])/(user|users|nfs|gpfs|lustre|scratch|work|data)/[^ \t\n\r;,\"]+")),
    ("tmp_path", re.compile(r"(^|[\s=:;,])/tmp/[^ \t\n\r;,\"]+")),
    ("windows_path", re.compile(r"[A-Za-z]:\\[^ \t\n\r;,\"]+")),
    ("file_uri", re.compile(r"file://[^ \t\n\r;,\"]+")),
    ("known_local_username", re.compile(r"\b(yul03|liuyuchen|ifrec)\b", re.IGNORECASE)),
]


def redact_private_value(category: str, value: str) -> str:
    if category == "known_local_username":
        return "[local username redacted]"
    return "[local path redacted]"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def matrix_summary(adata) -> dict[str, object]:
    matrix = adata.X
    if sparse.issparse(matrix):
        nnz = int(matrix.nnz)
        storage = type(matrix).__name__
    else:
        arr = np.asarray(matrix)
        nnz = int(np.count_nonzero(arr))
        storage = type(matrix).__name__
    total = int(adata.n_obs * adata.n_vars)
    return {
        "storage": storage,
        "nnz": nnz,
        "density": float(nnz / total) if total else 0.0,
    }


def list_keys(mapping: Any) -> list[str]:
    try:
        return [str(key) for key in mapping.keys()]
    except Exception:
        return []


def scalar_to_text(value: Any) -> str | None:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8", errors="replace")
    if isinstance(value, np.str_):
        return str(value)
    return None


def array_strings(values: Any, *, limit: int = 200_000) -> list[str]:
    arr = np.asarray(values)
    if arr.size > limit:
        return []
    if arr.dtype.kind not in {"O", "S", "U"}:
        return []
    out: list[str] = []
    for value in arr.ravel():
        text = scalar_to_text(value)
        if text is not None:
            out.append(text)
    return out


def detect_private_strings(path: Path) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []

    def inspect_text(location: str, text: str) -> None:
        if not text:
            return
        for category, pattern in PATH_PATTERNS:
            for match in pattern.finditer(text):
                value = match.group(0).strip(" \t\n\r=:;,")
                findings.append(
                    {
                        "location": location,
                        "category": category,
                        "value": redact_private_value(category, value),
                        "_raw_value": value,
                    }
                )

    def visitor(name: str, obj: Any) -> None:
        for attr_key, attr_value in obj.attrs.items():
            text = scalar_to_text(attr_value)
            if text is not None:
                inspect_text(f"{name}.attrs[{attr_key}]", text)
            else:
                for item in array_strings(attr_value):
                    inspect_text(f"{name}.attrs[{attr_key}]", item)
        if isinstance(obj, h5py.Dataset):
            if h5py.check_string_dtype(obj.dtype) is not None or obj.dtype.kind in {"O", "S", "U"}:
                try:
                    values = obj.asstr()[()] if h5py.check_string_dtype(obj.dtype) is not None else obj[()]
                except Exception:
                    return
                for item in array_strings(values):
                    inspect_text(name, item)

    with h5py.File(path, "r") as handle:
        handle.visititems(visitor)

    deduped: dict[tuple[str, str, str], dict[str, str]] = {}
    for finding in findings:
        raw_value = finding.get("_raw_value", finding["value"])
        key = (finding["location"], finding["category"], raw_value)
        deduped[key] = {
            "location": finding["location"],
            "category": finding["category"],
            "value": finding["value"],
        }
    return list(deduped.values())


def compression_summary(path: Path) -> dict[str, object]:
    total = 0
    compressed = 0
    compression_counter: Counter[str] = Counter()
    matrix_datasets: list[dict[str, object]] = []

    def visitor(name: str, obj: Any) -> None:
        nonlocal total, compressed
        if not isinstance(obj, h5py.Dataset):
            return
        total += 1
        if obj.compression:
            compressed += 1
            compression_counter[str(obj.compression)] += 1
        if name.endswith(("/X/data", "/X", "/layers/mappabilitynorm/data", "/layers/mappabilitynorm")):
            matrix_datasets.append(
                {
                    "dataset": name,
                    "shape": "x".join(map(str, obj.shape)),
                    "dtype": str(obj.dtype),
                    "compression": str(obj.compression or ""),
                    "chunks": "x".join(map(str, obj.chunks)) if obj.chunks else "",
                }
            )

    with h5py.File(path, "r") as handle:
        handle.visititems(visitor)
    return {
        "datasets": total,
        "compressed_datasets": compressed,
        "compression_algorithms": dict(sorted(compression_counter.items())),
        "matrix_datasets": matrix_datasets,
    }


def describe_modalities(mdata) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for name, adata in mdata.mod.items():
        out[str(name)] = {
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "obs_columns": list(map(str, adata.obs.columns)),
            "var_columns": list(map(str, adata.var.columns)),
            "obsm_keys": list_keys(adata.obsm),
            "varm_keys": list_keys(adata.varm),
            "layers": list_keys(adata.layers),
            "uns_keys": list_keys(adata.uns),
            "matrix": matrix_summary(adata),
            "obs_name_examples": list(map(str, adata.obs_names[:5])),
            "var_name_examples": list(map(str, adata.var_names[:5])),
        }
    return out


def public_cell_id_assessment(mdata, findings: list[dict[str, str]]) -> str:
    cell_values: list[str] = []
    for adata in mdata.mod.values():
        cell_values.extend(map(str, adata.obs_names[:20]))
    cell_text = "\n".join(cell_values)
    suspicious = [
        finding
        for finding in findings
        if "/obs/" in finding["location"] or finding["location"].endswith("/obs/_index")
    ]
    if suspicious:
        return "review_required"
    if any(pattern.search(cell_text) for _, pattern in PATH_PATTERNS):
        return "review_required"
    return "yes_public_accession_or_canonical_cell_ids"


def collect_summary(dataset: str, path: Path) -> dict[str, object]:
    spec = MANUSCRIPT_OBJECTS[dataset]
    if not path.exists():
        raise FileNotFoundError(path)
    mdata = load_mudata(path)
    modalities = list(map(str, mdata.mod.keys()))
    modality_shapes = {name: [int(mdata.mod[name].n_obs), int(mdata.mod[name].n_vars)] for name in modalities}
    shared_cells = len(common_obs_names(*(mdata.mod[name] for name in modalities))) if modalities else 0
    pairwise_shared: dict[str, int] = {}
    for idx, left in enumerate(modalities):
        for right in modalities[idx + 1 :]:
            pairwise_shared[f"{left}&{right}"] = len(common_obs_names(mdata.mod[left], mdata.mod[right]))

    circ_name = find_modality(mdata, CIRC_MODALITY_CANDIDATES, required=False)
    host_annotated = host_total = 0
    host_rate = 0.0
    host_sources: dict[str, int] = {}
    median_count = math.nan
    median_support = math.nan
    missing_host_fields: list[str] = []
    if circ_name is not None:
        circ = mdata.mod[circ_name]
        host_annotated, host_total, host_rate = host_gene_recovery(circ)
        host_sources = host_gene_source_counts(circ)
        median_count = float(get_or_compute_circ_metric(mdata, circ, "circRNA_count").median())
        median_support = float(get_or_compute_circ_metric(mdata, circ, "circRNA_total_support").median())
        missing_host_fields = [field for field in EXPECTED_HOST_GENE_FIELDS if field not in circ.var.columns]

    findings = detect_private_strings(path)
    compression = compression_summary(path)
    modalities_detail = describe_modalities(mdata)

    return {
        "dataset": dataset,
        "accession": spec["accession"],
        "relative_path": str(path),
        "file_size_bytes": int(path.stat().st_size),
        "sha256": sha256_file(path),
        "modalities": modalities,
        "modality_shapes": modality_shapes,
        "shared_cells": int(shared_cells),
        "pairwise_shared_cells": pairwise_shared,
        "source_data_type": spec["source_data_type"],
        "public_data_only": public_cell_id_assessment(mdata, findings),
        "local_paths_detected": int(len(findings)),
        "local_path_findings": findings,
        "intended_role": spec["intended_role"],
        "manuscript_release_status": spec["manuscript_release_status"],
        "mdata_obs_columns": list(map(str, mdata.obs.columns)),
        "mdata_obsm_keys": list_keys(mdata.obsm),
        "mdata_varm_keys": list_keys(mdata.varm),
        "mdata_uns_keys": list_keys(mdata.uns),
        "modalities_detail": modalities_detail,
        "host_gene_annotated": int(host_annotated),
        "host_gene_total": int(host_total),
        "host_gene_rate": float(host_rate),
        "host_gene_source_counts": host_sources,
        "missing_host_gene_fields": missing_host_fields,
        "median_circRNA_count": median_count,
        "median_circRNA_total_support": median_support,
        "compression": compression,
    }


def manifest_rows(summaries: list[dict[str, object]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for summary in summaries:
        rows.append(
            {
                "dataset": summary["dataset"],
                "accession": summary["accession"],
                "relative_path": summary["relative_path"],
                "file_size_bytes": summary["file_size_bytes"],
                "sha256": summary["sha256"],
                "modalities": ";".join(summary["modalities"]),
                "modality_shapes": json_dumps(summary["modality_shapes"]),
                "shared_cells": summary["shared_cells"],
                "source_data_type": summary["source_data_type"],
                "public_data_only": summary["public_data_only"],
                "local_paths_detected": summary["local_paths_detected"],
                "intended_role": summary["intended_role"],
                "manuscript_release_status": summary["manuscript_release_status"],
            }
        )
    return pd.DataFrame(rows)


def format_rate(value: float) -> str:
    return f"{value:.1%}"


def md_table(headers: list[str], rows: list[list[object]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(out)


def build_markdown(summaries: list[dict[str, object]]) -> str:
    lines: list[str] = [
        "# Manuscript Object Audit",
        "",
        "Read-only audit of the committed manuscript-scale MuData objects. No matrices or metadata were rewritten.",
        "",
        "## Summary",
        "",
        md_table(
            [
                "Dataset",
                "Path",
                "Size bytes",
                "SHA-256",
                "Modalities",
                "Shapes",
                "Shared cells",
                "Host-gene recovery",
                "Local/private findings",
            ],
            [
                [
                    s["dataset"],
                    s["relative_path"],
                    s["file_size_bytes"],
                    s["sha256"],
                    ";".join(s["modalities"]),
                    json_dumps(s["modality_shapes"]),
                    s["shared_cells"],
                    f"{s['host_gene_annotated']}/{s['host_gene_total']} ({format_rate(float(s['host_gene_rate']))})",
                    s["local_paths_detected"],
                ]
                for s in summaries
            ],
        ),
        "",
    ]

    for summary in summaries:
        lines.extend(
            [
                f"## {summary['dataset']}",
                "",
                f"- Accession/source: `{summary['accession']}`.",
                f"- Role: {summary['intended_role']}",
                f"- Public-data assessment: `{summary['public_data_only']}`.",
                f"- Pairwise shared cells: `{json_dumps(summary['pairwise_shared_cells'])}`.",
                f"- Median detected circRNAs per cell: `{summary['median_circRNA_count']}`.",
                f"- Median circRNA support per cell: `{summary['median_circRNA_total_support']}`.",
                f"- Host-gene source counts: `{json_dumps(summary['host_gene_source_counts'])}`.",
                f"- Missing host-gene provenance fields: `{', '.join(summary['missing_host_gene_fields']) or 'none'}`.",
                f"- Top-level obs columns: `{', '.join(summary['mdata_obs_columns']) or 'none'}`.",
                f"- Top-level obsm keys: `{', '.join(summary['mdata_obsm_keys']) or 'none'}`.",
                f"- Top-level varm keys: `{', '.join(summary['mdata_varm_keys']) or 'none'}`.",
                f"- Top-level uns keys: `{', '.join(summary['mdata_uns_keys']) or 'none'}`.",
                "",
                "### Modalities",
                "",
            ]
        )
        for name, detail in summary["modalities_detail"].items():
            matrix = detail["matrix"]
            lines.extend(
                [
                    f"- `{name}` shape `{detail['shape'][0]} x {detail['shape'][1]}`; "
                    f"X storage `{matrix['storage']}`, nnz `{matrix['nnz']}`, density `{matrix['density']:.6g}`.",
                    f"  - obs columns: `{', '.join(detail['obs_columns']) or 'none'}`.",
                    f"  - var columns: `{', '.join(detail['var_columns']) or 'none'}`.",
                    f"  - obsm keys: `{', '.join(detail['obsm_keys']) or 'none'}`; "
                    f"varm keys: `{', '.join(detail['varm_keys']) or 'none'}`; "
                    f"layers: `{', '.join(detail['layers']) or 'none'}`; "
                    f"uns keys: `{', '.join(detail['uns_keys']) or 'none'}`.",
                    f"  - example cells: `{', '.join(detail['obs_name_examples'])}`.",
                    f"  - example features: `{', '.join(detail['var_name_examples'])}`.",
                ]
            )
        comp = summary["compression"]
        lines.extend(
            [
                "",
                "### Compression And Stored Matrices",
                "",
                f"- HDF5 datasets: `{comp['datasets']}`; compressed datasets: `{comp['compressed_datasets']}`; "
                f"algorithms: `{json_dumps(comp['compression_algorithms'])}`.",
                f"- Matrix datasets: `{json_dumps(comp['matrix_datasets'])}`.",
                "- Manuscript necessity assessment: RNA `X`, circRNA `X`, and RT/CNV `X` are directly needed for manuscript validation. "
                "Any CNV `mappabilitynorm` layer is useful provenance for processed public CNV summaries. No object rewrite was performed.",
                "",
                "### Local/Private Metadata Scan",
                "",
            ]
        )
        findings = summary["local_path_findings"]
        if findings:
            lines.append(
                md_table(
                    ["Location", "Category", "Value"],
                    [[item["location"], item["category"], item["value"]] for item in findings[:100]],
                )
            )
            if len(findings) > 100:
                lines.append(f"\nOnly the first 100 of {len(findings)} findings are shown.")
        else:
            lines.append("No absolute local paths, local usernames, server paths, or private-looking sample identifiers were detected by the configured scan.")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Read-only audit of manuscript-scale MuData objects.")
    parser.add_argument("--smartseq3", type=Path, default=MANUSCRIPT_OBJECTS["Smart-seq3"]["default_path"])
    parser.add_argument("--hap1", type=Path, default=MANUSCRIPT_OBJECTS["HAP1"]["default_path"])
    parser.add_argument("--imr90", type=Path, default=MANUSCRIPT_OBJECTS["IMR90"]["default_path"])
    parser.add_argument("--manifest-out", type=Path, default=Path("manuscript/results/manuscript_object_manifest.tsv"))
    parser.add_argument("--report-out", type=Path, default=Path("manuscript/results/manuscript_object_audit.md"))
    parser.add_argument("--json-out", type=Path, default=Path("manuscript/results/manuscript_object_audit.json"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    inputs = {
        "Smart-seq3": args.smartseq3,
        "HAP1": args.hap1,
        "IMR90": args.imr90,
    }
    summaries = [collect_summary(dataset, path) for dataset, path in inputs.items()]
    write_table(manifest_rows(summaries), args.manifest_out)
    args.report_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.write_text(build_markdown(summaries), encoding="utf-8")
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summaries, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
