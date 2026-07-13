#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from _mudata_utils import write_table
from manuscript_object_specs import MANUSCRIPT_OBJECTS


def load_audit(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {str(row["dataset"]): row for row in data}


def dataset_inventory(audit: dict[str, dict[str, object]]) -> pd.DataFrame:
    rows = []
    for dataset, spec in MANUSCRIPT_OBJECTS.items():
        observed = audit.get(dataset, {})
        rows.append(
            {
                "dataset": dataset,
                "accession": spec["accession"],
                "modalities": ";".join(observed.get("modalities", spec["expected_shapes"].keys())),
                "shared_cells": observed.get("shared_cells", spec["expected_shared_cells"]),
                "source_data_type": spec["source_data_type"],
                "manuscript_role": spec["intended_role"],
                "object_path": observed.get("relative_path", spec["default_path"]),
                "sha256": observed.get("sha256", ""),
            }
        )
    return pd.DataFrame(rows)


def workflow_status() -> pd.DataFrame:
    rows = [
        {
            "component": "CIRI3 Smart-seq3",
            "status": "validated",
            "evidence_scope": "E-MTAB-8735 Smart-seq3 RNA+circRNA object and workflow outputs",
            "caveat": "CIRI3 is the primary manuscript-validation detector; circyto does not introduce a new back-splice-junction algorithm.",
        },
        {
            "component": "CIRI3 single-end RamDA/scRR",
            "status": "validated",
            "evidence_scope": "IMR90 single-end scRR RNA-side workflow",
            "caveat": "Genome-state CNV modality comes from processed public summary tables.",
        },
        {
            "component": "CIRI3 paired-end RamDA/scRR",
            "status": "validated with caveat",
            "evidence_scope": "HAP1 paired-end scRR route executable with explicit opt-in",
            "caveat": "Use --allow-paired-ramda; paired-end RamDA remains a documented caveat in the CLI.",
        },
        {
            "component": "IMR90 RNA+circRNA+CNV",
            "status": "validated",
            "evidence_scope": "23-cell trimodal MuData from public RNA data and processed public CNV summaries",
            "caveat": "CNV is not raw DNA FASTQ reprocessing.",
        },
        {
            "component": "HAP1 RNA+circRNA+RT",
            "status": "validated",
            "evidence_scope": "56 shared trimodal cells from public RNA data and processed public RT/state summaries",
            "caveat": "RT is not raw DNA FASTQ reprocessing.",
        },
        {
            "component": "SComatic interoperability",
            "status": "exploratory",
            "evidence_scope": "RNA-derived candidate variant signal import/merge path",
            "caveat": "Not orthogonally confirmed DNA variation.",
        },
        {
            "component": "find-circ3 adapter",
            "status": "experimental adapter",
            "evidence_scope": "Adapter-level extensibility tests",
            "caveat": "Not the primary manuscript validation detector.",
        },
        {
            "component": "CIRI2/CIRI-full/CIRCexplorer2 adapters",
            "status": "experimental adapter",
            "evidence_scope": "Extensibility and normalization paths",
            "caveat": "Do not present as manuscript-scale validation unless additional evidence is recovered.",
        },
    ]
    return pd.DataFrame(rows)


def host_gene_recovery(audit: dict[str, dict[str, object]]) -> pd.DataFrame:
    rows = []
    for dataset, spec in MANUSCRIPT_OBJECTS.items():
        observed = audit.get(dataset, {})
        annotated = int(observed.get("host_gene_annotated", spec["expected_host_gene_annotated"]))
        total = int(observed.get("host_gene_total", spec["expected_host_gene_total"]))
        rows.append(
            {
                "dataset": dataset,
                "annotated_circRNAs": annotated,
                "total_circRNAs": total,
                "recovery_fraction": annotated / total if total else 0.0,
                "host_gene_source_counts": json.dumps(observed.get("host_gene_source_counts", {}), sort_keys=True),
            }
        )
    return pd.DataFrame(rows)


def software_dependencies() -> pd.DataFrame:
    rows = [
        {
            "name": "Python",
            "role": "Package runtime",
            "required_for_basic_tests": "yes",
            "install_scope": "environment.yml or pip",
            "license_or_note": ">=3.10 per pyproject.toml",
        },
        {
            "name": "mudata",
            "role": "Read/write manuscript MuData objects",
            "required_for_basic_tests": "yes in CI dev extra",
            "install_scope": "circyto[mudata] or circyto[dev]",
            "license_or_note": "Optional dependency in pyproject.toml",
        },
        {
            "name": "scanpy",
            "role": "Figure 1 RNA UMAP data export",
            "required_for_basic_tests": "no",
            "install_scope": "circyto[scanpy]",
            "license_or_note": "Only needed for downstream UMAP export",
        },
        {
            "name": "CIRI3",
            "role": "Primary circRNA detector for manuscript validation",
            "required_for_basic_tests": "no",
            "install_scope": "External tool or configured jar path",
            "license_or_note": "Local bundled license is GPL-2.0; redistribution should be confirmed before release packaging.",
        },
        {
            "name": "STAR",
            "role": "Smart-seq3 and paired-end alignment route",
            "required_for_basic_tests": "no",
            "install_scope": "External executable",
            "license_or_note": "Workflow-specific dependency",
        },
        {
            "name": "bwa",
            "role": "CIRI3 BWA/direct and rescue alignment route",
            "required_for_basic_tests": "no",
            "install_scope": "External executable",
            "license_or_note": "Workflow-specific dependency",
        },
        {
            "name": "samtools",
            "role": "Alignment handling",
            "required_for_basic_tests": "no",
            "install_scope": "External executable",
            "license_or_note": "Workflow-specific dependency",
        },
    ]
    return pd.DataFrame(rows)


def reproducibility_commands() -> pd.DataFrame:
    rows = [
        {
            "artifact": "Object audit and manifest",
            "command": "python scripts/manuscript/audit_manuscript_objects.py",
            "expected_outputs": "manuscript/results/manuscript_object_audit.md; manuscript/results/manuscript_object_manifest.tsv",
        },
        {
            "artifact": "Object validation",
            "command": "python scripts/manuscript/validate_manuscript_objects.py --outdir manuscript/results/validation_check",
            "expected_outputs": "manuscript/results/validation_check/validation_summary.tsv; validation_summary.json",
        },
        {
            "artifact": "Smart-seq3 Figure 1 data",
            "command": "python scripts/manuscript/export_smartseq3_figure1_data.py --outdir manuscript/results",
            "expected_outputs": "smartseq3_umap_cells.tsv; smartseq3_top_circRNA_candidates.tsv; smartseq3_top_hostgene_features.tsv; smartseq3_selected_feature_candidates.tsv",
        },
        {
            "artifact": "Supplement-ready small tables",
            "command": "python scripts/manuscript/generate_supplement_tables.py --outdir manuscript/results",
            "expected_outputs": "manuscript_dataset_inventory.tsv; detector_workflow_status.tsv; host_gene_annotation_recovery.tsv; software_external_dependency_table.tsv; reproducibility_command_table.tsv",
        },
        {
            "artifact": "Repository large-file guard",
            "command": "python scripts/check_repository_file_sizes.py --base-ref origin/main --fail-scope added",
            "expected_outputs": "stdout report; nonzero exit for unapproved newly added large files",
        },
    ]
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate small supplement-ready manuscript tables.")
    parser.add_argument("--audit-json", type=Path, default=Path("manuscript/results/manuscript_object_audit.json"))
    parser.add_argument("--outdir", type=Path, default=Path("manuscript/results"))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    audit = load_audit(args.audit_json)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_table(dataset_inventory(audit), args.outdir / "manuscript_dataset_inventory.tsv")
    write_table(workflow_status(), args.outdir / "detector_workflow_status.tsv")
    write_table(host_gene_recovery(audit), args.outdir / "host_gene_annotation_recovery.tsv")
    write_table(software_dependencies(), args.outdir / "software_external_dependency_table.tsv")
    write_table(reproducibility_commands(), args.outdir / "reproducibility_command_table.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
