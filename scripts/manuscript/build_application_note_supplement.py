#!/usr/bin/env python
"""Render small schema/workflow figures and Tables S1-S3 from audited evidence."""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from regenerate_application_note_results import (
    REPO_ROOT, package_versions, repository_provenance, sha256_file,
    validate_output_directory, write_json,
)

import pandas as pd


EVIDENCE = "manuscript/application_note_evidence.md"
WORKFLOWS = "docs/validated_workflows_summary.md"
BASELINE = "docs/manuscript_software_baseline.md"
STATUS = "docs/current_project_status.md"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    results = args.results.resolve()
    provenance = json.loads((results / "provenance.json").read_text())
    for filename, expected in provenance["generated_files"].items():
        path = (results / filename).resolve()
        if path.parent != results or sha256_file(path) != expected:
            raise SystemExit(f"Result checksum mismatch: {filename}")
    summary = json.loads((results / "regeneration_summary.json").read_text())
    smart, imr = summary["smartseq3"], summary["imr90"]
    output = validate_output_directory(args.output, [results])
    output.mkdir(parents=True)
    os.environ["MPLCONFIGDIR"] = str(output / ".cache" / "matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.hashsalt": "circyto-application-note-v2",
                                "font.family": "DejaVu Sans", "font.size": 10})
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    def save(fig, stem):
        for suffix in ("svg", "png"):
            fig.savefig(output / f"{stem}.{suffix}", dpi=300,
                        **({"metadata": {"Date": None}} if suffix == "svg" else {}))
        plt.close(fig)

    inventory = pd.DataFrame([
        {"dataset": "Smart-seq3", "accession": "E-MTAB-8735", "protocol": "Smart-seq3",
         "read_layout": "paired-end RNA + index reads", "route": "demultiplex; STAR; BWA rescue; CIRI3 STAR tuple",
         "reference": "hg38", "validation_cells": smart["cell_count"], "output_object": "RNA+circ MuData; circ AnnData",
         "circRNA_candidates": smart["circRNA_candidate_count"], "status": "real workflow; object regenerated",
         "role": "principal demonstration", "limitation": "detector candidates; no orthogonal circRNA validation", "source": WORKFLOWS + "; results_v2/regeneration_summary.json"},
        {"dataset": "IMR90", "accession": "GSE278958", "protocol": "scRR / RamDA-like",
         "read_layout": "single-end RNA; processed GEO CNV", "route": "BWA-MEM; CIRI3 direct SAM; cell remap; processed CNV merge",
         "reference": "hg38; processed 50-kb CNV bins", "validation_cells": imr["rna_shape"][0], "output_object": "RNA+circ+CNV MuData",
         "circRNA_candidates": imr["circRNA_shape"][1], "status": "real workflow; object regenerated",
         "role": "bounded multimodal interoperability", "limitation": "processed CNV only; no CNV biology inference", "source": WORKFLOWS + "; results_v2/regeneration_summary.json"},
        {"dataset": "HAP1", "accession": "GSE278952", "protocol": "scRR / RamDA-like",
         "read_layout": "paired-end RNA", "route": "STAR; BWA rescue; CIRI3 STAR tuple; --allow-paired-ramda",
         "reference": "hg38", "validation_cells": 10, "output_object": "RNA+circ MuData; circ AnnData",
         "circRNA_candidates": "not regenerated", "status": "real 10-cell workflow; committed report",
         "role": "protocol validation only", "limitation": "RT contract synthetic-tested; historical full object excluded", "source": STATUS},
    ])
    inventory.to_csv(output / "table_s1_dataset_workflow_inventory.tsv", sep="\t", index=False)

    rows = [
        ("Frozen full suite", "source tests", "345 passed; 8 skipped; 5 warnings", "committed frozen validation", BASELINE, "not rerun as a raw-data pipeline"),
        ("MuData transition warnings", "synchronization regression", "0 behavior-change FutureWarnings", "committed frozen validation", "docs/mudata_compatibility.md", "broad AnnData 0.13 migration deferred"),
        ("Distribution build", "wheel and sdist", "passed", "committed frozen validation", BASELINE, "version 0.10.0"),
        ("Clean wheel install", "declared dependencies; pip check; CLI", "passed", "committed frozen validation", BASELINE, "external detector executables separate"),
        ("Packaged detector resources", "installed resource lookup", "passed", "committed frozen validation", BASELINE, "not an accuracy benchmark"),
        ("H5MU round trip", "synthetic RNA+circ", "passed; warning-free", "committed frozen validation", BASELINE, "does not guarantee legacy global-index synchronization"),
        ("Smart-seq3 object", "checksum-matched processed object", f"{smart['cell_count']} cells; {smart['circRNA_candidate_count']} candidates; {smart['host_gene_annotated_count']} annotated ({100 * smart['host_gene_annotation_fraction']:.1f}%)", "regenerated", "results_v2/provenance.json", "host annotation is not independent circRNA validation"),
        ("IMR90 object", "checksum-matched processed object", f"{imr['trimodal_cell_overlap']} shared cells; {imr['circRNA_shape'][1]} candidates; {imr['host_gene_annotated_count']} annotated ({100 * imr['host_gene_annotation_fraction']:.1f}%)", "regenerated", "results_v2/provenance.json", "processed CNV interoperability only"),
        ("HAP1 RNA+circ", "real paired-end scRR", "10 cells; candidates/annotation not regenerated", "committed workflow validation", STATUS, "no historical full-object or real RT result used"),
        ("Optional RT", "synthetic import/merge fixtures", "contract tests passed", "synthetic only", WORKFLOWS, "real processed-file validation pending reconciliation"),
    ]
    pd.DataFrame(rows, columns=["check", "dataset_or_fixture", "result", "evidence_class", "source", "limitation"]).to_csv(
        output / "table_s2_software_reproducibility.tsv", sep="\t", index=False)

    versions = [("software", "circyto", "0.10.0 (source tree)", "frozen production baseline", BASELINE),
                ("software", "baseline commit", provenance["repository"]["frozen_baseline"], "frozen", BASELINE),
                ("regeneration", "repository commit", provenance["repository"]["commit"], "script digest disambiguates uncommitted work", "results_v2/provenance.json"),
                ("regeneration", "script SHA-256", provenance["script_sha256"], "executed source", "scripts/manuscript/regenerate_application_note_results.py"),
                ("environment", "Python", provenance["python"], "regeneration environment", "results_v2/provenance.json"),
                ("command", "results and Figure 1B/C", "python scripts/manuscript/regenerate_application_note_results.py --smartseq3 <SMART_H5MU> --imr90 <IMR90_H5MU> --output manuscript/results_v2_reproduced", "offline after objects are obtained", "manuscript/reproduce.md"),
                ("command", "schematics and supplement", "python scripts/manuscript/build_application_note_supplement.py --results manuscript/results_v2_reproduced --output manuscript/supplement_v2_reproduced", "offline; new directory required", "manuscript/reproduce.md"),
                ("source", "object archive commit", "c99cddaae2314d3af7d5232150c5b86bb6b52e25", "historical Git source; binaries excluded from current tree", EVIDENCE),
                ("external tools", "STAR / BWA / CIRI3", "not executed for processed-object regeneration", "historical executable versions not independently re-established in this pass", WORKFLOWS)]
    versions.extend(("environment", key, value, "regeneration environment", "results_v2/provenance.json") for key, value in provenance["package_versions"].items() if key != "circyto")
    versions.extend(("input SHA-256", item["dataset"], item["observed_sha256"], "verified before and after reading", item["path"]) for item in provenance["input_checksums"])
    versions.extend(("embedding", key, json.dumps(value), "RNA-only; historical method; seed 17", "results_v2/provenance.json") for key, value in provenance["umap_method"].items())
    pd.DataFrame(versions, columns=["category", "quantity", "value", "status", "source"]).to_csv(
        output / "table_s3_versions_commands_provenance.tsv", sep="\t", index=False)

    # S1 is a compact visual comparison; Table S1 retains accessions and sources.
    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    ax.axis("off")
    ax.set_title("S1  Protocol and workflow validation", loc="left", weight="bold", pad=18)
    table = ax.table(cellText=[
        ["Smart-seq3", "Paired-end RNA\n+ index reads", "Demultiplex → STAR\n→ BWA rescue → CIRI3", str(smart["cell_count"]), "RNA + circ\nMuData / AnnData", "Real data\nObject regenerated"],
        ["IMR90 scRR", "Single-end RNA\n+ processed CNV", "BWA-MEM → CIRI3\n→ cell remap + CNV", str(imr["rna_shape"][0]), "RNA + circ + CNV\nMuData", "Real data\nObject regenerated"],
        ["HAP1 scRR", "Paired-end RNA", "STAR → BWA rescue\n→ CIRI3", "10", "RNA + circ\nMuData / AnnData", "Real batch10\nCommitted report"],
    ], colLabels=["Protocol", "Read layout", "Alignment / detector route", "Cells", "Output object", "Status"],
       colWidths=[.12, .17, .25, .06, .19, .21], cellLoc="left", bbox=[0, .15, 1, .8])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#ccd5df")
        cell.set_facecolor("#e7eff6" if row == 0 else ("#f7f9fb" if row % 2 else "white"))
        if row == 0:
            cell.set_text_props(weight="bold")
    ax.text(0, .035, "Workflow execution and object interoperability; no independent circRNA accuracy claim.\nHAP1 uses the validated 10-cell RNA/circ route; real RT integration remains unresolved.", transform=ax.transAxes, fontsize=9)
    fig.subplots_adjust(left=.03, right=.97, top=.87, bottom=.07)
    save(fig, "figure_s1_protocol_workflows")

    def box(ax, x, y, w, h, title, detail, color="#e7eff6", dashed=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012", facecolor=color,
                                   edgecolor="#576b80", linestyle="--" if dashed else "-", linewidth=1))
        ax.text(x+w/2, y+h*.72, title, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x+w/2, y+h*.35, detail, ha="center", va="center", fontsize=9)

    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    ax.text(.02, .96, "S2  Modality-specific matrices, explicit cell identities", fontsize=14, weight="bold")
    box(ax, .03, .3, .28, .49, "RNA + circ", f"Smart-seq3: real processed object\nRNA: {smart['cell_count']} × {smart['rna_feature_count']:,}\ncirc: {smart['cell_count']} × {smart['circRNA_candidate_count']:,}\n{smart['cell_count']} shared cell IDs")
    box(ax, .36, .3, .28, .49, "RNA + circ + CNV", f"IMR90: real processed object\nRNA: {imr['rna_shape'][0]} × {imr['rna_shape'][1]:,}\ncirc: {imr['circRNA_shape'][0]} × {imr['circRNA_shape'][1]:,}\nCNV: {imr['cnv_shape'][0]} × {imr['cnv_shape'][1]:,}\n{imr['trimodal_cell_overlap']} shared cell IDs")
    box(ax, .69, .3, .28, .49, "Optional RNA + circ + RT", "Implemented integration contract\nSynthetic import / merge tests\nReal-file validation unresolved\nNo HAP1 RT result shown", color="#f5f2e8", dashed=True)
    ax.text(.5, .18, "Each modality: AnnData X[cell, feature] + obs + var  |  MuData: mod[rna, circ, …]", ha="center", fontsize=10)
    ax.text(.5, .07, "Cell mappings preserve identity; feature axes remain modality-specific. Processed CNV/RT are imported summaries.", ha="center", fontsize=9)
    fig.subplots_adjust(left=.01, right=.99, top=.97, bottom=.04)
    save(fig, "figure_s2_multimodal_schema")

    # No existing architecture asset exists in this branch; render its documented scaffold.
    fig, ax = plt.subplots(figsize=(10.5, 3.5))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    ax.text(.02, .94, "A  circyto connects established detectors to scverse", fontsize=14, weight="bold")
    steps = [(.02, "Full-length scRNA-seq", "Protocol-aware input\nCell identities"),
             (.27, "Established detectors", "Alignment + detector\norchestration"),
             (.52, "Matrix + annotation", "Cell × circRNA; QC\nHost genes; provenance"),
             (.77, "AnnData / MuData", "RNA + circ\nscverse analysis")]
    for x, title, detail in steps:
        box(ax, x, .39, .20, .37, title, detail)
    for x in (.23, .48, .73):
        ax.annotate("", xy=(x+.025, .58), xytext=(x-.005, .58), arrowprops={"arrowstyle": "->", "color": "#576b80"})
    ax.text(.5, .21, "Optional: processed CNV (validated IMR90) · RT (synthetic contract) · candidate variants (exploratory)", ha="center", fontsize=9, color="#56616b")
    ax.text(.5, .07, "circyto is not a new back-splice detection algorithm.", ha="center", fontsize=10)
    fig.subplots_adjust(left=.01, right=.99, top=.97, bottom=.04)
    save(fig, "figure1a_architecture")

    sources = [EVIDENCE, WORKFLOWS, BASELINE, STATUS, "docs/modality_schema.md", "docs/manuscript_figure_skeleton.md"]
    write_json(output / "provenance.json", {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/manuscript/build_application_note_supplement.py",
        "script_sha256": sha256_file(Path(__file__)),
        "repository": repository_provenance(), "package_versions": package_versions(),
        "results_provenance_sha256": sha256_file(results / "provenance.json"),
        "source_sha256": {path: sha256_file(REPO_ROOT / path) for path in sources},
        "generated_files": {path.name: sha256_file(path) for path in sorted(output.iterdir()) if path.is_file()},
        "randomness": "No random layout or sampling; fixed figure geometry; svg.hashsalt fixed",
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
