#!/usr/bin/env python3
"""Create deterministic, explicitly synthetic report inputs. No detector runs.

Usage: python examples/qc_report_demo.py --outdir work/qc_report_demo
Then:  circyto report --workdir work/qc_report_demo
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def create_demo(outdir: Path) -> None:
    if outdir.exists():
        raise ValueError("Use a new --outdir; the demo never overwrites an existing directory.")
    values = [
        [8, 3, 0, 0, 0, 0], [6, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0],
        [10, 2, 3, 1, 0, 0], [3, 0, 0, 0, 5, 0], [0, 5, 0, 0, 0, 0],
        [11, 1, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0], [1, 2, 0, 0, 0, 0],
        [2, 0, 1, 0, 0, 0], [0, 0, 0, 0, 2, 0], [5, 0, 0, 0, 0, 3],
    ]
    cells = [f"demo_{row}{col:02d}" for row in "AB" for col in range(1, 7)]
    candidates = [f"chr1:{1000 + i * 1000}|{1400 + i * 1000}" for i in range(6)]
    genes = ["DEMO_GENE_A", "DEMO_GENE_B", "", "DEMO_GENE_D", "DEMO_GENE_E", ""]
    for folder in ("qc", "matrix", "inputs"):
        (outdir / folder).mkdir(parents=True)
    with (outdir / "qc/cell_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cell_id", "circRNA_count", "total_circRNA_support", "detector_status", "alignment_status"])
        for cell, row in zip(cells, values):
            writer.writerow([cell, sum(value > 0 for value in row), sum(row), "success" if any(row) else "empty", "aligned"])
    with (outdir / "qc/circ_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["circ_id", "n_cells_detected", "host_gene", "total_support"])
        for index, candidate in enumerate(candidates):
            column = [row[index] for row in values]
            writer.writerow([candidate, sum(value > 0 for value in column), genes[index], sum(column)])
    entries = [(candidate + 1, cell + 1, value) for cell, row in enumerate(values) for candidate, value in enumerate(row) if value]
    (outdir / "matrix/circ_counts.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n% Synthetic documentation fixture; not biological data\n"
        + f"6 12 {len(entries)}\n" + "".join(f"{candidate} {cell} {value}\n" for candidate, cell, value in entries), encoding="utf-8")
    (outdir / "matrix/circ_index.txt").write_text("\n".join(candidates) + "\n", encoding="utf-8")
    (outdir / "matrix/cell_index.txt").write_text("\n".join(cells) + "\n", encoding="utf-8")
    summary = {
        "workflow_type": "full-length-circrna", "protocol": "ramda", "read_layout": "single-end",
        "detector_backend": "ciri3", "circyto_version": "synthetic fixture (no analysis run)",
        "started_at": "2026-01-01T00:00:00Z", "completed_at": "2026-01-01T00:01:00Z", "dry_run": False,
        "workflow_uuid": "synthetic-report-demonstration", "planned_cells": len(cells),
        "matrix": {"n_cells": 12, "n_circRNAs": 6},
        "alignment_status_counts": {"aligned": 12}, "detector_status_counts": {"success": 10, "empty": 2},
        "stage_graph": [{"stage": stage, "status": "completed"} for stage in ("alignment", "detector", "matrix", "summary_qc")],
        "warnings": ["SYNTHETIC DEMONSTRATION: counts, identifiers, annotation labels, statuses, and timestamps are illustrative. No biological analysis or CIRI3 execution was performed."],
    }
    (outdir / "workflow_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # Separate toy inputs make the documented workflow --dry-run executable.
    (outdir / "inputs/reads.fastq").write_text("@demo_read\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")
    (outdir / "inputs/genome.fa").write_text(">chr1\nACGTACGTACGTACGT\n", encoding="utf-8")
    (outdir / "inputs/genes.gtf").write_text('chr1\tdemo\texon\t1\t16\t.\t+\t.\tgene_id "DEMO_GENE"; transcript_id "DEMO_TX";\n', encoding="utf-8")
    (outdir / "inputs/manifest.tsv").write_text(
        "cell_id\tplatform\tread1\tread2\tbam\tlibrary_id\tn_input_reads\tprotocol\tstrandedness\tread_layout\n"
        f"demo_cell\tplate\t{(outdir / 'inputs/reads.fastq').resolve()}\t\t\tdemo\t1\tramda\tunstranded\tsingle\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    create_demo(args.outdir)
    print(f"Synthetic report inputs: {args.outdir}")
    print(f"circyto report --workdir {args.outdir}")


if __name__ == "__main__":
    main()
