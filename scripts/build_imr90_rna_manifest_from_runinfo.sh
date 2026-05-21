#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 <SraRunTable.csv> [output_dir]" >&2
  exit 1
fi

RUNINFO_CSV="$1"
OUTDIR="${2:-.}"

mkdir -p "${OUTDIR}"

MANIFEST_OUT="${OUTDIR}/manifest_imr90_rna_all.tsv"
INVENTORY_OUT="${OUTDIR}/imr90_rna_run_inventory.tsv"

python - "$RUNINFO_CSV" "$MANIFEST_OUT" "$INVENTORY_OUT" <<'PY'
import csv
import pathlib
import re
import sys
from collections import Counter

runinfo_csv = pathlib.Path(sys.argv[1])
manifest_out = pathlib.Path(sys.argv[2])
inventory_out = pathlib.Path(sys.argv[3])

if not runinfo_csv.exists():
    raise SystemExit(f"missing RunInfo CSV: {runinfo_csv}")


def normalize_key(key: str) -> str:
    return "".join(ch.lower() for ch in key if ch.isalnum())


def row_get(row, *names):
    normalized = {normalize_key(k): v for k, v in row.items()}
    for name in names:
        value = normalized.get(normalize_key(name), "")
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def infer_layout(row):
    layout = row_get(row, "LibraryLayout", "library_layout", "layout").upper()
    if "PAIRED" in layout:
        return "paired"
    if "SINGLE" in layout:
        return "single"
    return ""


def combined_metadata(row):
    fields = [
        row_get(row, "SampleName", "sample_name", "Sample Name"),
        row_get(row, "title", "Title"),
        row_get(row, "Experiment", "experiment"),
        row_get(row, "LibraryName", "library_name", "Library Name"),
        row_get(row, "GEO_Accession", "geo_accession", "GEO Accession"),
        row_get(row, "BioSample", "biosample"),
        row_get(row, "ScientificName", "scientific_name"),
        row_get(row, "source_name", "Source Name"),
    ]
    return " | ".join(v for v in fields if v)


def infer_treatment(text: str) -> str:
    lower = text.lower()
    if "aphidicolin" in lower:
        return "aphidicolin"
    return ""


def infer_cell_cycle_phase(text: str) -> str:
    lower = text.lower()
    if re.search(r"\bg1\b", lower):
        return "G1"
    if "mid-s" in lower or "mids" in lower or re.search(r"\bmid s\b", lower):
        return "mid-S"
    if re.search(r"\bs\b", lower):
        return "S"
    return ""


def infer_condition(text: str, treatment: str, phase: str) -> str:
    tokens = ["IMR90"]
    if treatment:
        tokens.append(treatment)
    if phase:
        tokens.append(phase)
    return "_".join(tokens)


def is_rna_side(row):
    assay_type = row_get(row, "Assay Type", "assay_type", "LibraryStrategy", "library_strategy").upper()
    library_source = row_get(row, "LibrarySource", "library_source").upper()
    library_selection = row_get(row, "LibrarySelection", "library_selection").lower()
    metadata = combined_metadata(row).lower()
    if assay_type != "RNA-SEQ":
        return False
    if "TRANSCRIPTOMIC" not in library_source:
        return False
    if library_selection != "cdna":
        return False
    if "dna" in metadata or "exome" in metadata:
        return False
    return True


rows = []
stats = Counter()
layout_counts = Counter()

with runinfo_csv.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        stats["total_rows"] += 1
        assay_type = row_get(row, "Assay Type", "assay_type", "LibraryStrategy", "library_strategy").upper()
        library_source = row_get(row, "LibrarySource", "library_source").upper()
        metadata = combined_metadata(row).lower()
        if assay_type == "OTHER" or library_source == "GENOMIC" or "dna" in metadata or "exome" in metadata:
            stats["excluded_dna_or_other"] += 1
        if not is_rna_side(row):
            continue
        srr = row_get(row, "Run", "run", "SRR", "srr")
        if not srr:
            continue
        read_layout = infer_layout(row) or "single"
        if read_layout != "single":
            continue
        meta = combined_metadata(row)
        treatment = row_get(row, "treatment", "Treatment") or infer_treatment(meta)
        cell_cycle_phase = row_get(row, "cell_cycle_phase", "CellCyclePhase", "cell cycle phase") or infer_cell_cycle_phase(meta)
        gsm = row_get(row, "GEO_Accession", "geo_accession", "GEO Accession")
        library_name = row_get(row, "LibraryName", "library_name", "Library Name")
        sample_id = library_name or gsm or f"{srr}_IMR90_scRR"
        rows.append(
            {
                "sample_id": sample_id,
                "srr": srr,
                "gsm": gsm,
                "library_name": library_name,
                "treatment": treatment,
                "condition": infer_condition(meta, treatment, cell_cycle_phase),
                "cell_cycle_phase": cell_cycle_phase,
                "biosample": row_get(row, "BioSample", "biosample"),
                "experiment": row_get(row, "Experiment", "experiment"),
                "read_layout": read_layout,
                "library_strategy": row_get(row, "LibraryStrategy", "library_strategy"),
                "library_source": row_get(row, "LibrarySource", "library_source"),
                "library_selection": row_get(row, "LibrarySelection", "library_selection"),
                "organism": row_get(row, "ScientificName", "scientific_name"),
                "sample_name": row_get(row, "SampleName", "sample_name", "Sample Name"),
                "title": row_get(row, "title", "Title"),
                "spots": row_get(row, "spots", "Spots"),
                "bases": row_get(row, "bases", "Bases"),
            }
        )
        stats["retained_rna_rows"] += 1
        layout_counts[read_layout] += 1

rows.sort(key=lambda r: r["srr"])

if not rows:
    raise SystemExit("no IMR90 RNA-side rows were found in the supplied RunInfo CSV")

manifest_fieldnames = [
    "sample_id",
    "fastq_1",
    "fastq_2",
    "protocol",
    "strandedness",
    "read_layout",
    "srr",
    "gsm",
    "condition",
    "cell_cycle_phase",
]

with manifest_out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=manifest_fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                "sample_id": row["sample_id"],
                "fastq_1": f"raw/{row['srr']}.fastq.gz",
                "fastq_2": "",
                "protocol": "ramda",
                "strandedness": "unstranded",
                "read_layout": "single",
                "srr": row["srr"],
                "gsm": row["gsm"],
                "condition": row["condition"],
                "cell_cycle_phase": row["cell_cycle_phase"],
            }
        )

with inventory_out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample_id",
            "srr",
            "gsm",
            "library_name",
            "treatment",
            "condition",
            "cell_cycle_phase",
            "biosample",
            "experiment",
            "read_layout",
            "library_strategy",
            "library_source",
            "library_selection",
            "organism",
            "sample_name",
            "title",
            "spots",
            "bases",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

with manifest_out.open(newline="", encoding="utf-8") as handle:
    manifest_reader = csv.DictReader(handle, delimiter="\t")
    for row in manifest_reader:
        sample_id = (row.get("sample_id") or "").strip()
        read_layout = (row.get("read_layout") or "").strip()
        if read_layout not in {"single", "single-end"}:
            raise SystemExit(
                f"manifest validation failed for {sample_id or '<unknown>'}: invalid read_layout {read_layout!r}"
            )
        if not (row.get("fastq_1") or "").strip():
            raise SystemExit(
                f"manifest validation failed for {sample_id or '<unknown>'}: fastq_1 is empty"
            )

print(f"[INFO] wrote {manifest_out}")
print(f"[INFO] wrote {inventory_out}")
print(f"[INFO] total RunInfo rows: {stats['total_rows']}")
print(f"[INFO] RNA-side rows retained: {stats['retained_rna_rows']}")
print(f"[INFO] DNA/genomic/OTHER rows excluded: {stats['excluded_dna_or_other']}")
print(f"[INFO] layout counts: {dict(layout_counts)}")
PY
