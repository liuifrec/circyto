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
import sys

runinfo_csv = pathlib.Path(sys.argv[1])
manifest_out = pathlib.Path(sys.argv[2])
inventory_out = pathlib.Path(sys.argv[3])

if not runinfo_csv.exists():
    raise SystemExit(f"missing RunInfo CSV: {runinfo_csv}")

TARGET = {
    "SRR30918126": {
        "gsm": "GSM8558852",
        "sample_id": "SRR30918126_IMR90_scRR",
        "condition": "IMR90_aphidicolin_treated",
        "cell_cycle_phase": "G1",
    },
    "SRR30918117": {
        "gsm": "GSM8558853",
        "sample_id": "SRR30918117_IMR90_scRR",
        "condition": "IMR90_aphidicolin_treated",
        "cell_cycle_phase": "G1",
    },
}

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
        row_get(row, "scientific_name", "ScientificName"),
        row_get(row, "source_name", "Source Name"),
    ]
    return " | ".join(v for v in fields if v)

rows = []
with runinfo_csv.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        srr = row_get(row, "Run", "run", "SRR", "srr")
        if not srr or srr not in TARGET:
            continue
        meta = combined_metadata(row).lower()
        if "dna" in meta or "exome" in meta:
            continue
        record = TARGET[srr].copy()
        record["srr"] = srr
        record["gsm"] = row_get(row, "GEO_Accession", "geo_accession", "GEO Accession") or record["gsm"]
        record["read_layout"] = infer_layout(row) or "single"
        record["library_strategy"] = row_get(row, "LibraryStrategy", "library_strategy")
        record["library_source"] = row_get(row, "LibrarySource", "library_source")
        record["organism"] = row_get(row, "ScientificName", "scientific_name")
        record["title"] = row_get(row, "title", "Title")
        record["sample_name"] = row_get(row, "SampleName", "sample_name", "Sample Name")
        record["spots"] = row_get(row, "spots", "Spots")
        record["bases"] = row_get(row, "bases", "Bases")
        rows.append(record)

rows.sort(key=lambda r: r["srr"])

if not rows:
    raise SystemExit("no validated IMR90 RNA-side SRRs were found in the supplied RunInfo CSV")

with manifest_out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
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
    ])
    for row in rows:
        writer.writerow([
            row["sample_id"],
            f"raw/{row['srr']}.fastq.gz",
            "",
            "ramda",
            "unstranded",
            row["read_layout"],
            row["srr"],
            row["gsm"],
            row["condition"],
            row["cell_cycle_phase"],
        ])

with inventory_out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "sample_id",
        "srr",
        "gsm",
        "condition",
        "cell_cycle_phase",
        "read_layout",
        "library_strategy",
        "library_source",
        "organism",
        "sample_name",
        "title",
        "spots",
        "bases",
    ])
    for row in rows:
        writer.writerow([
            row["sample_id"],
            row["srr"],
            row["gsm"],
            row["condition"],
            row["cell_cycle_phase"],
            row["read_layout"],
            row["library_strategy"],
            row["library_source"],
            row["organism"],
            row["sample_name"],
            row["title"],
            row["spots"],
            row["bases"],
        ])

print(f"[INFO] wrote {manifest_out}")
print(f"[INFO] wrote {inventory_out}")
print(f"[INFO] retained {len(rows)} validated IMR90 RNA-side SRRs")
PY
