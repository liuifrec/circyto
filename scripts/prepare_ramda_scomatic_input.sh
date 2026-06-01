#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
usage: prepare_ramda_scomatic_input.sh --alignment-manifest alignment_manifest.tsv --outdir OUTDIR [options]

Prepare one-cell-per-SRR RamDA/scRR alignments for SComatic's pooled BAM input.
This script does not run SComatic.

Required:
  --alignment-manifest PATH   circyto alignment_manifest.tsv with cell_id and bam or sam columns
  --outdir PATH               Output directory

Options:
  --sample-id ID              Sample prefix for merged BAM [default: ramda_scomatic]
  --threads N                 Threads for samtools sort/merge [default: 1]
  --cell-type-column COLUMN   Manifest column to use for SComatic Cell_type
  --default-cell-type LABEL   Fallback Cell_type for all rows. If omitted, Cell_type defaults to cell_id.
  -h, --help                  Show this help

Outputs:
  OUTDIR/per_cell/*.scomatic.sorted.bam
  OUTDIR/merged/<sample-id>.scomatic.bam
  OUTDIR/cell_annotations.tsv
  OUTDIR/adapter_summary.json
EOF
}

ALIGNMENT_MANIFEST=""
OUTDIR=""
SAMPLE_ID="ramda_scomatic"
THREADS="1"
CELL_TYPE_COLUMN=""
DEFAULT_CELL_TYPE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --alignment-manifest)
      ALIGNMENT_MANIFEST="${2:?--alignment-manifest requires a path}"
      shift 2
      ;;
    --outdir)
      OUTDIR="${2:?--outdir requires a path}"
      shift 2
      ;;
    --sample-id)
      SAMPLE_ID="${2:?--sample-id requires a value}"
      shift 2
      ;;
    --threads)
      THREADS="${2:?--threads requires a value}"
      shift 2
      ;;
    --cell-type-column)
      CELL_TYPE_COLUMN="${2:?--cell-type-column requires a value}"
      shift 2
      ;;
    --default-cell-type)
      DEFAULT_CELL_TYPE="${2:?--default-cell-type requires a value}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$ALIGNMENT_MANIFEST" || -z "$OUTDIR" ]]; then
  echo "[ERROR] --alignment-manifest and --outdir are required" >&2
  usage >&2
  exit 2
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "[ERROR] --threads must be a positive integer" >&2
  exit 2
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools not found on PATH" >&2
  exit 127
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "[ERROR] python3 not found on PATH" >&2
  exit 127
fi

mkdir -p "$OUTDIR"
ROWS_TSV="$OUTDIR/adapter_rows.tsv"
ANNOTATION_TSV="$OUTDIR/cell_annotations.tsv"
SUMMARY_JSON="$OUTDIR/adapter_summary.json"
PER_CELL_DIR="$OUTDIR/per_cell"
TAG_DIR="$OUTDIR/tag_reports"
MERGED_DIR="$OUTDIR/merged"
mkdir -p "$PER_CELL_DIR" "$TAG_DIR" "$MERGED_DIR"

python3 - "$ALIGNMENT_MANIFEST" "$OUTDIR" "$SAMPLE_ID" "$CELL_TYPE_COLUMN" "$DEFAULT_CELL_TYPE" <<'PY'
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
outdir = Path(sys.argv[2])
sample_id = sys.argv[3]
cell_type_column = sys.argv[4]
default_cell_type = sys.argv[5]

if not manifest.exists():
    raise SystemExit(f"[ERROR] alignment manifest not found: {manifest}")

rows_tsv = outdir / "adapter_rows.tsv"
annotation_tsv = outdir / "cell_annotations.tsv"
per_cell_dir = outdir / "per_cell"
tag_dir = outdir / "tag_reports"


def resolve_path(value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest.parent / path
    return path.resolve(strict=False)


def safe_file_id(value: str, fallback: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return safe or fallback


def safe_cell_type(value: str, fallback: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", value).strip("_-")
    safe = re.sub(r"_+", "_", safe)
    return safe or fallback


with manifest.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames or []
    if "cell_id" not in fieldnames:
        raise SystemExit(f"[ERROR] alignment manifest missing required column: cell_id ({manifest})")
    if "bam" not in fieldnames and "sam" not in fieldnames:
        raise SystemExit(f"[ERROR] alignment manifest must contain bam or sam column: {manifest}")
    if cell_type_column and cell_type_column not in fieldnames:
        raise SystemExit(f"[ERROR] --cell-type-column {cell_type_column!r} not found in {manifest}")

    seen_cell_ids: set[str] = set()
    seen_file_ids: dict[str, int] = {}
    parsed: list[dict[str, str]] = []
    annotations: list[dict[str, str]] = []

    for line_number, row in enumerate(reader, start=2):
        cell_id = (row.get("cell_id") or "").strip()
        if not cell_id:
            raise SystemExit(f"[ERROR] empty cell_id at {manifest}:{line_number}")
        if re.search(r"\s", cell_id):
            raise SystemExit(f"[ERROR] cell_id contains whitespace, which is unsafe for SComatic CB tags: {cell_id!r}")
        if cell_id in seen_cell_ids:
            raise SystemExit(f"[ERROR] duplicate cell_id in alignment manifest: {cell_id}")
        seen_cell_ids.add(cell_id)

        bam = (row.get("bam") or "").strip()
        sam = (row.get("sam") or "").strip()
        if bool(bam) == bool(sam):
            raise SystemExit(f"[ERROR] row for cell_id={cell_id} must provide exactly one of bam or sam")
        alignment_path = resolve_path(bam or sam)
        if not alignment_path.exists():
            raise SystemExit(f"[ERROR] alignment file not found for cell_id={cell_id}: {alignment_path}")

        raw_cell_type = ""
        if cell_type_column:
            raw_cell_type = (row.get(cell_type_column) or "").strip()
        if not raw_cell_type:
            raw_cell_type = default_cell_type.strip() if default_cell_type.strip() else cell_id
        cell_type = safe_cell_type(raw_cell_type, fallback=safe_cell_type(cell_id, fallback=f"cell_{line_number}"))

        base_file_id = safe_file_id(cell_id, fallback=f"cell_{line_number}")
        count = seen_file_ids.get(base_file_id, 0) + 1
        seen_file_ids[base_file_id] = count
        file_id = base_file_id if count == 1 else f"{base_file_id}_{count}"

        sorted_bam = per_cell_dir / f"{file_id}.scomatic.sorted.bam"
        tag_report = tag_dir / f"{file_id}.tags.tsv"
        parsed.append(
            {
                "file_id": file_id,
                "cell_id": cell_id,
                "alignment_path": str(alignment_path),
                "alignment_kind": "bam" if bam else "sam",
                "cell_type": cell_type,
                "sorted_bam": str(sorted_bam),
                "tag_report": str(tag_report),
            }
        )
        annotations.append({"Index": cell_id, "Cell_type": cell_type})

if not parsed:
    raise SystemExit(f"[ERROR] alignment manifest contains 0 rows: {manifest}")

with rows_tsv.open("w", encoding="utf-8", newline="") as handle:
    fieldnames = ["file_id", "cell_id", "alignment_path", "alignment_kind", "cell_type", "sorted_bam", "tag_report"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(parsed)

with annotation_tsv.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["Index", "Cell_type"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(annotations)
PY

SORTED_BAMS=()
while IFS=$'\t' read -r file_id cell_id alignment_path alignment_kind cell_type sorted_bam tag_report; do
  [[ "$file_id" == "file_id" ]] && continue
  echo "[INFO] preparing cell_id=${cell_id} input=${alignment_path}" >&2
  mkdir -p "$(dirname "$sorted_bam")" "$(dirname "$tag_report")"
  samtools view -h "$alignment_path" \
    | awk -v cb="$cell_id" -v report="$tag_report" '
      BEGIN {
        OFS = "\t";
        total = 0;
        missing_cb = 0;
        replaced_cb = 0;
        missing_nM = 0;
        missing_NH = 0;
        expected_cb = "CB:Z:" cb;
      }
      /^@/ { print; next }
      NF == 0 { print; next }
      {
        total++;
        has_cb = 0;
        has_nM = 0;
        has_NH = 0;
        for (i = 12; i <= NF; i++) {
          if ($i ~ /^CB:Z:/) {
            has_cb = 1;
            if ($i != expected_cb) {
              replaced_cb++;
              $i = expected_cb;
            }
          }
          if ($i ~ /^nM:/) { has_nM = 1; }
          if ($i ~ /^NH:/) { has_NH = 1; }
        }
        if (!has_cb) {
          missing_cb++;
          $(NF + 1) = expected_cb;
        }
        if (!has_nM) { missing_nM++; }
        if (!has_NH) { missing_NH++; }
        print;
      }
      END {
        print "total_reads", "missing_cb_before", "replaced_cb", "missing_nM", "missing_NH" > report;
        print total, missing_cb, replaced_cb, missing_nM, missing_NH >> report;
      }
    ' \
    | samtools view -b - \
    | samtools sort -@ "$THREADS" -o "$sorted_bam" -
  samtools index "$sorted_bam"
  SORTED_BAMS+=("$sorted_bam")
done < "$ROWS_TSV"

MERGED_BAM="$MERGED_DIR/${SAMPLE_ID}.scomatic.bam"
if [[ "${#SORTED_BAMS[@]}" -eq 1 ]]; then
  cp "${SORTED_BAMS[0]}" "$MERGED_BAM"
else
  samtools merge -@ "$THREADS" -f "$MERGED_BAM" "${SORTED_BAMS[@]}"
fi
samtools index "$MERGED_BAM"

python3 - "$ALIGNMENT_MANIFEST" "$OUTDIR" "$SAMPLE_ID" "$ROWS_TSV" "$ANNOTATION_TSV" "$MERGED_BAM" "$SUMMARY_JSON" <<'PY'
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

manifest = Path(sys.argv[1]).resolve()
outdir = Path(sys.argv[2]).resolve()
sample_id = sys.argv[3]
rows_tsv = Path(sys.argv[4])
annotation_tsv = Path(sys.argv[5]).resolve()
merged_bam = Path(sys.argv[6]).resolve()
summary_json = Path(sys.argv[7])

rows: list[dict[str, str]] = []
with rows_tsv.open("r", encoding="utf-8", newline="") as handle:
    rows.extend(csv.DictReader(handle, delimiter="\t"))

total_reads = 0
missing_cb = 0
replaced_cb = 0
missing_nM = 0
missing_NH = 0
for row in rows:
    report = Path(row["tag_report"])
    with report.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        record = next(reader)
    total_reads += int(record["total_reads"])
    missing_cb += int(record["missing_cb_before"])
    replaced_cb += int(record["replaced_cb"])
    missing_nM += int(record["missing_nM"])
    missing_NH += int(record["missing_NH"])

payload = {
    "command_name": "scripts/prepare_ramda_scomatic_input.sh",
    "description": "Prepared one-cell-per-SRR RamDA/scRR alignments for SComatic pooled BAM input.",
    "scomatic_result_terminology": "RNA-derived candidate variant signals",
    "sample_id": sample_id,
    "alignment_manifest": str(manifest),
    "outdir": str(outdir),
    "input_cells": len(rows),
    "merged_bam": str(merged_bam),
    "merged_bam_index": str(Path(str(merged_bam) + ".bai")),
    "cell_annotation_tsv": str(annotation_tsv),
    "per_cell_bams": [
        {
            "cell_id": row["cell_id"],
            "cell_type": row["cell_type"],
            "input_alignment": row["alignment_path"],
            "sorted_bam": row["sorted_bam"],
            "tag_report": row["tag_report"],
        }
        for row in rows
    ],
    "tag_summary": {
        "total_reads": total_reads,
        "missing_cb_before_injection": missing_cb,
        "replaced_cb_not_matching_cell_id": replaced_cb,
        "missing_nM": missing_nM,
        "missing_NH": missing_NH,
    },
    "scomatic_tag_contract": {
        "CB": "required by SComatic SplitBamCellTypes.py; injected or normalized to cell_id",
        "nM": "preserved and counted; required when running SComatic with --max_nM",
        "NH": "preserved and counted; required when running SComatic with --max_NH",
    },
    "runs_scomatic": False,
}
summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(str(summary_json))
PY
