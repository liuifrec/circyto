#!/usr/bin/env bash
set -euo pipefail

WORKDIR=""
REFERENCE=""
GTF=""
OUTDIR=""
SYNTHETIC=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)
      WORKDIR="$2"
      shift 2
      ;;
    --reference)
      REFERENCE="$2"
      shift 2
      ;;
    --gtf)
      GTF="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --synthetic)
      SYNTHETIC=1
      shift
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${WORKDIR}" || -z "${REFERENCE}" || -z "${GTF}" || -z "${OUTDIR}" ]]; then
  echo "usage: $0 --workdir WORKDIR --reference ref/chr21.fa --gtf ref/chr21.gtf --outdir OUTDIR [--synthetic]" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

if [[ "${SYNTHETIC}" -eq 1 ]]; then
  python - "$WORKDIR" "$REFERENCE" "$GTF" "$OUTDIR" <<'PY'
from pathlib import Path
import json
import sys

from circyto.pipeline.scomatic_chr21_poc import write_synthetic_scomatic_poc

summary = write_synthetic_scomatic_poc(
    workdir=Path(sys.argv[1]),
    reference=Path(sys.argv[2]),
    gtf=Path(sys.argv[3]),
    outdir=Path(sys.argv[4]),
)
print(json.dumps(summary, indent=2, sort_keys=True))
PY
  exit 0
fi

ALIGNMENT_PATH="$(find "${WORKDIR}" -type f \( -name '*.bam' -o -name '*.sam' \) | head -n 1 || true)"
if [[ -z "${ALIGNMENT_PATH}" ]]; then
  echo "[ERROR] no BAM/SAM discovered under ${WORKDIR}" >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools is required for non-synthetic chr21 SComatic POC when BAM/SAM conversion or indexing is needed" >&2
  exit 1
fi

if ! python - <<'PY' >/dev/null 2>&1
import importlib.util
import sys
sys.exit(0 if importlib.util.find_spec("SComatic") else 1)
PY
then
  echo "[ERROR] SComatic is not installed. Install it in the local env or rerun with --synthetic." >&2
  echo "[ERROR] Suggested install: pip install SComatic" >&2
  exit 1
fi

echo "[ERROR] Real local chr21 SComatic execution is scaffolded but not enabled in this script revision." >&2
echo "[ERROR] Use --synthetic to validate circyto import/join integration without genome-scale or external-tool execution." >&2
exit 1
