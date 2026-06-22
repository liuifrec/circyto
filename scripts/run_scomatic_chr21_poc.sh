#!/usr/bin/env bash
set -euo pipefail

WORKDIR=""
REFERENCE=""
GTF=""
OUTDIR=""
SYNTHETIC=0
REAL_SMOKE=0
SCOMATIC_DIR="/mnt/d/SComatic"
PYTHON_BIN="${PYTHON_BIN:-}"

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
    --real-smoke)
      REAL_SMOKE=1
      shift
      ;;
    --scomatic-dir)
      SCOMATIC_DIR="$2"
      shift 2
      ;;
    --python-bin)
      PYTHON_BIN="$2"
      shift 2
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${WORKDIR}" || -z "${REFERENCE}" || -z "${GTF}" || -z "${OUTDIR}" ]]; then
  echo "usage: $0 --workdir WORKDIR --reference ref/chr21.fa --gtf ref/chr21.gtf --outdir OUTDIR [--synthetic|--real-smoke] [--scomatic-dir /mnt/d/SComatic] [--python-bin /path/to/python]" >&2
  exit 1
fi

if [[ "${SYNTHETIC}" -eq 1 && "${REAL_SMOKE}" -eq 1 ]]; then
  echo "[ERROR] choose only one mode: --synthetic or --real-smoke" >&2
  exit 1
fi

if [[ "${SYNTHETIC}" -eq 0 && "${REAL_SMOKE}" -eq 0 ]]; then
  echo "[ERROR] explicitly choose one mode: --synthetic or --real-smoke" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
    PYTHON_BIN="${CONDA_PREFIX}/bin/python"
  elif [[ -x "$(cd "$(dirname "$0")/.." && pwd)/.venv/bin/python" ]]; then
    PYTHON_BIN="$(cd "$(dirname "$0")/.." && pwd)/.venv/bin/python"
  else
    PYTHON_BIN="$(command -v python)"
  fi
fi

if [[ "${SYNTHETIC}" -eq 1 ]]; then
  "${PYTHON_BIN}" - "$WORKDIR" "$REFERENCE" "$GTF" "$OUTDIR" <<'PY'
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

if [[ "${REAL_SMOKE}" -ne 1 ]]; then
  echo "[ERROR] non-synthetic mode currently requires --real-smoke" >&2
  exit 1
fi

ALIGNMENT_PATH="$(find "${WORKDIR}" -type f \( -name '*.bam' -o -name '*.sam' \) | head -n 1 || true)"
if [[ -z "${ALIGNMENT_PATH}" ]]; then
  echo "[ERROR] no BAM/SAM discovered under ${WORKDIR}" >&2
  exit 1
fi

if [[ ! -d "${SCOMATIC_DIR}" ]]; then
  echo "[ERROR] SComatic directory not found: ${SCOMATIC_DIR}" >&2
  echo "[ERROR] Clone the official repository or pass --scomatic-dir PATH." >&2
  echo "[ERROR] Suggested clone: git clone https://github.com/cortes-ciriano-lab/SComatic ${SCOMATIC_DIR}" >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools is required for real-smoke chr21 SComatic POC when BAM/SAM conversion or indexing is needed" >&2
  exit 1
fi

CHECK_SCRIPT="$(dirname "$0")/check_scomatic_environment.sh"
if [[ ! -f "${CHECK_SCRIPT}" ]]; then
  echo "[ERROR] environment checker not found: ${CHECK_SCRIPT}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}/logs"
ENV_LOG="${OUTDIR}/logs/scomatic_environment_check.log"
bash "${CHECK_SCRIPT}" --scomatic-dir "${SCOMATIC_DIR}" --python-bin "${PYTHON_BIN}" >"${ENV_LOG}" 2>&1 || {
  cat "${ENV_LOG}" >&2
  exit 1
}

PREPARED_BAM="${OUTDIR}/chr21_input.bam"
PREPARED_BAI="${PREPARED_BAM}.bai"
if [[ "${ALIGNMENT_PATH}" == *.sam ]]; then
  samtools view -bS "${ALIGNMENT_PATH}" | samtools sort -o "${PREPARED_BAM}"
else
  cp "${ALIGNMENT_PATH}" "${PREPARED_BAM}"
fi
samtools index "${PREPARED_BAM}"

BASECOUNTER_HELP_LOG="${OUTDIR}/logs/basecellcounter_help.log"
"${PYTHON_BIN}" "${SCOMATIC_DIR}/scripts/BaseCellCounter/BaseCellCounter.py" --help >"${BASECOUNTER_HELP_LOG}" 2>&1 || {
  echo "[ERROR] Failed to run SComatic BaseCellCounter help smoke test." >&2
  echo "[ERROR] See log: ${BASECOUNTER_HELP_LOG}" >&2
  exit 1
}

"${PYTHON_BIN}" - "$WORKDIR" "$REFERENCE" "$GTF" "$OUTDIR" "$SCOMATIC_DIR" "$ALIGNMENT_PATH" "$PREPARED_BAM" "$PREPARED_BAI" "$ENV_LOG" "$BASECOUNTER_HELP_LOG" <<'PY'
from pathlib import Path
import json
import sys

summary = {
    "mode": "real-smoke",
    "workdir": str(Path(sys.argv[1]).resolve()),
    "reference": str(Path(sys.argv[2]).resolve()),
    "gtf": str(Path(sys.argv[3]).resolve()),
    "outdir": str(Path(sys.argv[4]).resolve()),
    "scomatic_dir": str(Path(sys.argv[5]).resolve()),
    "discovered_alignment": str(Path(sys.argv[6]).resolve()),
    "prepared_bam": str(Path(sys.argv[7]).resolve()),
    "prepared_bai": str(Path(sys.argv[8]).resolve()),
    "environment_log": str(Path(sys.argv[9]).resolve()),
    "basecellcounter_help_log": str(Path(sys.argv[10]).resolve()),
    "note": "Environment smoke test only. No genome-scale SComatic calling was run.",
    "terminology_note": "Any future SComatic outputs must be treated as RNA-derived candidate variant signals unless orthogonal DNA validation exists.",
}
summary_path = Path(sys.argv[4]) / "scomatic_poc_summary.json"
summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2, sort_keys=True))
PY
