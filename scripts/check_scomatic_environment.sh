#!/usr/bin/env bash
set -euo pipefail

SCOMATIC_DIR="/mnt/d/SComatic"
PYTHON_BIN="${PYTHON_BIN:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scomatic-dir)
      SCOMATIC_DIR="$2"
      shift 2
      ;;
    --python-bin)
      PYTHON_BIN="$2"
      shift 2
      ;;
    -h|--help)
      echo "usage: $0 [--scomatic-dir /mnt/d/SComatic]" >&2
      exit 0
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
    PYTHON_BIN="${CONDA_PREFIX}/bin/python"
  else
    PYTHON_BIN="$(command -v python)"
  fi
fi

if [[ ! -d "${SCOMATIC_DIR}" ]]; then
  echo "[ERROR] SComatic directory not found: ${SCOMATIC_DIR}" >&2
  echo "[ERROR] Expected clone command:" >&2
  echo "        git clone https://github.com/cortes-ciriano-lab/SComatic ${SCOMATIC_DIR}" >&2
  exit 1
fi

PYTHON_VERSION="$("${PYTHON_BIN}" --version 2>&1)"

echo "SComatic directory: ${SCOMATIC_DIR}"
echo "Python: ${PYTHON_VERSION}"
echo "Python executable: ${PYTHON_BIN}"

"${PYTHON_BIN}" - <<'PY'
import sys
if sys.version_info >= (3, 11):
    raise SystemExit(
        "[ERROR] Python >=3.11 is not supported by the pinned SComatic requirements. "
        "Use a dedicated Python 3.10 conda environment."
    )
PY

if ! "${PYTHON_BIN}" -m pip --version >/dev/null 2>&1; then
  echo "[ERROR] pip is not available for ${PYTHON_BIN}" >&2
  exit 1
fi
if ! "${PYTHON_BIN}" -m pip show numpy >/dev/null 2>&1; then
  echo "[ERROR] numpy is not installed in the current environment" >&2
  echo "[ERROR] Install the pinned SComatic Python requirements in a Python 3.10 environment." >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools is missing" >&2
  exit 1
fi
if ! command -v bedtools >/dev/null 2>&1; then
  echo "[ERROR] bedtools is missing" >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript is missing" >&2
  exit 1
fi

echo "numpy:"
"${PYTHON_BIN}" -m pip show numpy | sed -n '1,5p'
echo "samtools: $(samtools --version | head -n 1)"
echo "bedtools: $(bedtools --version | head -n 1)"
if command -v datamash >/dev/null 2>&1; then
  echo "datamash: $(datamash --version | head -n 1)"
else
  echo "datamash: not found (optional)"
fi
echo "Rscript: $(Rscript --version 2>&1)"
Rscript -e 'library(VGAM); cat("VGAM OK\n")'

KEY_SCRIPTS=(
  "${SCOMATIC_DIR}/scripts/BaseCellCounter/BaseCellCounter.py"
  "${SCOMATIC_DIR}/scripts/BaseCellCalling/BaseCellCalling.step1.py"
  "${SCOMATIC_DIR}/scripts/SingleCellGenotype/SingleCellGenotype.py"
)

for script_path in "${KEY_SCRIPTS[@]}"; do
  if [[ ! -f "${script_path}" ]]; then
    echo "[ERROR] required SComatic script missing: ${script_path}" >&2
    exit 1
  fi
  echo "found script: ${script_path}"
done

echo "[INFO] SComatic environment check passed"
