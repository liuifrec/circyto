#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="scomatic"
SCOMATIC_DIR="/mnt/d/SComatic"
RUN_SETUP=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-name)
      ENV_NAME="$2"
      shift 2
      ;;
    --scomatic-dir)
      SCOMATIC_DIR="$2"
      shift 2
      ;;
    --run)
      RUN_SETUP=1
      shift
      ;;
    -h|--help)
      cat <<EOF
usage: $0 [--env-name NAME] [--scomatic-dir DIR] [--run]

Without --run, print the recommended local/WSL SComatic setup commands.
With --run, execute the conda-based setup in the current shell context.
EOF
      exit 0
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

print_plan() {
  cat <<EOF
Recommended local SComatic environment setup:

1. Create a dedicated conda environment pinned to Python 3.10:
   conda create -n ${ENV_NAME} python=3.10 -y

2. Activate the environment:
   conda activate ${ENV_NAME}

3. Install bioinformatics and R dependencies from conda-forge/bioconda:
   conda install -c conda-forge -c bioconda \\
     samtools bedtools datamash r-base r-vgam -y

4. Install SComatic Python requirements:
   pip install -r ${SCOMATIC_DIR}/requirements.txt

5. Confirm the R-side dependency is importable from the environment:
   Rscript -e 'library(VGAM); cat("VGAM OK\\n")'

6. Optional environment check:
   bash scripts/check_scomatic_environment.sh --scomatic-dir ${SCOMATIC_DIR}

If ${SCOMATIC_DIR} does not exist, clone the official repository first:
   git clone https://github.com/cortes-ciriano-lab/SComatic ${SCOMATIC_DIR}
EOF
}

if [[ "${RUN_SETUP}" -eq 0 ]]; then
  print_plan
  exit 0
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda is required for --run setup mode." >&2
  echo "[ERROR] Install Miniconda/Mambaforge first, or rerun without --run to print commands only." >&2
  exit 1
fi

if [[ ! -d "${SCOMATIC_DIR}" ]]; then
  echo "[ERROR] SComatic directory not found: ${SCOMATIC_DIR}" >&2
  echo "[ERROR] Clone the official repository first:" >&2
  echo "        git clone https://github.com/cortes-ciriano-lab/SComatic ${SCOMATIC_DIR}" >&2
  exit 1
fi

eval "$(conda shell.bash hook)"

if ! conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
  conda create -n "${ENV_NAME}" python=3.10 -y
fi

conda activate "${ENV_NAME}"
conda install -c conda-forge -c bioconda samtools bedtools datamash r-base r-vgam -y
pip install -r "${SCOMATIC_DIR}/requirements.txt"
Rscript -e 'library(VGAM); cat("VGAM OK\n")'

echo "[INFO] SComatic local environment setup completed for ${ENV_NAME}"
echo "[INFO] Verify with: bash scripts/check_scomatic_environment.sh --scomatic-dir ${SCOMATIC_DIR}"
