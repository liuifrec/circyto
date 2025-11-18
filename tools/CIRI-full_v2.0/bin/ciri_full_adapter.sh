#!/usr/bin/env bash
set -euo pipefail

R1="${R1:?}"; R2="${R2:-}"
REF_FA="${REF_FA:?}"; GTF="${GTF:?}"
OUT_TSV="${OUT_TSV:?}"
THREADS="${THREADS:-4}"
CIRI_EXTRA_FLAGS="${CIRI_EXTRA_FLAGS:-}"

OUT_PREFIX="${OUT_TSV%.tsv}"
OUT_DIR="$(dirname "${OUT_PREFIX}")"
OUT_BASENAME="$(basename "${OUT_PREFIX}")"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JAR=""
if [[ -f "${ROOT_DIR}/CIRI-full.jar" ]]; then
  JAR="${ROOT_DIR}/CIRI-full.jar"
elif [[ -f "${ROOT_DIR}/CIRI_Full.jar" ]]; then
  JAR="${ROOT_DIR}/CIRI_Full.jar"
else
  echo "ERROR: CIRI-full jar not found in ${ROOT_DIR}" >&2
  exit 2
fi

JAVA=${JAVA:-java}
LOG="${OUT_PREFIX}.ciri_full.log"
: > "${LOG}"

# --- sanity checks
echo ">>> Sanity" | tee -a "${LOG}"
for bin in bwa samtools; do
  if ! command -v "$bin" >/dev/null 2>&1; then
    echo "ERROR: '$bin' not found in PATH" | tee -a "${LOG}"
    exit 2
  fi
  echo "found $(which "$bin")" | tee -a "${LOG}"
done
for suf in "" .amb .ann .bwt .pac .sa; do
  f="${REF_FA}${suf}"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing reference/index file: $f" | tee -a "${LOG}"
    exit 2
  fi
done
[[ -s "${GTF}" ]] || { echo "ERROR: missing ${GTF}" | tee -a "${LOG}"; exit 2; }

echo ">>> CIRI-full Pipeline run" | tee -a "${LOG}"
echo "R1=${R1}" | tee -a "${LOG}"
[[ -n "${R2}" ]] && echo "R2=${R2}" | tee -a "${LOG}"
echo "REF_FA=${REF_FA}" | tee -a "${LOG}"
echo "GTF=${GTF}" | tee -a "${LOG}"
echo "OUT_DIR=${OUT_DIR} OUT_BASENAME=${OUT_BASENAME} THREADS=${THREADS}" | tee -a "${LOG}"
[[ -n "${CIRI_EXTRA_FLAGS}" ]] && echo "CIRI_EXTRA_FLAGS=${CIRI_EXTRA_FLAGS}" | tee -a "${LOG}"

RUN_DIR="${OUT_DIR}/${OUT_BASENAME}.ciri_full_run"
rm -rf "${RUN_DIR}" 2>/dev/null || true
# DO NOT mkdir RUN_DIR; the jar creates it.

TMP_INDIR="$(mktemp -d -p "${OUT_DIR}" "${OUT_BASENAME}.input.XXXX")"
trap 'rm -rf "${TMP_INDIR}" 2>/dev/null || true' EXIT

R1_IN="${R1}"; R2_IN="${R2}"
if [[ "${R1}" == *.gz ]]; then
  R1_IN="${TMP_INDIR}/R1.fq"
  zcat "${R1}" > "${R1_IN}"
  if [[ -n "${R2}" ]]; then
    R2_IN="${TMP_INDIR}/R2.fq"
    zcat "${R2}" > "${R2_IN}"
  fi
fi

cmd=( "${JAVA}" -Xmx8g -jar "${JAR}" Pipeline \
      -1 "${R1_IN}" \
      -d "${RUN_DIR}" \
      -o "${OUT_BASENAME}" \
      -r "${REF_FA}" \
      -a "${GTF}" \
      -t "${THREADS}" )

[[ -n "${R2_IN}" ]] && cmd+=( -2 "${R2_IN}" )

# Append extra flags for CIRI if provided (e.g. -0, -low, --no_strigency)
if [[ -n "${CIRI_EXTRA_FLAGS}" ]]; then
  # Let the shell split the flags
  cmd+=( ${CIRI_EXTRA_FLAGS} )
fi

echo ">>> CMD: ${cmd[*]}" | tee -a "${LOG}"
if ! bash -lc "${cmd[*]}" >> "${LOG}" 2>&1; then
  echo "ERROR: CIRI-full Pipeline failed; see ${LOG}" >&2
  echo ">>> RUN_DIR listing after failure:" | tee -a "${LOG}"
  (ls -lah "${RUN_DIR}" || true) | tee -a "${LOG}"
  exit 2
fi

echo ">>> RUN_DIR contents:" | tee -a "${LOG}"
(ls -lah "${RUN_DIR}" || true) | tee -a "${LOG}"

echo ">>> RUN_DIR recursive listing:" | tee -a "${LOG}"
(find "${RUN_DIR}" -maxdepth 3 -type f -printf "%P\t%p\n" || true) | tee -a "${LOG}"

# --- choose a primary output file

candidates=()

# 1) Prefer CIRI_output/*.ciri
while IFS= read -r -d '' f; do
  candidates+=("$f")
done < <(find "${RUN_DIR}" -maxdepth 3 -type f -name "${OUT_BASENAME}*.ciri" -print0)

# 2) If none, look for CIRI-full_output/ ro2/merge outputs
if [[ ${#candidates[@]} -eq 0 ]]; then
  while IFS= read -r -d '' f; do
    candidates+=("$f")
  done < <(find "${RUN_DIR}" -maxdepth 3 -type f \( \
      -name "${OUT_BASENAME}_ro2_info.list" -o \
      -name "${OUT_BASENAME}_merge_circRNA_detail.anno" \
    \) -print0)
fi

if [[ ${#candidates[@]} -eq 0 ]]; then
  echo "ERROR: No CIRI-full output files (.ciri / ro2 / merge) matched. See ${LOG} for full listing." | tee -a "${LOG}"
  exit 2
fi

MAIN_TXT="${candidates[0]}"
echo ">>> Using output: ${MAIN_TXT}" | tee -a "${LOG}"

# --- Normalize -> OUT_TSV
awk -v OFS='\t' '
BEGIN { print "circ_id","chr","start","end","strand","support" }
NR==1 {
  for(i=1;i<=NF;i++){
    k=$i; gsub(/[^A-Za-z0-9_]/,"",k); k=tolower(k); h[k]=i
  }
  next
}
{
  chr   = (h["chr"] ? $h["chr"] : (h["chrom"] ? $h["chrom"] : $1))
  start = (h["circrnastart"] ? $h["circrnastart"] : (h["start"] ? $h["start"] : $2))
  end   = (h["circrnaend"]   ? $h["circrnaend"]   : (h["end"]   ? $h["end"]   : $3))
  strand= (h["strand"] ? $h["strand"] : (h["str"] ? $h["str"] : "+"))
  supp  = 1
  if (h["junctionreads"]) supp=$h["junctionreads"]
  else if (h["readnum"])  supp=$h["readnum"]
  else if (h["support"])  supp=$h["support"]
  if (h["circrnaid"]) cid=$h["circrnaid"]; else cid=chr ":" start "|" end "|" strand
  print cid, chr, start, end, strand, supp
}
' "${MAIN_TXT}" > "${OUT_TSV}"

echo ">>> Wrote normalized TSV: ${OUT_TSV}" | tee -a "${LOG}"
