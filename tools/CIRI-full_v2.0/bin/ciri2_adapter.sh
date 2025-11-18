#!/usr/bin/env bash
set -euo pipefail

# Env vars from circyto:
#   R1, R2, REF_FA, GTF, OUT_TSV, THREADS
R1="${R1:?}"; R2="${R2:-}"
REF_FA="${REF_FA:?}"; GTF="${GTF:?}"
OUT_TSV="${OUT_TSV:?}"
THREADS="${THREADS:-4}"

OUT_PREFIX="${OUT_TSV%.tsv}"
OUT_DIR="$(dirname "${OUT_PREFIX}")"
OUT_BASENAME="$(basename "${OUT_PREFIX}")"

# Layout:
#   tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CIRI2_DIR="${ROOT_DIR}/CIRI_v2.0.6"
CIRI2_PL="${CIRI2_DIR}/CIRI2.pl"

if [[ ! -x "${CIRI2_PL}" && ! -f "${CIRI2_PL}" ]]; then
  echo "ERROR: CIRI2.pl not found at ${CIRI2_PL}" >&2
  exit 2
fi

# Log + run dir
RUN_DIR="${OUT_DIR}/${OUT_BASENAME}.ciri2_run"
LOG="${OUT_PREFIX}.ciri2.log"
mkdir -p "${RUN_DIR}"
: > "${LOG}"

echo ">>> CIRI2 adapter" | tee -a "${LOG}"
echo "R1=${R1}" | tee -a "${LOG}"
[[ -n "${R2}" ]] && echo "R2=${R2}" | tee -a "${LOG}"
echo "REF_FA=${REF_FA}" | tee -a "${LOG}"
echo "GTF=${GTF}" | tee -a "${LOG}"
echo "OUT_TSV=${OUT_TSV}" | tee -a "${LOG}"
echo "THREADS=${THREADS}" | tee -a "${LOG}"

# --- sanity: perl & bwa
for bin in perl bwa samtools; do
  if ! command -v "$bin" >/dev/null 2>&1; then
    echo "ERROR: '$bin' not found in PATH" | tee -a "${LOG}"
    exit 2
  fi
done

# ensure index exists (bwa index will no-op if already built)
if [[ ! -s "${REF_FA}.bwt" ]]; then
  echo ">>> building BWA index for ${REF_FA}" | tee -a "${LOG}"
  bwa index "${REF_FA}" >> "${LOG}" 2>&1
fi

# decompress FASTQ if needed
TMP_INDIR="$(mktemp -d -p "${OUT_DIR}" "${OUT_BASENAME}.ciri2_input.XXXX")"
trap 'rm -rf "${TMP_INDIR}" 2>/dev/null || true' EXIT

R1_IN="${R1}"
R2_IN="${R2}"

if [[ "${R1}" == *.gz ]]; then
  R1_IN="${TMP_INDIR}/R1.fq"
  zcat "${R1}" > "${R1_IN}"
  if [[ -n "${R2}" ]]; then
    R2_IN="${TMP_INDIR}/R2.fq"
    zcat "${R2}" > "${R2_IN}"
  fi
fi

SAM="${RUN_DIR}/${OUT_BASENAME}.sam"

echo ">>> running bwa mem" | tee -a "${LOG}"
if [[ -n "${R2_IN}" ]]; then
  bwa mem -T 19 -t "${THREADS}" "${REF_FA}" "${R1_IN}" "${R2_IN}" > "${SAM}" 2>> "${LOG}"
else
  bwa mem -T 19 -t "${THREADS}" "${REF_FA}" "${R1_IN}" > "${SAM}" 2>> "${LOG}"
fi

# --- CIRI2 flags
# default: low/zero stringency (-0); allow override via CIRI2_FLAGS
CIRI2_FLAGS_DEFAULT="-0"
CIRI2_FLAGS="${CIRI2_FLAGS:-${CIRI2_FLAGS_DEFAULT}}"

OUT_TXT="${RUN_DIR}/${OUT_BASENAME}.ciri2.txt"

echo ">>> running CIRI2.pl on ${SAM}" | tee -a "${LOG}"
echo "    perl ${CIRI2_PL} -I ${SAM} -O ${OUT_TXT} -F ${REF_FA} -A ${GTF} -T ${THREADS} ${CIRI2_FLAGS}" | tee -a "${LOG}"

if ! perl "${CIRI2_PL}" \
    -I "${SAM}" \
    -O "${OUT_TXT}" \
    -F "${REF_FA}" \
    -A "${GTF}" \
    -T "${THREADS}" \
    ${CIRI2_FLAGS} >> "${LOG}" 2>&1; then
  echo "ERROR: CIRI2.pl failed; see ${LOG}" >&2
  exit 2
fi

if [[ ! -s "${OUT_TXT}" ]]; then
  echo ">>> CIRI2 output empty; writing header-only TSV ${OUT_TSV}" | tee -a "${LOG}"
  printf "circ_id\tchr\tstart\tend\tstrand\tsupport\n" > "${OUT_TSV}"
  exit 0
fi

# --- Normalize CIRI2 output -> OUT_TSV
echo ">>> normalizing CIRI2 output to ${OUT_TSV}" | tee -a "${LOG}"

awk -v OFS='\t' '
BEGIN {
  print "circ_id","chr","start","end","strand","support"
}
NR==1 {
  # header mapping
  for (i=1; i<=NF; i++) {
    key=$i
    gsub(/[^A-Za-z0-9_]/,"",key)
    key=tolower(key)
    h[key]=i
  }
  next
}
{
  chr   = (h["chr"] ? $h["chr"] : $2)
  start = (h["circrnastart"] ? $h["circrnastart"] : (h["start"] ? $h["start"] : $3))
  end   = (h["circrnaend"]   ? $h["circrnaend"]   : (h["end"]   ? $h["end"]   : $4))
  strand= (h["strand"] ? $h["strand"] : "+")
  supp  = 1
  if (h["junctionreads"]) supp=$h["junctionreads"]
  else if (h["readnum"])  supp=$h["readnum"]
  else if (h["support"])  supp=$h["support"]

  if (h["circrnaid"]) cid=$h["circrnaid"];
  else if (h["circrna_id"]) cid=$h["circrna_id"];
  else cid=chr ":" start "|" end "|" strand

  print cid, chr, start, end, strand, supp
}
' "${OUT_TXT}" > "${OUT_TSV}"

echo ">>> Wrote normalized TSV: ${OUT_TSV}" | tee -a "${LOG}"
