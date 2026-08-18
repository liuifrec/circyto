#!/usr/bin/env bash
set -euo pipefail

R1="${R1:?}"; R2="${R2:-}"
REF_FA="${REF_FA:?}"; GTF="${GTF:?}"
OUT_TSV="${OUT_TSV:?}"
THREADS="${THREADS:-4}"
CIRI_EXTRA_FLAGS="${CIRI_EXTRA_FLAGS:-}"
READ_LAYOUT_ENV="${CIRCYTO_READ_LAYOUT:-}"
CIRI2_BWA_MEM_FLAGS="${CIRI2_BWA_MEM_FLAGS:-}"

OUT_PREFIX="${OUT_TSV%.tsv}"
OUT_DIR="$(dirname "${OUT_PREFIX}")"
OUT_BASENAME="$(basename "${OUT_PREFIX}")"
RUN_DIR="${OUT_DIR}/${OUT_BASENAME}.ciri_full_run"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${ROOT_DIR}/bin"
JAR="${CIRCYTO_CIRI_FULL_JAR:-${CIRI_FULL_JAR:-}}"
if [[ -n "${JAR}" && -f "${JAR}" ]]; then
  :
elif [[ -f "${ROOT_DIR}/CIRI-full.jar" ]]; then
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

detect_layout() {
  if [[ "${READ_LAYOUT_ENV}" == "single" || "${READ_LAYOUT_ENV}" == "paired" ]]; then
    echo "${READ_LAYOUT_ENV}"
    return
  fi
  if [[ -n "${R2}" ]]; then
    echo "paired"
  else
    echo "single"
  fi
}

tail_log() {
  tail -n 40 "${LOG}" 2>/dev/null || true
}

print_cmd() {
  printf '%q ' "$@"
  printf '\n'
}

normalize_primary_output() {
  local main_txt="$1"
  awk -v OFS='\t' '
BEGIN {
  print "circ_id","chr","start","end","strand","support"
}
NR==1 {
  for (i=1; i<=NF; i++) {
    key = $i
    gsub(/[^A-Za-z0-9_]/, "", key)
    key = tolower(key)
    h[key] = i
  }
  next
}
{
  chr = (h["chr"] ? $h["chr"] :
         (h["chrom"] ? $h["chrom"] :
          $2))

  start = (h["circrna_start"] ? $h["circrna_start"] :
           (h["circrnastart"] ? $h["circrnastart"] :
           (h["start"] ? $h["start"] :
            $3)))

  end = (h["circrna_end"] ? $h["circrna_end"] :
         (h["circrnaend"] ? $h["circrnaend"] :
         (h["end"] ? $h["end"] :
          $4)))

  strand = (h["strand"] ? $h["strand"] :
            (h["str"] ? $h["str"] : "+"))

  supp = 1
  if (h["junction_reads"])      supp = $h["junction_reads"]
  else if (h["junctionreads"])  supp = $h["junctionreads"]
  else if (h["readnum"])        supp = $h["readnum"]
  else if (h["support"])        supp = $h["support"]

  if (h["circrna_id"])          cid = $h["circrna_id"]
  else if (h["circrnaid"])      cid = $h["circrnaid"]
  else                          cid = chr ":" start "|" end "|" strand

  print cid, chr, start, end, strand, supp
}
' "${main_txt}" > "${OUT_TSV}"
}

run_single_end_fallback() {
  local ciri2_pl="${BIN_DIR}/CIRI_v2.0.6/CIRI2.pl"
  local tmp_indir r1_in sam out_txt cmd_string

  if [[ ! -f "${ciri2_pl}" ]]; then
    echo "ERROR: CIRI2.pl not found at ${ciri2_pl}" | tee -a "${LOG}"
    exit 2
  fi

  mkdir -p "${RUN_DIR}"
  tmp_indir="$(mktemp -d -p "${OUT_DIR}" "${OUT_BASENAME}.input.XXXX")"
  trap 'rm -rf "'"${tmp_indir}"'" 2>/dev/null || true' EXIT

  r1_in="${R1}"
  if [[ "${R1}" == *.gz ]]; then
    r1_in="${tmp_indir}/R1.fq"
    echo ">>> zcat+sanitize R1 -> ${r1_in}" | tee -a "${LOG}"
    zcat "${R1}" | awk '
      NR%4==1 { print $1; next }
      NR%4==3 { print "+"; next }
      { print }
    ' > "${r1_in}"
  fi

  sam="${RUN_DIR}/${OUT_BASENAME}.sam"
  out_txt="${RUN_DIR}/${OUT_BASENAME}.ciri2.txt"

  if [[ -z "${CIRI2_BWA_MEM_FLAGS}" ]]; then
    CIRI2_BWA_MEM_FLAGS="-T 19"
  fi

  echo ">>> Execution mode: single-end fallback via bundled CIRI2" | tee -a "${LOG}"
  echo ">>> bwa inputs: R1_IN=${r1_in}" | tee -a "${LOG}"
  echo ">>> bwa mem flags: ${CIRI2_BWA_MEM_FLAGS}" | tee -a "${LOG}"

  cmd_string="$(print_cmd bwa mem ${CIRI2_BWA_MEM_FLAGS} -t "${THREADS}" "${REF_FA}" "${r1_in}")"
  echo ">>> CMD[bwa]: ${cmd_string}" | tee -a "${LOG}"
  if ! bwa mem ${CIRI2_BWA_MEM_FLAGS} -t "${THREADS}" "${REF_FA}" "${r1_in}" > "${sam}" 2>> "${LOG}"; then
    echo "ERROR: single-end bwa mem failed" | tee -a "${LOG}"
    exit 2
  fi

  cmd_string="$(print_cmd perl "${ciri2_pl}" -I "${sam}" -O "${out_txt}" -F "${REF_FA}" -A "${GTF}" -T 1 -0 -U 15)"
  echo ">>> CMD[ciri2]: ${cmd_string}" | tee -a "${LOG}"
  if ! perl "${ciri2_pl}" \
      -I "${sam}" \
      -O "${out_txt}" \
      -F "${REF_FA}" \
      -A "${GTF}" \
      -T 1 \
      -0 \
      -U 15 >> "${LOG}" 2>&1; then
    echo "ERROR: bundled CIRI2 single-end fallback failed" | tee -a "${LOG}"
    exit 2
  fi

  if [[ ! -s "${out_txt}" ]]; then
    echo ">>> CIRI2 single-end fallback output empty; writing header-only TSV ${OUT_TSV}" | tee -a "${LOG}"
    printf "circ_id\tchr\tstart\tend\tstrand\tsupport\n" > "${OUT_TSV}"
    return
  fi

  echo ">>> Using output: ${out_txt}" | tee -a "${LOG}"
  normalize_primary_output "${out_txt}"
}

run_paired_end_pipeline() {
  local tmp_indir r1_in r2_in
  local -a cmd
  local main_txt

  rm -rf "${RUN_DIR}" 2>/dev/null || true
  tmp_indir="$(mktemp -d -p "${OUT_DIR}" "${OUT_BASENAME}.input.XXXX")"
  trap 'rm -rf "'"${tmp_indir}"'" 2>/dev/null || true' EXIT

  r1_in="${R1}"
  r2_in="${R2}"
  if [[ "${R1}" == *.gz ]]; then
    r1_in="${tmp_indir}/R1.fq"
    zcat "${R1}" > "${r1_in}"
  fi
  if [[ "${R2}" == *.gz ]]; then
    r2_in="${tmp_indir}/R2.fq"
    zcat "${R2}" > "${r2_in}"
  fi

  cmd=( "${JAVA}" -Xmx8g -jar "${JAR}" Pipeline
        -1 "${r1_in}"
        -d "${RUN_DIR}"
        -o "${OUT_BASENAME}"
        -r "${REF_FA}"
        -a "${GTF}"
        -t "${THREADS}"
        -2 "${r2_in}" )
  if [[ -n "${CIRI_EXTRA_FLAGS}" ]]; then
    cmd+=( ${CIRI_EXTRA_FLAGS} )
  fi

  echo ">>> Execution mode: paired-end CIRI-full Pipeline" | tee -a "${LOG}"
  echo ">>> CMD[pipeline]: $(print_cmd "${cmd[@]}")" | tee -a "${LOG}"
  if ! "${cmd[@]}" >> "${LOG}" 2>&1; then
    echo "ERROR: CIRI-full Pipeline failed" | tee -a "${LOG}"
    echo ">>> RUN_DIR listing after failure:" | tee -a "${LOG}"
    (ls -lah "${RUN_DIR}" || true) | tee -a "${LOG}"
    exit 2
  fi

  echo ">>> RUN_DIR contents:" | tee -a "${LOG}"
  (ls -lah "${RUN_DIR}" || true) | tee -a "${LOG}"

  echo ">>> RUN_DIR recursive listing:" | tee -a "${LOG}"
  (find "${RUN_DIR}" -maxdepth 3 -type f -printf "%P\t%p\n" || true) | tee -a "${LOG}"

  main_txt="$(find "${RUN_DIR}" -maxdepth 3 -type f -name "${OUT_BASENAME}*.ciri" | head -n 1 || true)"
  if [[ -z "${main_txt}" ]]; then
    main_txt="$(find "${RUN_DIR}" -maxdepth 3 -type f \( \
      -name "${OUT_BASENAME}_ro2_info.list" -o \
      -name "${OUT_BASENAME}_merge_circRNA_detail.anno" \
    \) | head -n 1 || true)"
  fi

  if [[ -z "${main_txt}" ]]; then
    echo "ERROR: No CIRI-full output files (.ciri / ro2 / merge) matched." | tee -a "${LOG}"
    exit 2
  fi

  echo ">>> Using output: ${main_txt}" | tee -a "${LOG}"
  normalize_primary_output "${main_txt}"
}

echo ">>> Sanity" | tee -a "${LOG}"
for bin in bwa samtools perl awk zcat; do
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

READ_LAYOUT="$(detect_layout)"
echo ">>> read layout: ${READ_LAYOUT}" | tee -a "${LOG}"
echo "R1=${R1}" | tee -a "${LOG}"
[[ -n "${R2}" ]] && echo "R2=${R2}" | tee -a "${LOG}"
echo "REF_FA=${REF_FA}" | tee -a "${LOG}"
echo "GTF=${GTF}" | tee -a "${LOG}"
echo "OUT_DIR=${OUT_DIR} OUT_BASENAME=${OUT_BASENAME} THREADS=${THREADS}" | tee -a "${LOG}"
[[ -n "${CIRI_EXTRA_FLAGS}" ]] && echo "CIRI_EXTRA_FLAGS=${CIRI_EXTRA_FLAGS}" | tee -a "${LOG}"

if [[ "${READ_LAYOUT}" == "paired" ]]; then
  if [[ -z "${R2}" ]]; then
    echo "ERROR: read layout is paired but R2 is missing." | tee -a "${LOG}"
    exit 2
  fi
  run_paired_end_pipeline
else
  run_single_end_fallback
fi

echo ">>> Wrote normalized TSV: ${OUT_TSV}" | tee -a "${LOG}"
