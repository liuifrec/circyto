#!/usr/bin/env bash
set -euo pipefail

# ===================== CONFIG (override with env) =====================
FASTQ_DIR="${FASTQ_DIR:-fastq/E-MTAB-6072}"
OUTROOT="${OUTROOT:-work_smartseq2/ciri_full_out}"

# Use the small chr21 reference you prepared
REF_FA="${REF_FA:-ref/chr21.fa}"
GTF="${GTF:-ref/chr21.gtf}"          # set empty to skip annotation: GTF=""

THREADS="${THREADS:-8}"
WRAPPER="${WRAPPER:-tools/CIRI-full}"
MANIFEST="${MANIFEST:-tmp/E-MTAB-6072_manifest.tsv}"
# =================================================================

die(){ echo "[run_tests] $*" >&2; exit 1; }
note(){ echo -e "\n[run_tests] $*\n"; }

# --- sanity checks ---
command -v java >/dev/null 2>&1 || die "java not found in PATH"
command -v bwa  >/dev/null 2>&1 || die "bwa not found in PATH"
[[ -x "$WRAPPER" ]] || die "wrapper not executable: $WRAPPER"
[[ -d "$FASTQ_DIR" ]] || die "FASTQ dir not found: $FASTQ_DIR"
[[ -f "$REF_FA" ]] || die "Reference FASTA not found: $REF_FA"
[[ -n "${GTF:-}" ]] && [[ ! -f "$GTF" ]] && die "GTF not found: $GTF"

mkdir -p "$(dirname "$MANIFEST")" "$OUTROOT"

# --- index reference if needed ---
if [[ ! -f "${REF_FA}.bwt" ]]; then
  note "Indexing reference with bwa..."
  bwa index "$REF_FA"
fi

# --- build manifest from FASTQ pairs ---
note "Scanning $FASTQ_DIR for R1/R2 pairs..."
tmp_manifest="$(mktemp)"
echo -e "sample\tr1\tr2" > "$tmp_manifest"

pairs_found=0
# patterns: *_R1*.fastq.gz, *_R1*.fq.gz, *_1.fastq.gz, *_1.fq.gz, *.R1.fastq.gz
shopt -s nullglob
for r1 in \
  "$FASTQ_DIR"/*_R1*.fastq* \
  "$FASTQ_DIR"/*_R1*.fq*    \
  "$FASTQ_DIR"/*_1.fastq*   \
  "$FASTQ_DIR"/*_1.fq*      \
  "$FASTQ_DIR"/*.R1.fastq*  \
  "$FASTQ_DIR"/*.R1.fq*     ;
do
  [[ -e "$r1" ]] || continue
  base="$(basename "$r1")"

  # candidate R2 names
  r2=""
  for c in \
    "${r1/_R1/_R2}" \
    "${r1/_1./_2.}" \
    "${r1/.R1./.R2.}" \
    "${r1/_R1_/_R2_}" \
    "${r1/_1_/ _2_}";
  do
    [[ -f "$c" ]] && { r2="$c"; break; }
  done
  [[ -z "$r2" ]] && { echo "[run_tests] warn: no R2 for $r1 — skipping" >&2; continue; }

  # sample name = filename minus R1 marker and extensions
  s="$base"
  s="${s/_R1/}" ; s="${s/.R1/}" ; s="${s/_1/}"
  s="${s%%.fastq*}" ; s="${s%%.fq*}"

  echo -e "${s}\t${r1}\t${r2}" >> "$tmp_manifest"
  ((pairs_found++))
done
shopt -u nullglob


[[ $pairs_found -gt 0 ]] || die "No R1/R2 pairs detected under $FASTQ_DIR"

mv "$tmp_manifest" "$MANIFEST"
note "Manifest written to $MANIFEST (pairs: $pairs_found)"
# Uncomment to view:
# column -t -s$'\t' "$MANIFEST" | sed -e 's/^/[manifest] /'

# --- run wrapper in manifest mode ---
note "Running CIRI-full on E-MTAB-6072 (manifest mode)…"
cmd=( "$WRAPPER" manifest "$MANIFEST" --ref "$REF_FA" --outdir "$OUTROOT" --threads "$THREADS" )
[[ -n "${GTF:-}" ]] && cmd+=( --gtf "$GTF" )
echo "[run_tests] ${cmd[*]}"
"${cmd[@]}"

# --- validate outputs for each sample ---
note "Validating outputs…"
fails=0; oks=0
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample r1 r2; do
  sdir="${OUTROOT}/${sample}"
  for d in CIRI_output CIRI-AS_output CIRI-full_output sam; do
    [[ -d "${sdir}/${d}" ]] || { echo "[run_tests] MISSING dir: ${sdir}/${d}"; fails=$((fails+1)); continue; }
  done
  detail="${sdir}/CIRI-full_output/${sample}_merge_circRNA_detail.anno"
  if [[ -s "$detail" ]]; then
    echo "[run_tests] ✅ ${sample}: ok"
    oks=$((oks+1))
  else
    echo "[run_tests] ❌ ${sample}: empty/missing ${detail}"
    fails=$((fails+1))
  fi
done

echo
echo "[run_tests] Summary → OK: ${oks} | FAIL: ${fails}"
[[ "${fails}" -eq 0 ]] || exit 1

echo "[run_tests] ✅ All E-MTAB-6072 samples completed successfully."
