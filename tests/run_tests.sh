#!/usr/bin/env bash
# tests/run_tests.sh
export PATH="tools:$PATH"
command -v CIRI-full >/dev/null || { echo "CIRI-full not found (wrapper expected at tools/CIRI-full)"; exit 1; }

set -euo pipefail

# --------------------------------------------------------------------
# 0. Self-bootstrap: ensure CIRI-full binary exists and is executable
# --------------------------------------------------------------------
CIRI_DIR="tools/CIRI-full"
CIRI_BIN="$CIRI_DIR/CIRI-full"

if ! command -v CIRI-full >/dev/null 2>&1; then
  echo "==> CIRI-full not found; bootstrapping local copy under $CIRI_DIR ..."
  mkdir -p "$CIRI_DIR"
  cd "$CIRI_DIR"
  # Download official package if missing
  if [[ ! -f "CIRI-full.pl" ]]; then
    echo "Downloading CIRI-full release from GitHub..."
    curl -L -O https://github.com/YangLab/CIRI-full/archive/refs/heads/master.zip
    unzip master.zip >/dev/null 2>&1 || true
    mv CIRI-full-master/* . || true
    rm -rf CIRI-full-master master.zip
  fi

  # Build wrapper launcher
  cat > CIRI-full <<'PL'
#!/usr/bin/env bash
# Simple wrapper to invoke CIRI-full.pl in this directory
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec perl "$DIR/CIRI-full.pl" "$@"
PL
  chmod +x CIRI-full

  cd - >/dev/null
  export PATH="$PWD/$CIRI_DIR:$PATH"
  echo "CIRI-full installed locally at $CIRI_BIN"
else
  echo "==> CIRI-full already available at $(command -v CIRI-full)"
fi

# --------------------------------------------------------------------
# 1. Paths and configuration
# --------------------------------------------------------------------
FASTQ_DIR="fastq/E-MTAB-6072"
REF_FA="ref/genome.fa"
GTF="ref/genes.gtf"
OUTDIR="work_smartseq2"
THREADS="${THREADS:-4}"
TARGET_LIST="tests/target_cells.txt"
CHR21_FA="ref/chr21.fa"
CHR21_GTF="ref/chr21.gtf"
CIRCSC_LIST="tests/circsc_chr21_list.txt"

# Choose references (chr21 if present)
if [[ -f "$CHR21_FA" && -f "$CHR21_GTF" ]]; then
  USE_FA="$CHR21_FA"; USE_GTF="$CHR21_GTF"
  echo ">> Using chr21-only references: $USE_FA , $USE_GTF"
else
  USE_FA="$REF_FA";  USE_GTF="$GTF"
  echo ">> Using full references: $USE_FA , $USE_GTF"
fi

[[ -f "$USE_FA" ]] || { echo "Missing FASTA: $USE_FA"; exit 1; }
[[ -f "$USE_GTF" ]] || { echo "Missing GTF:   $USE_GTF"; exit 1; }

command -v CIRI-full >/dev/null 2>&1 || {
  echo "ERROR: CIRI-full not on PATH. Make sure it’s installed and callable as 'CIRI-full'."; exit 1;
}

# ---------------------------
# Build a manifest from FASTQ_DIR
# ---------------------------
echo "==> Scanning FASTQs in $FASTQ_DIR"
mkdir -p "$OUTDIR"
MAN="$OUTDIR/manifest.tsv"
echo -e "cell_id\tr1\tr2" > "$MAN"

# Find R1 files and pair with R2
while IFS= read -r R1; do
  base="$(basename "$R1" | sed 's/_1\.fastq\.gz$//')"
  R2="$FASTQ_DIR/${base}_2.fastq.gz"
  [[ -f "$R2" ]] || { echo "WARN: Missing R2 for $R1, skipping"; continue; }
  cell_id="$base"  # Example: ERR2140608
  echo -e "${cell_id}\t${R1}\t${R2}" >> "$MAN"
done < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*_1.fastq.gz" | sort)

echo "Manifest has $(($(wc -l < "$MAN")-1)) rows."

# ---------------------------
# Optionally restrict to target ERR IDs
# ---------------------------
if [[ -s "$TARGET_LIST" ]]; then
  echo "==> Restricting to targets listed in $TARGET_LIST"
  awk 'NR==FNR{t[$1]=1;next} NR==1 || t[$1]' "$TARGET_LIST" "$MAN" > "$MAN.tmp" && mv "$MAN.tmp" "$MAN"
  echo "Filtered manifest now has $(($(wc -l < "$MAN")-1)) rows."
fi

# Quick bail if nothing to run
if [[ $(wc -l < "$MAN") -le 1 ]]; then
  echo "No rows in manifest after filtering. Nothing to do."
  exit 0
fi

# ---------------------------
# Run real CIRI-full via circyto
# ---------------------------
echo "==> Running CIRI-full via circyto (manifest mode)"
# Real template: change -p if you want per-job multithreading, but circyto already parallelizes at batch level later.
CMD_TPL='CIRI-full -M CIRI_Full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv} -p 1'

circyto make \
  --manifest "$MAN" \
  --ref-fa "$USE_FA" \
  --gtf "$USE_GTF" \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  --cmd-template "$CMD_TPL"

# ---------------------------
# Collect matrix
# ---------------------------
echo "==> Collecting circ×cell matrix"
MTX="tests/smartseq2.mtx"
CIRC_IDX="tests/smartseq2_circ.txt"
CELL_IDX="tests/smartseq2_cell.txt"
mkdir -p tests

circyto collect \
  --cirifull-dir "$OUTDIR/cirifull_out" \
  --matrix "$MTX" \
  --circ-index "$CIRC_IDX" \
  --cell-index "$CELL_IDX"

# ---------------------------
# Optional: overlap with circSC list (if provided)
# ---------------------------
if [[ -s "$CIRCSC_LIST" ]]; then
  echo "==> Checking overlap with circSC list: $CIRCSC_LIST"
  OVER="tests/overlap_with_circsc.txt"
  # Grep all circ IDs reported by circyto in outputs and intersect with your list
  grep -h -v '^circ_id' "$OUTDIR"/cirifull_out/*.tsv | cut -f1 | sort -u > tests/detected_circ_ids.txt
  comm -12 <(sort -u "$CIRCSC_LIST") tests/detected_circ_ids.txt > "$OVER" || true
  echo "Overlap circ IDs written to $OVER (lines=$(wc -l < "$OVER"))"
fi

# ---------------------------
# Summary
# ---------------------------
echo "==> Summary"
echo "MTX:      $MTX"
echo "CIRC IDX: $CIRC_IDX  (#=$(wc -l < "$CIRC_IDX" 2>/dev/null || echo 0))"
echo "CELL IDX: $CELL_IDX  (#=$(wc -l < "$CELL_IDX" 2>/dev/null || echo 0))"
if [[ -f tests/overlap_with_circsc.txt ]]; then
  echo "Overlap lines: $(wc -l < tests/overlap_with_circsc.txt)"
fi

echo "Done."
