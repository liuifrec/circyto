#!/usr/bin/env bash
set -euo pipefail

# ---- paths
REF_FA="ref/genome.fa"
GTF="ref/genes.gtf"

# sanity check
[ -f "$REF_FA" ] || { echo "Missing $REF_FA"; exit 1; }
[ -f "$GTF" ]    || { echo "Missing $GTF"; exit 1; }

mkdir -p data/E-MTAB-6072 data/PBMC_v3 work_smartseq2 work_pbmc_v3 tests

echo "==> Building manifests (idempotent)"
# Smart-seq2 (5 cells)
cat > data/E-MTAB-6072/manifest.tsv <<'TSV'
cell_id	r1	r2
cell1	data/E-MTAB-6072/ERR2140608_1.fastq.gz	data/E-MTAB-6072/ERR2140608_2.fastq.gz
cell2	data/E-MTAB-6072/ERR2141250_1.fastq.gz	data/E-MTAB-6072/ERR2141250_2.fastq.gz
cell3	data/E-MTAB-6072/ERR2141294_1.fastq.gz	data/E-MTAB-6072/ERR2141294_2.fastq.gz
cell4	data/E-MTAB-6072/ERR2141406_1.fastq.gz	data/E-MTAB-6072/ERR2141406_2.fastq.gz
cell5	data/E-MTAB-6072/ERR2142696_1.fastq.gz	data/E-MTAB-6072/ERR2142696_2.fastq.gz
TSV

# PBMC v3 (10x 5′ mini subset)
cat > data/PBMC_v3/manifest.tsv <<'TSV'
cell_id	r1	r2
PBMCv3_subset	data/PBMC_v3/PBMCv3_R1_subset.fastq.gz	data/PBMC_v3/PBMCv3_R2_subset.fastq.gz
TSV

echo "==> Running smoke test: Smart-seq2 (plate/full-length)"
circyto make \
  --manifest data/E-MTAB-6072/manifest.tsv \
  --ref-fa "$REF_FA" \
  --gtf "$GTF" \
  --outdir work_smartseq2 \
  --threads 2 \
  --cmd-template 'bash -lc "printf \"circ_id\tchr\tstart\tend\tstrand\tsupport\n\" > {out_tsv}"'

echo "==> Running smoke test: PBMC v3 (10x 5′)"
circyto make \
  --manifest data/PBMC_v3/manifest.tsv \
  --ref-fa "$REF_FA" \
  --gtf "$GTF" \
  --outdir work_pbmc_v3 \
  --threads 2 \
  --cmd-template 'bash -lc "printf \"circ_id\tchr\tstart\tend\tstrand\tsupport\n\" > {out_tsv}"'

echo "==> Collecting matrix from Smart-seq2 outputs"
circyto collect \
  --cirifull-dir work_smartseq2 \
  --matrix tests/smartseq2.mtx \
  --circ-index tests/smartseq2_circ.txt \
  --cell-index tests/smartseq2_cell.txt

echo "==> Converting to loom (will write an .empty.tsv if matrix is empty)"
circyto convert \
  --matrix tests/smartseq2.mtx \
  --circ-index tests/smartseq2_circ.txt \
  --cell-index tests/smartseq2_cell.txt \
  --loom tests/smartseq2.loom || true

echo "All done ✅"
