#!/usr/bin/env bash
set -euo pipefail

echo "🔍 Testing circyto CLI basic commands..."

# Check main help
circyto --help > /dev/null

# Check subcommands exist
circyto run --help > /dev/null
circyto collect --help > /dev/null
circyto convert --help > /dev/null

# Create dummy FASTQ files and manifest
mkdir -p tests/out_ciri-full
printf "circ_id\tchr\tstart\tend\tstrand\tsupport\ncircA\tchr1\t100\t200\t+\t3\n" > tests/out_ciri-full/Batch1.tsv

# Run collection
circyto collect \
  --cirifull-dir tests/out_ciri-full \
  --matrix tests/circ.mtx \
  --circ-index tests/circ_ids.txt \
  --cell-index tests/cell_ids.txt

# Run conversion
circyto convert \
  --matrix tests/circ.mtx \
  --circ-index tests/circ_ids.txt \
  --cell-index tests/cell_ids.txt \
  --loom tests/circ_full.loom

echo "✅ circyto CLI test finished OK!"
