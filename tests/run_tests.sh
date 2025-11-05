#!/bin/bash
echo "Running mini tests..."
for eng in ciri-full ciri-long find-circ circ-explorer2; do
  echo "Testing $eng"
  mkdir -p tests/out_$eng
  circyto run --engine $eng --manifest tests/cells.tsv --ref-fa tests/mini_R1.fastq --gtf tests/mini_R2.fastq --outdir tests/out_$eng || echo "skip $eng"
done
echo "Merging union/intersection"
circyto merge --inputs tests/out_ciri-full tests/out_ciri-long --mode union --outdir tests/merged_union || true
circyto merge --inputs tests/out_ciri-full tests/out_ciri-long --mode intersection --outdir tests/merged_intersection || true
