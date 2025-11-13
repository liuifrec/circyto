# CIRI-full Integration (chr21 Smart-seq2 Example)

## Overview
This guide explains how to run CIRI-full v2.x through `circyto` using a chr21-only Smart-seq2 subset from E-MTAB-6072. It demonstrates:

1. Preparing a manifest
2. Running CIRI-full using `run-manifest`
3. Normalizing outputs into TSV
4. Building a circRNA × cell matrix using `collect`
5. Integration testing

## File Structure Example
```
ref/
  chr21.fa
  chr21.fa.{amb,ann,bwt,pac,sa}
  chr21.gtf

fastq/E-MTAB-6072/
  ERRxxxxxx_1.fastq.gz
  ERRxxxxxx_2.fastq.gz

tools/CIRI-full_v2.0/
  CIRI-full.jar
  bin/ciri_full_adapter.sh

manifest.tsv
```

## Manifest Format
```
cell_id    r1    r2
ERR2139486 fastq/E-MTAB-6072/ERR2139486_1.fastq.gz fastq/E-MTAB-6072/ERR2139486_2.fastq.gz
ERR2139559 fastq/E-MTAB-6072/ERR2139559_1.fastq.gz fastq/E-MTAB-6072/ERR2139559_2.fastq.gz
```

## Running CIRI-full via run-manifest
```
circyto run-manifest \
  --manifest manifest.tsv \
  --outdir work_smartseq2/ciri_full_chr21_all \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --cmd-template 'bash -lc "
    export R1={r1} R2={r2} REF_FA={ref_fa} GTF={gtf} OUT_TSV={out_tsv} THREADS=8 ;
    tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh
  "'
```

## Collecting into a Matrix
```
circyto collect \
  --cirifull-dir work_smartseq2/ciri_full_chr21_all \
  --matrix work_smartseq2/circ_chr21_all.mtx \
  --circ-index work_smartseq2/circ_chr21_all_ids.txt \
  --cell-index work_smartseq2/cell_chr21_all_ids.txt
```

Matrix example:
```
12 16 13
1 1 1
2 2 1
3 3 1
4 4 1
```

## Integration Test
Located at:
```
tests/test_cirifull_chr21_integration.py
```

Run:
```
pytest -vv
```

## TODO
- MatrixMarket header already uses `general`.
- Add `.h5ad` writer
- Add detector comparison table
- Normalize circRNA ID format
