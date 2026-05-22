# Local SComatic chr21 POC

This page defines a technical proof-of-concept path for integrating tiny chr21-scale SComatic-like outputs into `circyto`.

Scope:

- local mini-reference only
- no full hg38 analysis
- no full DNA variant calling
- no genome-scale SComatic execution inside `circyto`

## Command shape

```bash
bash scripts/run_scomatic_chr21_poc.sh \
  --workdir WORKDIR \
  --reference ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --outdir OUTDIR \
  --synthetic
```

## Required local inputs

- a completed local chr21 mini `circyto` workdir
- `ref/chr21.fa`
- `ref/chr21.gtf`

If a real BAM/SAM-backed local POC is attempted later:

- `samtools` must be installed
- `SComatic` must be installed explicitly

## Current behavior

### `--synthetic`

Creates:

- `scomatic_candidate_summary.tsv`
- `scomatic_poc_summary.json`

The candidate table is already normalized into the current `circyto` schema:

- `variant_id`
- `cell_id`
- `chrom`
- `pos`
- `ref`
- `alt`
- `gene`
- `filter_status`
- `candidate_variant_class`
- `read_support`
- `vaf`
- `caller`

### Non-synthetic mode

The script currently:

- discovers an existing BAM/SAM if present
- checks for `samtools`
- checks for `SComatic`
- fails clearly if the environment is not ready

Real tiny chr21 SComatic execution is intentionally not enabled in this revision. That keeps the current POC technical and local, without pretending to provide validated biological variant calls.

## Integration with circyto

Import the normalized candidate table:

```bash
circyto import-dna-snv-summary \
  --workdir WORKDIR \
  --dna-cell-summary dna_cell_summary.tsv \
  --scomatic-candidate-summary OUTDIR/scomatic_candidate_summary.tsv
```

Then generate the joined summary:

```bash
circyto summarize-dna-rna-circ \
  --workdir WORKDIR \
  --write-summary
```

## Scientific guardrail

These outputs should be described as:

- `RNA-derived candidate variant signals`

They should not be described as validated somatic mutations without orthogonal DNA evidence.
