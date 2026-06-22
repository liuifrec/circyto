# Local SComatic chr21 POC

This page defines a technical proof-of-concept path for integrating tiny chr21-scale SComatic-like outputs into `circyto`.

Scope:

- local mini-reference only
- no full hg38 analysis
- no full DNA variant calling
- no genome-scale SComatic execution inside `circyto`

## A. Synthetic integration test

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
- `SComatic` must be installed explicitly in a dedicated Python 3.10 environment
- `VGAM` must be available to `Rscript`

See [scomatic_environment_setup.md](./scomatic_environment_setup.md).

### Current WSL boundary

Real local WSL SComatic smoke mode is currently blocked on this machine by native conda package instability:

- repeated `Bus error (core dumped)` during minimal dependency installation

For now, use WSL only for:

- `--synthetic`
- `circyto import-dna-snv-summary`
- `circyto summarize-dna-rna-circ`

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

For environment/protocol smoke testing, use:

```bash
bash scripts/run_scomatic_chr21_poc.sh \
  --workdir WORKDIR \
  --reference ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --outdir OUTDIR \
  --scomatic-dir /mnt/d/SComatic \
  --real-smoke
```

## B. Real SComatic environment smoke test

The real-smoke path:

- discovers an existing BAM/SAM if present
- runs the local SComatic environment checker
- rejects Python 3.11+ explicitly
- rejects missing `VGAM` explicitly
- prepares a tiny local BAM/index under the POC output directory
- runs a minimal `BaseCellCounter.py --help` launch smoke test
- writes `scomatic_poc_summary.json`

This is still not a production variant-calling run. It is only an environment/protocol compatibility check.

At present, this real-smoke path should be attempted on:

- an HPC/server conda environment
- or a container / `mamba` / `micromamba` environment

Do not keep retrying local WSL native installs after a repeated `Bus error`.

## C. Future production genome-scale run

Genome-scale SComatic use should remain clearly separate from this local POC. The HAP1 batch10 server smoke validated the technical path through BaseCellCounter, Step1, Step2, and result normalization, but candidate interpretation remains exploratory:

- use the dedicated SComatic environment
- use the full reference resources expected by SComatic
- validate resource compatibility carefully
- prefer HPC/server or containerized execution if native WSL package installation is unstable
- keep terminology conservative:
  - `RNA-derived candidate variant signals`
  - RNA-derived candidate variant signals without orthogonal DNA evidence

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

They should not be described as orthogonally confirmed somatic variants without orthogonal DNA evidence.
