# Reproducibility Checklist

This checklist captures the recommended release-hygiene artifacts for `circyto` benchmark and manuscript work.

## Recommended supplementary artifacts

- per-workdir `workflow_summary.json`
- benchmark table from `circyto summarize-benchmark`
- `anndata/circ_counts.h5ad`
- `mudata/full_length.h5mu` when available
- `dna/cnv.h5ad` for scRR CNV imports
- `dna/scrr_cnv_import_summary.json` or the actual `scrr_cnv_import_summary.json` path emitted by `import-scrr-cnv`
- `dna_rt/rt.h5ad` for scRR replication timing/state imports
- `dna_rt/scrr_rt_import_summary.json` or the actual `scrr_rt_import_summary.json` path emitted by `import-scrr-rt`
- scRR GSM-to-cell mapping TSV from `build-scrr-cell-map`
- tri-modal RNA+circ+CNV `.h5mu` when produced by `merge-scrr-cnv`
- tri-modal RNA+circ+RT `.h5mu` when produced by `merge-scrr-rt`
- `qc/` summaries
- manifest TSVs used for execution
- command logs for public benchmark runs

## Recommended archive structure

- `benchmarks/`
  - `benchmark_summary.tsv`
  - `benchmark_summary.json`
- `workdirs/`
  - one subdirectory per benchmark
- `manifests/`
  - execution manifests
  - public-download manifests
- `docs/`
  - workflow and protocol notes

## Suggested Zenodo release contents

- versioned source archive
- benchmark summary tables
- selected completed workdir summaries
- representative `h5ad` / `h5mu` artifacts
- processed scRR CNV or RT import summaries and mapping TSVs when DNA modalities are included
- manifest files
- README describing dataset, protocol, and reference assumptions

## Expected benchmark table outputs

From `circyto summarize-benchmark`:

- dataset/workdir identity
- workflow type
- protocol and read layout
- cell counts
- RNA and circ feature counts
- median RNA and circ QC metrics
- `h5mu` presence and size
- CNV or RT modality presence and feature resolution for scRR tri-modal datasets
- cleanup and success status
