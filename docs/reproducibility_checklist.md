# Reproducibility Checklist

This checklist captures the recommended release-hygiene artifacts for `circyto` benchmark and manuscript work.

## Recommended supplementary artifacts

- per-workdir `workflow_summary.json`
- benchmark table from `circyto summarize-benchmark`
- `anndata/circ_counts.h5ad`
- `mudata/full_length.h5mu` when available
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
- cleanup and success status
