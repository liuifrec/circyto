# scRR Multimodal Roadmap

This note sketches the next `circyto` architecture step for scRR-seq-style multimodal outputs.

Current validated foundation:

- `rna`: post-hoc or workflow-time simple-overlap RNA profiling
- `circ`: circRNA matrix, QC, `h5ad`, and `MuData` export
- `cnv`: processed scRR GEO CNV state import and IMR90 full23 tri-modal MuData
- `rt`: processed scRR replication timing/state import for HAP1-style DNA tables, with synthetic test coverage

Current target:

- integrate RNA, circRNA, processed scRR DNA CNV or RT/state profiles, replication-state metadata when available, and exploratory candidate variant summaries under one shared cell axis

## Proposed modalities

- `mdata["rna"]`
  - gene counts
  - future velocity-compatible layers
- `mdata["circ"]`
  - circRNA counts
  - circ feature coordinates and host-gene annotations
- `mdata["cnv"]`
  - processed scRR GEO copy-number profiles
  - bin-level CNV features
- `mdata["rt"]`
  - processed scRR replication timing/state profiles
  - genomic or gene-intersect RT features
- `mdata["candidate_snv"]`
  - future exploratory candidate variant modality
  - likely sparse cell-by-variant or summarized burden layers

## Shared `obs` metadata targets

The long-term shared cell metadata layer should support:

- `replication_score`
- `cell_cycle_phase`
- `cnv_burden`
- `snv_burden`
- `scomatic_pass_variants`
- `radiation_condition`
- `treatment`
- `clone_id`

These fields should stay optional and provenance-aware.

## Current architectural stance

`circyto` should:

- preserve one shared cell axis across modalities
- export clean, inspectable `MuData`
- keep modality-specific feature tables explicit
- keep RNA/circ summaries lightweight and reproducible
- add RNA-derived candidate variant signals modalities only as explicit layers, not hidden workflow side effects

`circyto` should not yet:

- run full DNA calling internally
- claim orthogonally confirmed somatic variant inference from RNA alone
- collapse DNA, RNA, and circ evidence into a single overinterpreted score

## Near-term helper layer

To prepare for that multimodal future, `circyto` now adds lightweight workdir helpers:

- `circyto inspect-workdir`
  - summarize available modalities and workflow artifacts in a completed directory
- `circyto summarize-circ-host-genes`
  - summarize host-gene recurrence across detected circRNAs
- `circyto export-circ-bed`
  - export circ coordinates and support for genome-browser interoperability

These helpers are intentionally lightweight and do not imply biological claims by themselves.

## Integration path

1. Keep `rna` and `circ` outputs stable.
2. Import processed CNV summaries with `circyto import-scrr-cnv` where source files are CNV.
3. Import processed replication timing/state summaries with `circyto import-scrr-rt` where source files are RT/state.
4. Add optional replication-state metadata into shared `obs` when derived summaries are available.
5. Remap GSM IDs and merge tri-modal MuData with `remap-scrr-mudata-obs` plus `merge-scrr-cnv` or `merge-scrr-rt`.
6. Add RNA-derived candidate variant signals burden summaries as optional exploratory `candidate_snv`.
