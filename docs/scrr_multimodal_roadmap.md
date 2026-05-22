# scRR Multimodal Roadmap

This note sketches the next `circyto` architecture step for scRR-seq-style multimodal outputs.

Current validated foundation:

- `rna`: post-hoc or workflow-time simple-overlap RNA profiling
- `circ`: circRNA matrix, QC, `h5ad`, and `MuData` export

Future target:

- integrate RNA, circRNA, replication-state metadata, DNA copy-number summaries, and exploratory candidate variant summaries under one shared cell axis

## Proposed future modalities

- `mdata["rna"]`
  - gene counts
  - future velocity-compatible layers
- `mdata["circ"]`
  - circRNA counts
  - circ feature coordinates and host-gene annotations
- `mdata["dna_cnv"]`
  - future scDNA or inferred copy-number profiles
  - segment/bin-level or gene-level CNV features
- `mdata["dna_snv"]`
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
- add DNA/SNV modalities only as explicit future layers, not hidden workflow side effects

`circyto` should not yet:

- run full DNA calling internally
- claim validated somatic mutation inference from RNA alone
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

## Future integration path

1. Keep `rna` and `circ` outputs stable.
2. Add optional replication-state metadata imports into shared `obs`.
3. Add future CNV burden summaries as `dna_cnv`.
4. Add future candidate SNV burden summaries as `dna_snv`.
5. Export richer `MuData` while preserving current circ-focused outputs.
