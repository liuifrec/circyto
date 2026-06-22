# circyto Manuscript Workspace

Working title:

> Genome-state-associated circular RNA programs revealed by multimodal full-length single-cell sequencing with circyto

This directory organizes manuscript-facing plans and reproducibility notes for `circyto` `v0.10.0`.

`circyto` should be described as a CLI/scverse-compatible framework for single-cell circular RNA detection, annotation, and multimodal integration from full-length single-cell RNA-seq and full-length single-cell multi-omic data. The framework connects RNA, circRNA, replication timing, CNV, and exploratory RNA-derived candidate variant signal layers in AnnData and MuData objects.

## Contents

- `analysis_plan.md`: planned analyses and decision points.
- `figures_plan.md`: proposed figure structure.
- `tables_plan.md`: proposed manuscript and supplement tables.
- `methods_commands.md`: command shapes for reproducing processed objects and summaries.
- `data_inventory.md`: expected processed datasets and current manuscript numbers.
- `results_checklist.md`: results needed before submission.
- `caveats.md`: conservative interpretation language.
- `rerf_radiation_positioning.md`: internal RERF/RP positioning and limits.

Companion reusable scripts live under `scripts/manuscript/`. Example command files live under `examples/manuscript/`.

## Reproducibility Boundary

The repository does not include the large FASTQ, reference genome, GEO processed tables, circAtlas tables, or completed `.h5mu` objects used for manuscript-scale analyses. Documentation uses placeholder paths such as `/path/to/data/...` for private or local resources.

Committed materials should be sufficient to understand the workflow, rerun scripts on equivalent local inputs, and test script behavior on synthetic objects.
