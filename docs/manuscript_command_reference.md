# Manuscript Command Reference

This page collects the command families most relevant to manuscript reproducibility. See `docs/manuscript_reproducibility.md` for full command shapes.

## Full-Length Workflows

- `circyto workflow smartseq3-ciri3`: pooled Smart-seq3 FASTQs to circRNA outputs, optional RNA and MuData export.
- `circyto workflow full-length-circrna`: manifest-driven RamDA, Shin-RamDA, or scRR full-length RNA circRNA workflow.

## Annotation and Repair

- `circyto annotate-circs`: known circRNA annotation against circAtlas-style or compatible tables.
- `circyto repair-host-genes`: repair or add host-gene provenance fields to `.h5ad` or `.h5mu` objects.

## scRR Genome-State Integration

- `circyto build-scrr-cell-map`: build RNA/DNA cell mapping from GEO SOFT metadata.
- `circyto remap-scrr-mudata-obs`: remap RNA/circRNA MuData observations to canonical biological-cell IDs.
- `circyto import-scrr-cnv`: import processed scRR CNV summaries as a CNV AnnData modality.
- `circyto merge-scrr-cnv`: merge RNA+circRNA MuData with a CNV modality.
- `circyto import-scrr-rt`: import processed scRR replication timing/state summaries as an RT AnnData modality.
- `circyto merge-scrr-rt`: merge RNA+circRNA MuData with an RT modality.

## Candidate Variant Signal Interoperability

- `circyto prepare-scomatic-input`: prepare full-length RNA alignments and cell annotations for external SComatic use.
- `circyto run-scomatic`: write and optionally execute external SComatic command plans.
- `circyto import-scomatic`: normalize SComatic candidate tables.
- `circyto merge-scomatic`: merge exploratory RNA-derived candidate variant signals as `candidate_snv`.

## Manuscript Scripts

- `scripts/manuscript/summarize_mudata_inventory.py`
- `scripts/manuscript/hap1_rt_circ_analysis.py`
- `scripts/manuscript/imr90_cnv_circ_analysis.py`
- `scripts/manuscript/cross_dataset_host_overlap.py`
- `scripts/manuscript/known_novel_circ_summary.py`
