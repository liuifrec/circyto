# Release Notes

## Unreleased

### Full-Length RNA SComatic Adapter

- Added `circyto prepare-scomatic-input` for Smart-seq2, Smart-seq3, RamDA-seq, ShinRamDA, and scRR RNA-arm SComatic input preparation.
- Added `circyto run-scomatic` to write an external BaseCellCounter/Step1/Step2 command plan and optionally execute it with explicit `--execute`.
- Added `circyto import-scomatic` as a conservative import boundary for SComatic Step1/Step2 or candidate tables.
- Added `circyto merge-scomatic` to export normalized RNA-derived candidate variant signals as an exploratory `candidate_snv` MuData modality.
- SComatic terminology remains `RNA-derived candidate variant signals`, not validated somatic mutations.
- `run-scomatic` defaults to planning-only; use split environments and prefer server/container/HPC execution because local WSL SComatic execution may be unstable.
- Validation status: HAP1 batch10 technical smoke and real Step1/Step2 normalization are validated; `candidate_snv` MuData export is partially validated with synthetic tests; clinically validated variants and validated somatic mutation discovery are not yet validated.
- Lessons learned: HAP1 confirms the technical adapter path but also shows homogeneous cell-type designs may yield no Step2 candidates; IMR90 reinforces that processed DNA state, such as CNV, is the primary scRR DNA evidence while SComatic remains an optional RNA-derived candidate variant signals sidecar.

### HAP1 scRR Replication Timing/State Import

- Added `circyto import-scrr-rt` for processed scRR DNA replication timing/state tables.
- Added `circyto merge-scrr-rt` for RNA+circ+RT MuData construction.
- The `rt` modality is separate from IMR90 `cnv`; HAP1 binarized replication-state and average RT files should not be treated as CNV by default.
- `rt.h5ad` uses cells as `obs`, processed genomic or gene-intersect features as `var`, and replication-state/RT values in `X`.
- Optional average RT bedGraph values are stored as `var["avg_rt"]` only when coordinates match the RT table features.
- Synthetic tests cover binary RT parsing, canonical HAP1 cell IDs, no-h5ad mode, optional h5ad writing, and RT MuData merge when `mudata` is available.

Known limitation: the named GSE278952 HAP1 processed files were not present locally during implementation, so real-file import validation remains pending. The importer uses processed GEO-style tables and does not rerun raw DNA FASTQ/scRepli-seq processing.

## v0.10.0

Recommended tag: `v0.10.0`

This release is a minor milestone because it extends circyto from validated RNA+circ workflows into validated scRR RNA+circ+CNV integration using processed public scRR DNA outputs.

### RNA/circ Workflow Validation

- IMR90 full 23-cell scRR-seq RNA/circ workflow completed.
- HAP1 10-cell paired-end scRR-seq batch workflow completed.
- RNA + circ MuData export works for completed full-length workflows.
- Cleanup-aware resume and stale-workdir safeguards remain part of the validated workflow hygiene layer.

### scRR CNV Import

- `circyto import-scrr-cnv` imports processed GEO CNV state tables.
- `X` stores integer CNV state parsed from labels such as `2-somy`.
- `layers["mappabilitynorm"]` stores the optional mappability-normalized CNV signal.
- IMR90 processed DNA CNV data were imported for 23 cells at 50 kb resolution.
- The CNV import summary is written as `scrr_cnv_import_summary.json`.

### GSM-to-Biological-Cell Mapping

- `circyto build-scrr-cell-map` parses GEO `family.soft` or `family.soft.gz` sample blocks.
- The mapping bridges:

```text
GSM8558852 -> RNA_IMR90_A_100 -> DNA_IMR90_A_100 -> IMR90_A_100
```

- Output columns are `gsm_id`, `rna_cell_id`, `dna_cell_id`, `canonical_cell_id`, `sample_title`, `molecule`, `treatment`, and `source_name`.

### Tri-Modal RNA+circ+CNV MuData

- `circyto remap-scrr-mudata-obs` remaps RNA/circ MuData obs IDs from GEO GSM IDs to canonical scRR biological-cell IDs.
- `circyto merge-scrr-cnv` merges remapped RNA/circ MuData with CNV AnnData.
- IMR90 full23 tri-modal MuData has been validated:
  - `rna`: 23 x 63187
  - `circ`: 23 x 2443
  - `cnv`: 23 x 60607
  - trimodal overlap: 23

### SComatic Interoperability

- RamDA/scRR to SComatic adapter behavior is validated.
- CB tags are injected or normalized.
- NH/nM tags are preserved when present.
- HAP1 batch10 BaseCellCounter, Step1, and Step2 execution has been validated as a technical smoke.
- `circyto normalize-scomatic-results` supports real SComatic Step1/Step2 outputs.
- SComatic results remain `RNA-derived candidate variant signals`, not validated somatic mutations.

### Cleanup/Resume Safeguards

- Cleanup-aware stale-workdir protection is part of the validated workflow hygiene layer.
- `cleanup-workflow` reduced the HAP1 batch10 working directory from about 132G to about 9.2G while preserving final matrices, AnnData/MuData artifacts, QC, provenance, and SComatic-ready artifacts.
- Alignment cleanup must preserve final matrix, `h5ad`/`h5mu`, QC, provenance, and SComatic adapter artifacts.

### Known Limitations

- HAP1 full set is still pending full FASTQ download and full workflow run.
- IMR90 SComatic whole-genome candidate calling is still exploratory/in progress.
- `candidate_snv` is RNA-derived and should not be described as validated somatic mutation.
- CNV import currently uses processed GEO summary tables, not raw DNA FASTQ reprocessing.
- MuData may emit an upstream `FutureWarning` from `mudata` update behavior; this is not currently a circyto failure.
- CNV segmentation and raw scRepliseq DNA reprocessing are not implemented inside circyto.
