# RamDA/scRR RNA-Derived Candidate Variant Signal Master Plan

This plan defines the end-to-end path for integrating RamDA/scRR RNA-derived candidate variant signals into `circyto`.

Scope:

- do not implement a variant caller in `circyto`
- do not treat SComatic RNA outputs as validated somatic mutations
- keep real SComatic execution outside local WSL development
- import and integrate external candidate-variant outputs in a reproducible, scverse-compatible way

Preferred terminology:

- `RNA-derived candidate variant signals`
- not `validated somatic mutations`, unless orthogonal DNA evidence is available

## 1. Current State

### Completed Components

`circyto` currently supports:

- manifest-driven RamDA/scRR full-length circRNA workflows
- BWA-MEM single-end CIRI3 route for IMR90-like one-SRR-per-cell scRR/RamDA inputs
- STAR plus BWA rescue paired-end CIRI3 route for HAP1-like paired-end scRR/RamDA inputs
- RNA snapshot import from external gene-count tables
- lightweight RNA expression profiling from existing alignments with `simple-overlap`
- circRNA matrix collection and QC
- RNA+circ joined summaries
- `h5ad` export for circRNA matrices
- MuData export for RNA+circ objects
- processed scRR CNV import with `circyto import-scrr-cnv`
- processed scRR replication timing/state import with `circyto import-scrr-rt`
- GSM-to-biological-cell mapping with `circyto build-scrr-cell-map`
- RNA/circ MuData remapping with `circyto remap-scrr-mudata-obs`
- tri-modal RNA+circ+CNV MuData merging with `circyto merge-scrr-cnv`
- tri-modal RNA+circ+RT MuData merging with `circyto merge-scrr-rt`
- DNA/SNV scaffold tables under `WORKDIR/dna/`
- synthetic SComatic candidate-table import and RNA/circ/DNA summary joins
- RamDA/scRR to SComatic pooled-BAM adapter
- cleanup-aware and stale-workdir resume guards for full-length workflows
- benchmark/report/snapshot helpers for local and server workflow outputs

### Validated Components

Validated locally with synthetic or small chr21-style fixtures:

- full-length workflow orchestration and reporting
- RNA+circ integration bookkeeping
- MuData RNA+circ export
- DNA/SNV import schema
- synthetic and real SComatic candidate table normalization
- `import-dna-snv-summary`
- `summarize-dna-rna-circ`
- RamDA/scRR SComatic adapter behavior:
  - one-cell-per-BAM/SAM input handling
  - `CB` tag injection or normalization
  - SComatic `Index`/`Cell_type` annotation writing
  - missing alignment failure
- real SComatic Step1/Step2 output normalization
- processed GEO scRR CNV import into `cnv.h5ad`
- processed scRR replication timing/state import into `rt.h5ad` with synthetic tests
- IMR90 full23 tri-modal RNA+circ+CNV MuData

Validated on real public workflow inputs:

- Smart-seq3 original proof-of-concept benchmark path
- IMR90 scRepli-RamDA-seq two-cell pilot
- HAP1 scRepli-RamDA-seq three-cell pilot
- IMR90 full 23-cell RNA+circ workflow and processed CNV integration
- HAP1 batch10 RNA+circ workflow and SComatic technical smoke

Planned server-scale runs:

- HAP1 future full set run after remaining FASTQ download

### Exploratory / Optional Components

The following remain exploratory or optional:

- imported `scomatic_candidate_summary.tsv`
- per-cell candidate signal counts in joined summaries
- SComatic candidate generation from RamDA/scRR alignments
- `candidate_snv` MuData representation
- DNA SNV interpretation without orthogonal DNA validation

## 2. Real SComatic Execution Plan

Real SComatic execution should happen on HPC, a server, or a containerized environment. Local WSL native execution is not the target because repeated SComatic dependency installation attempts failed with `Bus error (core dumped)`.

### HPC Execution Stages

1. Complete a `circyto workflow full-length-circrna` run.
2. Confirm the workflow has not been cleanup-pruned in a way that removes required alignments.
3. Run the RamDA/scRR SComatic adapter to produce pooled BAM input.
4. Transfer or stage the adapter output into the SComatic environment.
5. Run SComatic Step 1 to split the pooled BAM by `Cell_type`.
6. Run SComatic base counting by chromosome or interval.
7. Merge base-count outputs.
8. Run SComatic candidate calling and filtering with a compatible reference and Panel of Normals.
9. Normalize final SComatic pass/fail outputs into `scomatic_candidate_summary.tsv`.
10. Import normalized candidate tables into `circyto`.
11. Generate joined RNA/circ/candidate signal summaries.
12. Export or refresh multimodal analysis objects.

### Required Inputs

From `circyto`:

- `WORKDIR/align/alignment_manifest.tsv`
- underlying SAM/BAM files referenced by the alignment manifest
- `WORKDIR/matrix/circ_counts.mtx`
- `WORKDIR/matrix/circ_index.txt`
- `WORKDIR/matrix/cell_index.txt`
- `WORKDIR/rna/gene_counts.tsv` when RNA modality is present
- `WORKDIR/rna/gene_feature_table.tsv`
- `WORKDIR/qc/*` workflow summaries and joined RNA/circ summaries
- `workflow_summary.json`

From the adapter:

- `scomatic_adapter/merged/<sample-id>.scomatic.bam`
- `scomatic_adapter/merged/<sample-id>.scomatic.bam.bai`
- `scomatic_adapter/cell_annotations.tsv`
- `scomatic_adapter/adapter_summary.json`

From SComatic environment:

- SComatic checkout or installed package
- compatible Python and R environment
- `samtools`
- `bedtools`
- reference FASTA matching the circyto alignments
- SComatic Panel of Normals appropriate for genome build and assay
- RNA-editing resources if used in filtering

### Expected Outputs

External SComatic outputs should be archived as raw provenance, then normalized into:

- `scomatic_candidate_summary.tsv`
  - one row per RNA-derived candidate variant signal or candidate signal observation
  - required by the current `circyto import-dna-snv-summary` path
- optional raw SComatic step outputs:
  - split BAM reports
  - base-cell count tables
  - merged count tables
  - SComatic calling step outputs
  - filter reports
  - command logs

Imported `circyto` outputs:

- `dna/scomatic_candidate_summary.tsv`
- `dna/dna_snv_import_summary.json`
- `qc/dna_rna_circ_cell_summary.tsv`
- `qc/dna_rna_circ_summary.json`

### Risks

- WSL local dependency instability for SComatic.
- Reference mismatch between circyto alignments, adapter BAMs, SComatic reference FASTA, and PoN files.
- `nM`/`NH` tag availability. SComatic can run without filters that require them, but using `--max_nM` or `--max_NH` requires those tags.
- RNA-derived signals are affected by expression, allele-specific expression, RNA editing, mapping artifacts, and transcript coverage.
- One-cell-per-SRR labels can create many very small SComatic groups; grouping strategy must be chosen deliberately.
- Cleanup may remove alignment files required by the adapter; use a fresh run or preserve alignments for SComatic preparation.
- SComatic output schemas may vary by version and filtering path; normalization must be version-aware.

## 3. RamDA/scRR Adapter Role

The adapter bridges circyto's one-cell-per-SRR alignment model to SComatic's pooled single-cell BAM expectation.

Adapter script:

```bash
bash scripts/prepare_ramda_scomatic_input.sh \
  --alignment-manifest WORKDIR/align/alignment_manifest.tsv \
  --outdir WORKDIR/scomatic_adapter \
  --sample-id ramda_scrr
```

Responsibilities:

- read `alignment_manifest.tsv`
- use manifest `cell_id` as the synthetic single-cell barcode
- convert SAM or BAM records into coordinate-sorted per-cell BAMs
- inject `CB:Z:<cell_id>` if `CB` is missing
- normalize non-matching existing `CB` tags to the manifest `cell_id`
- preserve `nM` and `NH` tags when present
- report missing `nM` and `NH` counts
- merge all prepared BAMs into one pooled coordinate-sorted BAM
- index the pooled BAM
- write SComatic metadata as `Index`/`Cell_type`
- write `adapter_summary.json`

Default metadata:

```text
Index	Cell_type
cellA	cellA
cellB	cellB
```

For one-cell-per-SRR RamDA/scRR, this default preserves single-cell granularity. For broader SComatic grouping, use a manifest column or `--default-cell-type`, but document the biological rationale because grouping affects SComatic filtering and interpretation.

## 4. DNA-Side Integration Plan

`circyto` should continue to treat DNA-side state as externally generated evidence and metadata.

### `dna_cell_summary.tsv`

Per-cell DNA-state table. Current contract:

- `cell_id`
- `dna_library_id`
- `cnv_burden`
- `replication_score`
- `cell_cycle_phase`
- `dna_variant_count`
- `notes`

Near-term additions can remain optional until promoted into a versioned contract:

- `aneuploidy_score`
- `replication_timing_state`
- `clone_id`
- `dna_source`
- `dna_pipeline`
- `dna_qc_status`

### `dna_variant_summary.tsv`

Externally generated DNA variant table. Current contract:

- `variant_id`
- `cell_id`
- `chrom`
- `pos`
- `ref`
- `alt`
- `gene`
- `consequence`
- `evidence_type`
- `caller`
- `filter_status`

This table is for future external DNA pipelines. It should not be populated by RNA-only SComatic results unless explicitly labeled as RNA-derived candidate evidence in separate columns or in the existing SComatic candidate table.

### `scomatic_candidate_summary.tsv`

RNA-derived candidate signal table. Current contract:

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

`candidate_variant_class` should use:

```text
RNA-derived candidate variant signal
```

### CNV Burden and Replication State

For scRR datasets, DNA replication context is biologically central. The integration layer should support:

- `replication_score`
- `cell_cycle_phase`
- `cnv_burden`
- `clone_id`
- `radiation_condition`
- `treatment`

For scRR GEO processed CNV summaries, bin-level state now enters through:

- `circyto import-scrr-cnv`
- `cnv.h5ad`
- `scrr_cnv_import_summary.json`

For HAP1 processed replication timing/state summaries, RT state now enters through:

- `circyto import-scrr-rt`
- `rt.h5ad`
- `scrr_rt_import_summary.json`

Per-cell burden or replication-state summaries may still enter through `dna_cell_summary.tsv` or external DNA/CNV/RT postprocessing when they are not available as a first-class AnnData modality.

### Future External Pipelines

Future inputs can come from:

- DNA SNV callers
- DNA CNV callers
- scRR replication timing tools
- SComatic real runs over RNA alignments
- custom server/HPC postprocessing scripts

`circyto` should normalize and join these outputs, not replace specialized callers.

## 5. Integrated Analysis Outputs

### RNA

Existing RNA modality:

- `rna/gene_counts.tsv`
- `rna/gene_feature_table.tsv`
- RNA QC fields:
  - `total_rna_counts`
  - `detected_genes`
  - mitochondrial/ribosomal fractions when available

### circRNA

Existing circRNA modality:

- `matrix/circ_counts.mtx`
- `matrix/circ_index.txt`
- `matrix/cell_index.txt`
- `matrix/circ_feature_table.tsv`
- cell-level circ fields:
  - `circRNA_count`
  - `circRNA_total_support`

### Candidate Variants

Imported RNA-derived candidate signal modality:

- `dna/scomatic_candidate_summary.tsv`
- per-cell summary:
  - `scomatic_candidate_count`
  - future optional read-support and VAF summaries

This is not a DNA modality by itself. It is currently stored under `dna/` because it belongs to the candidate-signal integration scaffold, but manuscript language should call it RNA-derived candidate variant evidence.

### DNA State

External DNA/scRR state:

- `dna/cnv.h5ad` when imported with `circyto import-scrr-cnv`
- `dna_rt/rt.h5ad` when imported with `circyto import-scrr-rt`
- `scrr_cnv_import_summary.json`
- `scrr_rt_import_summary.json`
- `dna/dna_cell_summary.tsv`
- `cnv_burden`
- `replication_score`
- `cell_cycle_phase`
- `clone_id` when available
- future DNA SNV and per-cell CNV summaries

### MuData Structure

Current RNA/circ MuData:

- `mdata["rna"]`
- `mdata["circ"]`
- shared `mdata.obs`

Current scRR tri-modal MuData after remapping and CNV merge:

- `mdata["rna"]`
  - gene counts
- `mdata["circ"]`
  - circRNA counts
- `mdata["cnv"]`
  - processed GEO CNV bin states in `X`
  - optional `layers["mappabilitynorm"]`
- `mdata.obs`
  - shared cell metadata:
    - RNA/circ membership
    - `circRNA_count`
    - `canonical_cell_id`
    - original GSM/RNA/DNA identifiers when available

Parallel HAP1 replication timing/state MuData:

- `mdata["rna"]`
- `mdata["circ"]`
- `mdata["rt"]`
  - processed replication-state or RT values in `X`
  - optional matching average RT values in `var["avg_rt"]`

Planned optional candidate-signal representation:

- `candidate_snv`
  - RNA-derived SComatic candidate signals only
  - not validated DNA somatic mutations

## 6. Manuscript Figures Enabled

This architecture enables methods-focused figures without overclaiming mutation biology.

Candidate figure panels:

- Workflow schematic:
  - FASTQ/SRR manifest to alignment to circRNA detection to RNA/circ/DNA-state integration
  - external SComatic execution shown as an optional sidecar path
- Dataset benchmark table:
  - Smart-seq3, IMR90 scRR, HAP1 scRR
  - cells, layout, aligner route, circRNA outputs, RNA snapshots, candidate-signal readiness
- Workdir integrity and reproducibility figure:
  - manifest-driven execution
  - cleanup-aware resume protection
  - snapshotable outputs
- RNA/circ integration figure:
  - per-cell RNA depth versus circRNA count
  - host-gene-linked circRNA summaries
- scRR biological context figure:
  - circRNA burden versus replication score or cell-cycle phase when DNA-state metadata is available
- Candidate signal integration figure:
  - RNA-derived candidate variant signal counts overlaid with RNA/circ and replication-state metadata
  - explicitly labeled exploratory and RNA-derived
- MuData interoperability figure:
  - `rna`, `circ`, validated `cnv`, implemented `rt`, optional exploratory `candidate_snv`
  - shared `obs` schema

These figures support a framework and benchmark manuscript. They do not require claiming validated somatic mutations.

## 7. Remaining Blockers

Execution blockers:

- real SComatic environment must be validated on HPC/server/container
- WSL native SComatic installation remains unstable
- full HAP1 set still needs remaining FASTQ download and full workflow run
- IMR90 whole-genome SComatic candidate calling remains exploratory/in progress

Data blockers:

- HAP1 full public dataset FASTQs/alignments must complete
- reference and PoN compatibility must be documented per run
- cell grouping strategy for SComatic must be chosen:
  - one cell per `Cell_type`
  - grouped by sample condition
  - grouped by inferred state
- candidate output normalization must be finalized against real SComatic output files

Interpretation blockers:

- RNA-derived signals require RNA-editing and mapping-artifact caution
- DNA validation is absent unless supplied by external DNA pipelines
- candidate signal burden may correlate with expression/coverage rather than mutation burden
- HAP1/IMR90 biological conclusions should remain exploratory until validated

Implementation blockers:

- no stable `candidate_snv` MuData modality export yet
- no formal provenance schema for external SComatic command logs and resource versions
- no server-side run template for SComatic over adapter output yet

## Recommended Next Milestone

Keep using the read-only SComatic result normalizer:

```bash
circyto normalize-scomatic-results \
  --scomatic-output SCOMATIC_STEP_OUTPUT.tsv \
  --cell-annotations WORKDIR/scomatic_adapter/cell_annotations.tsv \
  --outdir WORKDIR/scomatic_normalized
```

Requirements:

- parse real SComatic Step1/Step2 or final pass/fail tables
- preserve raw SComatic columns in provenance or sidecar files
- emit the existing `scomatic_candidate_summary.tsv` contract
- include resource provenance:
  - SComatic commit/version
  - reference FASTA
  - PoN file
  - RNA-editing resource
  - command line
- never run SComatic internally

This milestone converts external candidate generation into a reproducible `circyto` import boundary while keeping variant calling outside the package.
