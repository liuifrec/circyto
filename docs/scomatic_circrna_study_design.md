# SComatic-Integrated circRNA Study Design

This document defines a future study design for linking `circyto` full-length single-cell circRNA outputs with analysis of external SComatic-based RNA-derived candidate variant signals.

It is a design and terminology document, not a workflow execution guide.

## A. Study rationale

`circyto` now supports biologically relevant full-length single-cell RNA workflows across:

- pooled Smart-seq3
- single-end human scRR / RamDA-like data
- paired-end human scRR data

This creates a realistic foundation for an integrated study of:

- circRNA burden
- full-length single-cell transcriptome state
- RNA-derived candidate variant signals
- replication stress and genome instability biology
- future radiation-response and genome-instability themes relevant to RERF

The scientific aim is hypothesis generation, prioritization, and structured follow-up.

The scientific aim is not direct proof of DNA variant causality from scRNA-seq alone.

## B. Dataset tiers

### Tier 1: validated public human full-length datasets

- `GSE278958`
  - IMR90
  - single-end scRR
  - validated `BWA+CIRI3` route
- `GSE278952`
  - HAP1
  - paired-end scRR
  - validated `STAR+CIRI3` route with explicit paired-end opt-in

These are the first public human datasets suitable for circRNA plus candidate variant signal pilots.

### Tier 2: Smart-seq3 benchmark

- `E-MTAB-8735`
  - pooled Smart-seq3 benchmark
  - validated end-to-end in `circyto`

This tier is useful for workflow benchmarking, schema stabilization, and exploratory cross-protocol comparisons.

### Tier 3: future RERF-aligned internal datasets

Examples if scientifically and politically appropriate:

- radiation-response cohorts
- clonal hematopoiesis-linked full-length single-cell datasets
- RP1-23-1-like internal projects

These should be treated as future internal pilots only after public-dataset contracts are stable.

## C. Required inputs

Minimum study inputs:

- raw FASTQ or validated alignment inputs for `circyto`
- reference FASTA and GTF used for circRNA calling
- per-cell circRNA outputs from `circyto`
- circRNA matrix outputs
- `h5ad` outputs with cell-level QC
- per-cell BAM manifest or equivalent alignment manifest for external SComatic preparation
- cell metadata table with experimental condition, cell type, and replication-state labels if available
- external SComatic RNA-derived candidate variant signals table

Optional but important inputs:

- host-gene annotations for circRNAs
- replication-state assignments for scRR datasets
- gene sets related to DNA damage response, replication stress, and radiation biology
- orthogonal DNA validation data if available

## D. Workflow architecture

Proposed study architecture:

```text
FASTQ/BAM
-> circyto circRNA detection
-> circRNA x cell matrix + h5ad
-> export clean SComatic-compatible input tables
-> run SComatic externally under explicit user control
-> import RNA-derived candidate variant signals table
-> per-cell burden of RNA-derived candidate variant signals
-> per-cell circRNA burden
-> host-gene circRNA features
-> circRNA / candidate variant association summaries
```

Role boundaries:

`circyto` should:

- export clean SComatic inputs
- collect SComatic outputs
- join RNA-derived candidate variant signals summaries with circRNA matrices
- produce interpretable summary tables
- preserve provenance

`circyto` should not yet:

- claim validated somatic mutation calling
- run full SComatic automatically without explicit user control
- mix DNA and RNA evidence without explicit labels
- overfit biological conclusions from 1-2 cell pilots

## E. Output tables

Proposed future study outputs:

- `circ_snv_cell_summary.tsv`
  - one row per cell
  - circRNA burden
  - burden of RNA-derived candidate variant signals
  - QC and metadata joins
- `circ_snv_host_gene_summary.tsv`
  - host-gene aggregation across circRNA features and candidate variant signals
- `circ_snv_candidate_variant_summary.tsv`
  - RNA-derived candidate variant signals summary table
  - recurrence, gene context, cell distribution, and filtering provenance
- `circ_snv_recurrent_circRNA_summary.tsv`
  - recurrent circRNAs stratified by burden of RNA-derived candidate variant signals or condition
- `circ_snv_provenance.json`
  - references, command shapes, input paths, versions, terminology flags, and join assumptions

Current scaffold status:

- already present:
  - `circ_snv_cell_summary.tsv`
  - `circ_snv_host_gene_summary.tsv`
- not yet implemented:
  - `circ_snv_candidate_variant_summary.tsv`
  - `circ_snv_recurrent_circRNA_summary.tsv`
  - `circ_snv_provenance.json`

## F. Statistical analyses

Recommended first-pass analyses:

- per-cell circRNA count vs burden of RNA-derived candidate variant signals
- per-cell circRNA read-support burden vs burden of RNA-derived candidate variant signals
- host-gene overlap enrichment between circRNA host genes and candidate variant genes
- recurrence of circRNAs in high-candidate-variant-burden cells
- replication-state-aware comparisons for scRR datasets
- annotation of candidate-variant-associated genes against genome instability, replication stress, and radiation-response gene sets

Recommended statistical style:

- descriptive summaries first
- robust effect-size estimates
- nonparametric comparisons where cell counts are small
- permutation or bootstrap strategies only after sample sizes justify them
- no mechanistic or causal language without validation

## G. Biological interpretation limits

Use conservative terms:

- “RNA-derived candidate variant signals”
- “candidate variant signals”
- “exploratory RNA-derived candidate variant signals”

Do not use:

- “validated somatic mutations”
- “DNA mutations”
- “mutation burden”

unless orthogonal DNA validation exists and the tables are explicitly labeled as such.

Key limitations:

- transcriptional noise
- RNA editing
- allele-specific expression
- mapping artifacts
- protocol-specific bias
- low cell counts
- mismatch between RNA signal and DNA state

## H. Validation Status

Validated:

- public full-length RNA/circ workflow routes for IMR90 and HAP1
- HAP1 batch10 SComatic technical smoke through BaseCellCounter, Step1, Step2, and normalization
- terminology and schema guardrails for RNA-derived candidate variant signals

Partially validated:

- candidate table joins with circRNA summaries
- `candidate_snv` MuData representation for RNA-derived candidate variant signals
- split-environment preparation/import workflow

Not yet validated:

- clinically validated variants
- validated somatic mutation discovery
- biological interpretation of genome-scale SComatic outputs
- stable candidate behavior in homogeneous datasets or datasets without meaningful cell-type/group contrasts

## I. Lessons Learned from HAP1 and IMR90

- HAP1 batch10 shows that the technical path is feasible but that homogeneous grouping can yield no Step2 candidates.
- IMR90 shows that processed DNA state integration, especially CNV, is the stronger validated scRR DNA branch.
- Both datasets require explicit separation between DNA modalities and RNA-derived candidate variant signals.

## J. Validation strategy

Validation should proceed in layers:

1. terminology and schema validation
2. public-data contract validation
3. exploratory biological replication across public datasets
4. orthogonal validation where possible

Concrete validation steps:

- verify export/import schema stability on Tier 1 public datasets
- confirm that cell identifiers remain consistent across circRNA and SComatic tables
- benchmark whether RNA-derived candidate variant signals summaries are stable under filtering perturbations
- cross-check recurrent genes against known replication-stress / DNA damage response biology
- require orthogonal DNA or targeted validation before any true somatic-mutation claim

## K. RERF / internal-politics positioning

Suggested positioning:

- emphasize method development for integrated single-cell transcriptome plus circRNA plus RNA-derived candidate variant signals analysis
- emphasize hypothesis generation for genome instability and radiation biology
- avoid claiming validated somatic mutation discovery from scRNA-seq alone
- position public human scRR datasets as the external methodological bridge before internal dataset use
- frame internal pilots as carefully labeled exploratory studies aligned with long-term mutation and transcriptome themes

This positioning is safer scientifically and institutionally than claiming validated somatic mutation discovery from RNA-only evidence.

## L. Implementation roadmap

### Phase 0: design and terminology

- freeze conservative language
- document scope limits
- define intended output schema names

### Phase 1: input/output contracts only

- stabilize BAM manifest export
- stabilize cell metadata requirements
- document expected RNA-derived candidate variant signals input schema

### Phase 2: SComatic export/import wrappers

- improve `export-scomatic-inputs`
- improve `join-circ-snv-summary`
- add provenance capture

### Phase 3: RNA-derived candidate variant signals + circRNA summary joins

- add RNA-derived candidate variant signals summary table
- add recurrent circRNA summary table
- add replication-state-aware stratification hooks

### Phase 4: public scRR pilot

- apply the scaffold to:
  - `GSE278958`
  - `GSE278952`
- evaluate descriptive signal patterns only

### Phase 5: RERF-aligned internal pilot

- only after public contracts are stable
- only with explicit exploratory labeling
- only with careful review of politics, wording, and validation boundaries

## Immediate conclusion

The current `circyto` SComatic scaffold is sufficient for exploratory export and basic circRNA/RNA-derived candidate variant signals summary joins.

It is not yet sufficient for a full integrated study layer without:

- richer summary tables
- provenance tracking
- clearer replication-state-aware metadata contracts
- stronger safeguards around interpretation language
