# Bioinformatics Application Note figure skeleton

The main paper has one figure. This scaffold supersedes the former five-figure
genome-state biology design.

## Figure 1. circyto workflow and Smart-seq3 outputs

### 1A. Architecture

```text
full-length scRNA-seq
    -> protocol-aware alignment / circRNA detection
    -> cell x circRNA matrix + QC
    -> host-gene / known-circRNA annotation
    -> AnnData / MuData
```

Design rules:

- distinguish upstream detectors from circyto's orchestration, collection,
  annotation, QC, provenance, and object layers;
- emphasize RNA and circRNA as the core modalities;
- show processed CNV, processed RT, and exploratory candidate variant signals
  as small optional branches;
- do not imply that every optional branch is biologically validated.

### 1B. Smart-seq3 circRNA burden

- E-MTAB-8735, 192 cells;
- RNA-derived UMAP with fixed preprocessing and seed;
- per-cell detected circRNA candidates as an overlay;
- median 12 detected circRNAs/cell may be stated only after regeneration from
  the checksum-matched manuscript object.

### 1C. Representative MAN1A2-associated candidate

- same RNA-derived coordinates as 1B;
- binary detection overlay for `chr1:117402186|117420649`;
- annotate the feature as MAN1A2-associated;
- historical values (6/192 cells, total support 15) require regeneration;
- no functional or mechanistic interpretation.

## Supplementary Figure S1. Protocol/workflow validation

Compare the completed Smart-seq3, IMR90 single-end scRR, and HAP1 paired-end
scRR routes. Report inputs, alignment/detector path, outputs, validation scale,
and limitations. Use the validated HAP1 10-cell RNA/circ route.

## Supplementary Figure S2. Multimodal object/schema

Show RNA+circ outputs, validated IMR90 RNA+circ+processed-CNV integration with
23 shared cells, and the optional RT import/merge contract. Label the RT
contract as synthetic-test validated with real processed-file validation
pending reconciliation.

## Deferred biological follow-up

HAP1 RT regression, IMR90 CNV associations, local CNV-at-circRNA analyses,
cross-dataset host-gene programs, and radiation applications are not current
Application Note figures.
