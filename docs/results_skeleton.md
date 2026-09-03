# Bioinformatics Application Note results skeleton

This is a writing scaffold, not an evidence source. Resolve every quantitative
statement through `manuscript/application_note_evidence.md`.

## 1. Introduction

- Full-length single-cell protocols can retain RNA species that are not
  represented by standard gene-expression workflows.
- Existing circRNA detectors produce detector-specific outputs but do not by
  themselves supply a common cell-by-feature/scverse analysis layer.
- circyto fills that workflow and interoperability gap; it is not a new
  back-splice detection algorithm.

## 2. circyto workflow

- Describe protocol-aware Smart-seq3 and manifest-driven RamDA/scRR routes.
- Follow detector output through cell-level collection, matrix construction,
  QC, host-gene/known-circRNA annotation, provenance, and AnnData/MuData export.
- Introduce Figure 1A and keep optional processed CNV, processed RT, and
  candidate-variant integrations visually and verbally secondary.

## 3. Single-cell and multimodal outputs

- Use E-MTAB-8735 Smart-seq3 as the primary real-data demonstration.
- Report the regenerated object shape, annotation recovery, and per-cell
  summaries.
- Describe Figure 1B as an RNA-derived UMAP with circRNA burden overlaid.
- Describe Figure 1C as an illustrative MAN1A2-associated candidate detection
  overlay without functional interpretation.
- Briefly state that the 23-cell IMR90 RNA+circ+processed-CNV object
  demonstrates multimodal interoperability.
- Point workflow breadth, object schemas, software validation, commands, and
  provenance to Supplementary Figures S1-S2 and Tables S1-S3.

## 4. Conclusions

- Reiterate the bridge from detector calls to interoperable single-cell
  objects.
- Keep biological discovery, RT/CNV association analysis, long-read validation,
  and predictive biogenesis outside the central claim.
- State relevant limitations: detector-derived candidates, full-length protocol
  focus, small multimodal example, and optional integration boundaries.

## Required wording guardrails

- “circRNA candidate” or “detected circRNA,” unless orthogonally validated.
- “processed CNV integration,” not CNV inference from raw DNA reads.
- “RNA-derived candidate variant signals,” not validated DNA variants.
- No finished HAP1 RT-circRNA discovery, radiation result, single-cell Nanopore
  circRNA validation, or predictive biogenesis claim.
