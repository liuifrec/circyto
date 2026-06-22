# Analysis Plan

## Core Framing

The manuscript evaluates `circyto` as a framework for full-length single-cell circRNA detection, annotation, and multimodal genome-state integration. It is not framed as only a circRNA caller.

Main analysis layers:

- Manifest-driven full-length single-cell workflow.
- Multi-detector circRNA architecture, currently validated mainly with CIRI3.
- AnnData and MuData export for scverse-compatible downstream analysis.
- Host-gene annotation with explicit provenance.
- Known circRNA annotation through circAtlas-style tables.
- scRR GEO SOFT RNA/DNA cell mapping.
- scRR CNV and replication timing modality import and merge.
- SComatic interoperability for RNA-derived candidate variant signals.

## Dataset Analyses

### Smart-seq3 / E-MTAB-8735

Use `mtab8735_smartseq3/full_length.hostgene_fixed.h5mu` for general circRNA detection, annotation benchmarking, and host-gene landscape.

Expected object:

- RNA: 192 cells x 63,187 genes
- circRNA: 192 cells x 2,503 circRNAs
- host-gene annotation: 2,379 / 2,503 circRNAs, 95.0% recovery
- `host_gene_source` mostly `gtf`

### HAP1 scRR-seq

Use `scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu` for RNA+circRNA+RT integration.

Expected object:

- RNA: 63 cells x 63,187 genes
- circRNA: 63 cells x 3,209 circRNAs
- RT: 56 cells x 56,881 genomic bins
- RNA/circRNA overlap: 63 cells
- RNA/circRNA/RT overlap: 56 cells
- host-gene annotation: 3,117 / 3,209 circRNAs, 97.1% recovery

Primary test:

```text
circRNA_count ~ detected_genes + frac_rt_pos
```

This tests whether replication-timing state explains circRNA burden beyond detected gene number / transcriptional complexity.

### IMR90 scRR-seq

Use `scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu` for RNA+circRNA+CNV integration.

Expected object:

- RNA: 23 cells x 63,187 genes
- circRNA: 23 cells x 2,443 circRNAs
- CNV: 23 cells x 60,607 genomic bins
- RNA/circRNA/CNV overlap: 23 cells
- host-gene annotation: 2,429 / 2,443 circRNAs, 99.4% recovery

Primary interpretation:

IMR90 does not show a strong global CNV-burden to circRNA-burden relationship, but CNV-high cells show biologically coherent fibroblast, ECM, and stress-associated host-gene programs.

Local CNV analysis should test whether circRNA abundance follows copy number at the circRNA locus for host genes such as `COL1A1`, `FN1`, `VIM`, `COL6A3`, and `P4HB`.

## Cross-Dataset Analyses

- Pairwise and three-way circRNA host-gene overlap across Smart-seq3, HAP1, and IMR90.
- Shared positive host-gene programs between HAP1 RT-high and IMR90 CNV-high analyses when enrichment tables are available.
- Known versus novel circRNA summaries if circAtlas annotation fields exist in the processed objects.

## Interpretation Rules

- HAP1 and IMR90 public scRR analyses are proof-of-concept genome-state integrations, not exposure-specific conclusions.
- The current manuscript should emphasize framework plus proof-of-concept discovery, not causality.
- SComatic-derived layers must be described as RNA-derived candidate variant signals.
