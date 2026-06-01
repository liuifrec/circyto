# Manuscript Figure Skeleton

This page is a planning scaffold for a methods-oriented `circyto` manuscript.

## Figure 1. circyto workflow architecture

Purpose:

- show the workflow hierarchy from FASTQ/BAM to circRNA matrices, RNA profiles, MuData, and QC summaries

Panels:

- `1A` command surface:
  - `workflow smartseq3-ciri3`
  - `workflow full-length-circrna`
  - post-hoc profiling and summary commands
- `1B` single-end RamDA/scRR route:
  - `BWA-MEM -> direct SAM -> CIRI3`
- `1C` paired-end RamDA/scRR route:
  - `STAR -> BWA rescue -> CIRI3`
- `1D` output objects:
  - matrix
  - `h5ad`
  - `h5mu`
  - QC/provenance

## Figure 2. Public dataset reconstruction and benchmark table

Purpose:

- show dataset structure and benchmarked public runs without overclaiming biology

Panels:

- `2A` dataset structure:
  - Smart-seq3 pooled/demux-required
  - scRR/RamDA one-library-per-run
- `2B` benchmark summary table:
  - all192
  - IMR90 pilot / full set
  - HAP1 pilot / batch set
- `2C` workflow reproducibility and integrity outputs

## Figure 3. RNA+circ MuData integration

Purpose:

- show the `scverse`-style integration layer

Panels:

- `3A` `MuData` schema:
  - `rna`
  - `circ`
  - shared `obs`
- `3B` RNA QC and circRNA QC summary fields
- `3C` RNA-only / circ-only / shared cell bookkeeping

## Figure 4. scRR multimodal roadmap: DNA / RNA / circ / SComatic

Purpose:

- present the roadmap without claiming full DNA/SNV integration is already validated

Panels:

- `4A` future modality map:
  - `rna`
  - `circ`
  - `dna_cnv`
  - `dna_snv`
- `4B` candidate integration contracts:
  - DNA cell summaries
  - DNA variant summaries
  - SComatic candidate summaries
- `4C` terminology guardrail:
  - RNA-derived candidate variant signals
  - not validated somatic mutations

## Figure 5. Public scRR biological pilot analysis

Purpose:

- provide cautious pilot biological summaries once server runs complete

Candidate panels:

- `5A` IMR90 RNA+circ summary
- `5B` HAP1 RNA+circ summary
- `5C` exploratory MuData QC/UMAP views
- `5D` future bridge to DNA/RNA/circ integration

## General note

The figure set should clearly separate:

- validated completed analyses
- running server analyses
- exploratory future layers
