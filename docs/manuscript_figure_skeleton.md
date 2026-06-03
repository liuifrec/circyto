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

## Figure 3. RNA+circ and scRR DNA MuData integration

Purpose:

- show the `scverse`-style integration layer from RNA+circ through processed scRR DNA modalities

Panels:

- `3A` `MuData` schema:
  - `rna`
  - `circ`
  - `cnv`
  - `rt`
  - shared `obs`
- `3B` RNA QC and circRNA QC summary fields
- `3C` GSM-to-biological-cell remapping:
  - `GSM8558852 -> RNA_IMR90_A_100 -> IMR90_A_100`
- `3D` IMR90 full23 tri-modal overlap:
  - `rna`: 23 x 63187
  - `circ`: 23 x 2443
  - `cnv`: 23 x 60607
  - trimodal overlap: 23

## Figure 4. scRR DNA branch and SComatic sidecar

Purpose:

- present processed scRR DNA state modalities and SComatic as exploratory candidate-signal interoperability

Panels:

- `4A` modality map:
  - `rna`
  - `circ`
  - `cnv`
  - `rt`
  - optional exploratory `candidate_snv`
- `4B` CNV integration contracts:
  - processed GEO `summary_CNV_states_*`
  - processed GEO `summary_CNV_mappabilitynorm_*`
  - `cnv.h5ad`
- `4C` RT/state integration contracts:
  - processed HAP1 binarized replication-state tables
  - average RT bedGraph when coordinates match
  - `rt.h5ad`
- `4D` candidate signal contract:
  - SComatic candidate summaries
  - RNA-derived candidate variant signals
  - not validated somatic mutations

## Figure 5. Public scRR biological pilot analysis

Purpose:

- provide cautious pilot biological summaries once server runs complete

Candidate panels:

- `5A` IMR90 RNA+circ+CNV tri-modal summary
- `5B` HAP1 batch10 RNA+circ and SComatic technical smoke summary
- `5C` exploratory MuData QC/UMAP views
- `5D` HAP1 full pending-run placeholder

## General note

The figure set should clearly separate:

- validated completed analyses
- pending full-data analyses
- exploratory future layers
