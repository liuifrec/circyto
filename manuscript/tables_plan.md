# Tables Plan

## Main Tables

Table 1: Dataset and modality inventory.

- dataset
- source accession
- protocol
- RNA cells/features
- circRNA cells/features
- DNA-derived modality cells/features
- shared-cell overlap
- processed object path

Table 2: Host-gene annotation recovery.

- dataset
- circRNA count
- annotated circRNAs
- recovery fraction
- `host_gene_source` counts
- dominant source

Table 3: HAP1 RT-circRNA association.

- metric pair
- number of paired cells
- Pearson r and p-value
- Spearman r and p-value
- OLS terms for `circRNA_count ~ detected_genes + frac_rt_pos`

Table 4: IMR90 CNV-high host-gene programs.

- host gene
- number of circRNAs
- CNV-high mean support
- CNV-low mean support
- difference and log2 fold-change
- gene-program label when curated

## Supplement Tables

Supplement Table S1: Full `summarize_mudata_inventory.py` output.

Supplement Table S2: Known versus novel circRNA summary.

Supplement Table S3: Cross-dataset host-gene overlap.

Supplement Table S4: Optional local CNV-at-circRNA-locus summary.

Supplement Table S5: Reproducibility commands and software versions.
