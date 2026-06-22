# Figures Plan

## Figure 1: circyto Architecture

Planned panels:

- FASTQ/full-length scRNA input.
- circRNA detector layer.
- host-gene annotation.
- known circRNA annotation.
- AnnData/MuData export.
- integration with RNA, CNV, RT, and `candidate_snv` RNA-derived candidate variant signal layers.

## Figure 2: Dataset Validation and Host-Gene Annotation

Planned panels:

- Smart-seq3, HAP1, and IMR90 inventory.
- circRNA counts by dataset.
- host-gene annotation recovery.
- `host_gene_source` breakdown.

## Figure 3: HAP1 RT-circRNA Association

Planned panels:

- RT-positive fraction versus circRNA count.
- RT-positive fraction versus circRNA total support.
- OLS / partial regression result for `circRNA_count ~ detected_genes + frac_rt_pos`.
- RT-high enriched host genes.

## Figure 4: IMR90 CNV-circRNA Analysis

Planned panels:

- global CNV burden versus circRNA count/support.
- CNV-high versus CNV-low host-gene enrichment.
- ECM/stress-associated host genes, including `COL1A1`, `FN1`, `VIM`, `COL6A3`, `P4HB`, `COL1A2`, `CRIM1`, `SPARC`, `LTBP1`, `CORO1C`, and `ERBIN`.
- local CNV-at-circRNA-locus results when coordinates support the analysis.

## Figure 5: Cross-Dataset Host-Gene Program

Planned panels:

- pairwise and three-way host-gene overlaps.
- shared HAP1 RT-high and IMR90 CNV-high positive host genes.
- candidate recurrent genome-state-associated circRNA host-gene program.
