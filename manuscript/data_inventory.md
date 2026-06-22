# Data Inventory

The processed objects below are expected local manuscript inputs. They are not committed to the repository.

## Smart-seq3 / E-MTAB-8735

Expected object:

```text
mtab8735_smartseq3/full_length.hostgene_fixed.h5mu
```

Expected modalities:

- RNA: 192 cells x 63,187 genes
- circRNA: 192 cells x 2,503 circRNAs

Expected host-gene annotation:

- 2,379 / 2,503 circRNAs annotated
- 95.0% annotation recovery
- `host_gene_source` mostly `gtf`

Primary use:

- general circRNA detection and annotation benchmarking
- host-gene landscape
- known versus novel circRNA summary when circAtlas fields are present

## HAP1 scRR-seq

Expected object:

```text
scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu
```

Expected modalities:

- RNA: 63 cells x 63,187 genes
- circRNA: 63 cells x 3,209 circRNAs
- RT: 56 cells x 56,881 genomic bins

Expected overlap:

- RNA/circRNA overlap: 63 cells
- RNA/circRNA/RT overlap: 56 cells

Expected host-gene annotation:

- 3,117 / 3,209 circRNAs annotated
- 97.1% annotation recovery
- `host_gene_source` mostly `gtf`

Known correlations, n = 56 paired cells:

- `frac_rt_pos` vs `circRNA_count`: Pearson r = 0.3778, p = 0.0041; Spearman r = 0.2777, p = 0.0383
- `frac_rt_pos` vs `circRNA_total_support`: Pearson r = 0.3052, p = 0.0222; Spearman r = 0.3081, p = 0.0209
- `frac_rt_neg` vs `circRNA_count`: Pearson r = -0.3887, p = 0.0031; Spearman r = -0.2852, p = 0.0331
- `frac_rt_neg` vs `circRNA_total_support`: Pearson r = -0.3075, p = 0.0211; Spearman r = -0.3105, p = 0.0198
- `detected_genes` vs `circRNA_count`: Pearson r = 0.7357, p = 1.05e-10; Spearman r = 0.6698, p = 1.65e-8

Required next analysis:

```text
circRNA_count ~ detected_genes + frac_rt_pos
```

## IMR90 scRR-seq

Expected object:

```text
scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu
```

Expected modalities:

- RNA: 23 cells x 63,187 genes
- circRNA: 23 cells x 2,443 circRNAs
- CNV: 23 cells x 60,607 genomic bins

Expected overlap:

- RNA/circRNA/CNV overlap: 23 cells

Expected host-gene annotation:

- 2,429 / 2,443 circRNAs annotated
- 99.4% annotation recovery
- `host_gene_source` mostly `gtf`

Expected CNV-high fibroblast/ECM/stress signal genes:

- `COL1A1`
- `FN1`
- `VIM`
- `COL6A3`
- `P4HB`
- `COL1A2`
- `CRIM1`
- `SPARC`
- `LTBP1`
- `CORO1C`
- `ERBIN`
