# Manuscript Consistency Report

## Objects inspected

- Smart-seq3: `/mnt/d/circyto/load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu`
- HAP1: `/mnt/d/circyto/load_work/scrr_hap1/full_length_rna_circ_rt.hostgene_fixed.h5mu`
- IMR90: `/mnt/d/circyto/load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu`

## Object inventory

| dataset | modalities | modality_n_obs | modality_n_vars | all_modality_cell_overlap | pairwise_cell_overlap |
| --- | --- | --- | --- | --- | --- |
| Smart-seq3 | circ;rna | {"circ":192,"rna":192} | {"circ":2503,"rna":63187} | 192 | {"circ&rna":192} |
| HAP1 | circ;rna;rt | {"circ":63,"rna":63,"rt":56} | {"circ":3209,"rna":63187,"rt":56881} | 56 | {"circ&rna":63,"circ&rt":56,"rna&rt":56} |
| IMR90 | circ;cnv;rna | {"circ":23,"cnv":23,"rna":23} | {"circ":2443,"cnv":60607,"rna":63187} | 23 | {"circ&cnv":23,"circ&rna":23,"cnv&rna":23} |

## Host-gene annotation

| dataset | n_circRNAs | host_gene_annotated | host_gene_unannotated | host_gene_annotation_rate | n_unique_host_genes | host_gene_source_counts |
| --- | --- | --- | --- | --- | --- | --- |
| Smart-seq3 | 2503 | 2379 | 124 | 0.950459448661606 | 2199 | {"gtf":2379,"missing":124} |
| HAP1 | 3209 | 3117 | 92 | 0.9713306325958242 | 2746 | {"gtf":3117,"missing":92} |
| IMR90 | 2443 | 2429 | 14 | 0.994269340974212 | 2005 | {"gtf":2429,"missing":14} |

## Manuscript number checks

| dataset | metric | manuscript_value | object_value | status |
| --- | --- | --- | --- | --- |
| Smart-seq3 | cells | 192 | 192 | match |
| Smart-seq3 | RNA features | 63187 | 63187 | match |
| Smart-seq3 | circRNAs | 2503 | 2503 | match |
| Smart-seq3 | host annotation rate | 95.0% | 95.0% (2379/2503) | match |
| Smart-seq3 | median circRNAs/cell | not documented | 12.0 | needs_manuscript_value |
| Smart-seq3 | median support/cell | not documented | 22.5 | needs_manuscript_value |
| HAP1 | cells | 63 | 63 | match |
| HAP1 | circRNAs | 3209 | 3209 | match |
| HAP1 | RT overlap | 56 | 56 | match |
| HAP1 | host annotation rate | 97.1% | 97.1% (3117/3209) | match |
| IMR90 | cells | 23 | 23 | match |
| IMR90 | circRNAs | 2443 | 2443 | match |
| IMR90 | CNV overlap | 23 | 23 | match |
| IMR90 | host annotation rate | 99.4% | 99.4% (2429/2443) | match |

## Smart-seq3 Figure 1 data

- RNA UMAP in object: `no`
- UMAP coordinate source: `computed_internal_exact_knn_umap`
- Final figure design and plotting were not performed.

## Candidate Figure 1 panel C features

| role | feature_type | feature_id | host_gene | n_cells_detected | total_support | combined_score | umap_detected_spread_ratio | umap_detected_area_fraction | selection_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| individual_circRNA | circRNA | chr7:158869855\|158876691 | DYNC2I1 | 10 | 28.0 | 0.9992209348781462 | 0.9383623694923711 | 0.5255601849366683 | Top individual circRNA by combined recurrence/support score with at least 6 detected cells. |
| host_gene_signal | host_gene | MAN1A2 | MAN1A2 | 12 | 31.0 | 0.9921100500227377 | 1.103220365227557 | 0.8477904816167132 | Top host-gene-level circRNA signal by recurrence/support/feature-count score with at least 6 detected cells. |
| backup_candidate | circRNA | chr1:31915895\|31919658 | PTP4A2;ENSG00000269967;ENSG00000288678 | 14 | 26.0 | 0.9990011985617259 | 0.7572557980000445 | 0.5898777277746798 | Second-ranked individual circRNA candidate by the same data-driven score. |

## Metadata columns inspected

### Smart-seq3

- obs columns by modality: `{"circ":["membership","total_rna_counts","detected_genes","mitochondrial_fraction","ribosomal_fraction","circRNA_count","circRNA_total_support"],"rna":["membership","total_rna_counts","detected_genes","mitochondrial_fraction","ribosomal_fraction","circRNA_count","circRNA_total_support"]}`
- var columns by modality: `{"circ":["circ_id","chrom","start","end","strand","host_gene","host_gene_from_gtf","host_gene_id","host_genes_multi","host_gene_ids_multi","host_gene_n","host_gene_source","host_gene_from_circatlas","host_gene_from_circatlas_id"],"rna":["gene_id","gene_name","chrom","start","end","strand","gene_biotype"]}`

### HAP1

- obs columns by modality: `{"circ":["membership","total_rna_counts","detected_genes","mitochondrial_fraction","ribosomal_fraction","circRNA_count","circRNA_total_support","original_obs_id","gsm_id","canonical_cell_id","rna_cell_id","dna_cell_id","sample_title","molecule","original_circ_obs_id"],"rna":["membership","total_rna_counts","detected_genes","mitochondrial_fraction","ribosomal_fraction","circRNA_count","circRNA_total_support","original_obs_id","gsm_id","canonical_cell_id","rna_cell_id","dna_cell_id","sample_title","molecule","original_rna_obs_id"],"rt":["cell_id","canonical_cell_id","dna_cell_id","rna_cell_id","source_cell_id","source_molecule","pairing_strategy"]}`
- var columns by modality: `{"circ":["circ_id","chrom","start","end","strand","host_gene","host_gene_from_gtf","host_gene_id","host_genes_multi","host_gene_ids_multi","host_gene_n","host_gene_source","host_gene_from_circatlas","host_gene_from_circatlas_id"],"rna":["gene_id","gene_name","chrom","start","end","strand","gene_biotype"],"rt":["chr","start","end","seqname","bin_size","feature_id","feature_type"]}`

### IMR90

- obs columns by modality: `{"circ":["membership","total_rna_counts","detected_genes","mitochondrial_fraction","ribosomal_fraction","circRNA_count","circRNA_total_support","original_obs_id","gsm_id","canonical_cell_id","rna_cell_id","dna_cell_id","sample_title","molecule","original_circ_obs_id"],"cnv":["cell_id","canonical_cell_id","dna_cell_id","rna_cell_id","pairing_strategy"],"rna":["membership","total_rna_counts","detected_genes","mitochondrial_fraction","ribosomal_fraction","circRNA_count","circRNA_total_support","original_obs_id","gsm_id","canonical_cell_id","rna_cell_id","dna_cell_id","sample_title","molecule","original_rna_obs_id"]}`
- var columns by modality: `{"circ":["circ_id","chrom","start","end","strand","host_gene","host_gene_from_gtf","host_gene_id","host_genes_multi","host_gene_ids_multi","host_gene_n","host_gene_source","host_gene_from_circatlas","host_gene_from_circatlas_id"],"cnv":["bin_id","seqname","start","end","bin_size"],"rna":["gene_id","gene_name","chrom","start","end","strand","gene_biotype"]}`

## Summary

No documented count mismatches were found against the inspected objects.
2 Smart-seq3 median value(s) were not documented in the current manuscript markdown and should be inserted from this report/table if needed.
