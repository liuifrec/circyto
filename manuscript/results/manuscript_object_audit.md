# Manuscript Object Audit

Read-only audit of the committed manuscript-scale MuData objects. No matrices or metadata were rewritten.

## Summary

| Dataset | Path | Size bytes | SHA-256 | Modalities | Shapes | Shared cells | Host-gene recovery | Local/private findings |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Smart-seq3 | load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu | 30911142 | 0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339 | rna;circ | {"circ":[192,2503],"rna":[192,63187]} | 192 | 2379/2503 (95.0%) | 9 |
| HAP1 | load_work/scrr_hap1/full_length_rna_circ_rt.hostgene_fixed.h5mu | 41952586 | 6d3f460371000dde2f307c67887cfb5939ea591527c4e813d8daf2ef2b1bece5 | rna;circ;rt | {"circ":[63,3209],"rna":[63,63187],"rt":[56,56881]} | 56 | 3117/3209 (97.1%) | 28 |
| IMR90 | load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu | 41036354 | bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8 | rna;circ;cnv | {"circ":[23,2443],"cnv":[23,60607],"rna":[23,63187]} | 23 | 2429/2443 (99.4%) | 30 |

## Smart-seq3

- Accession/source: `E-MTAB-8735`.
- Role: Manuscript-scale RNA+circRNA reproducibility object.
- Public-data assessment: `yes_public_accession_or_canonical_cell_ids`.
- Pairwise shared cells: `{"rna&circ":192}`.
- Median detected circRNAs per cell: `12.0`.
- Median circRNA support per cell: `22.5`.
- Host-gene source counts: `{"gtf":2379,"missing":124}`.
- Missing host-gene provenance fields: `none`.
- Top-level obs columns: `rna:membership, rna:total_rna_counts, rna:detected_genes, rna:mitochondrial_fraction, rna:ribosomal_fraction, rna:circRNA_count, rna:circRNA_total_support, circ:membership, circ:total_rna_counts, circ:detected_genes, circ:mitochondrial_fraction, circ:ribosomal_fraction, circ:circRNA_count, circ:circRNA_total_support, membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support`.
- Top-level obsm keys: `circ, rna`.
- Top-level varm keys: `circ, rna`.
- Top-level uns keys: `circyto`.

### Modalities

- `rna` shape `192 x 63187`; X storage `csr_matrix`, nnz `1695807`, density `0.139781`.
  - obs columns: `membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support`.
  - var columns: `gene_id, gene_name, chrom, start, end, strand, gene_biotype`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `none`.
  - example cells: `DIYHEK_1, DIYHEK_10, DIYHEK_100, DIYHEK_101, DIYHEK_102`.
  - example features: `ENSG00000290825.1, ENSG00000223972.6, ENSG00000227232.6, ENSG00000278267.1, ENSG00000243485.5`.
- `circ` shape `192 x 2503`; X storage `csr_matrix`, nnz `2659`, density `0.00553294`.
  - obs columns: `membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support`.
  - var columns: `circ_id, chrom, start, end, strand, host_gene, host_gene_from_gtf, host_gene_id, host_genes_multi, host_gene_ids_multi, host_gene_n, host_gene_source, host_gene_from_circatlas, host_gene_from_circatlas_id`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `none`.
  - example cells: `DIYHEK_1, DIYHEK_10, DIYHEK_100, DIYHEK_101, DIYHEK_102`.
  - example features: `chr19:46910586|46911183, chr9:33953285|33996333, chr9:33960826|33989126, chrX:101049107|101050336, chr15:65011401|65012591`.

### Compression And Stored Matrices

- HDF5 datasets: `101`; compressed datasets: `0`; algorithms: `{}`.
- Matrix datasets: `[{"chunks":"2659","compression":"","dataset":"mod/circ/X/data","dtype":"int32","shape":"2659"},{"chunks":"6625","compression":"","dataset":"mod/rna/X/data","dtype":"int32","shape":"1695807"}]`.
- Manuscript necessity assessment: RNA `X`, circRNA `X`, and RT/CNV `X` are directly needed for manuscript validation. Any CNV `mappabilitynorm` layer is useful provenance for processed public CNV summaries. No object rewrite was performed.

### Local/Private Metadata Scan

| Location | Category | Value |
| --- | --- | --- |
| uns/circyto/rna_circ_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_circ_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_import_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_import_summary_json | known_local_username | [local username redacted] |
| uns/circyto/source_workdir | server_path | [local path redacted] |
| uns/circyto/source_workdir | known_local_username | [local username redacted] |
| uns/circyto/source_workdir | known_local_username | [local username redacted] |
| uns/circyto/workflow_summary_json | known_local_username | [local username redacted] |
| uns/circyto/workflow_summary_json | known_local_username | [local username redacted] |

## HAP1

- Accession/source: `GSE278952`.
- Role: Manuscript-scale RNA+circRNA+RT reproducibility object.
- Public-data assessment: `yes_public_accession_or_canonical_cell_ids`.
- Pairwise shared cells: `{"circ&rt":56,"rna&circ":63,"rna&rt":56}`.
- Median detected circRNAs per cell: `50.0`.
- Median circRNA support per cell: `371.0`.
- Host-gene source counts: `{"gtf":3117,"missing":92}`.
- Missing host-gene provenance fields: `none`.
- Top-level obs columns: `rna:membership, rna:total_rna_counts, rna:detected_genes, rna:mitochondrial_fraction, rna:ribosomal_fraction, rna:circRNA_count, rna:circRNA_total_support, rna:original_obs_id, rna:gsm_id, rna:canonical_cell_id, rna:rna_cell_id, rna:dna_cell_id, rna:sample_title, rna:molecule, rna:original_rna_obs_id, circ:membership, circ:total_rna_counts, circ:detected_genes, circ:mitochondrial_fraction, circ:ribosomal_fraction, circ:circRNA_count, circ:circRNA_total_support, circ:original_obs_id, circ:gsm_id, circ:canonical_cell_id, circ:rna_cell_id, circ:dna_cell_id, circ:sample_title, circ:molecule, circ:original_circ_obs_id, rt:cell_id, rt:canonical_cell_id, rt:dna_cell_id, rt:rna_cell_id, rt:source_cell_id, rt:source_molecule, rt:pairing_strategy`.
- Top-level obsm keys: `circ, rna, rt`.
- Top-level varm keys: `circ, rna, rt`.
- Top-level uns keys: `circyto, circyto_scrr_remap, circyto_scrr_rt_trimodal`.

### Modalities

- `rna` shape `63 x 63187`; X storage `csr_matrix`, nnz `661042`, density `0.166058`.
  - obs columns: `membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support, original_obs_id, gsm_id, canonical_cell_id, rna_cell_id, dna_cell_id, sample_title, molecule, original_rna_obs_id`.
  - var columns: `gene_id, gene_name, chrom, start, end, strand, gene_biotype`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `none`.
  - example cells: `HAP1_scRR1_MidS_104, HAP1_scRR1_MidS_103, HAP1_scRR1_MidS_102, HAP1_scRR1_MidS_101, HAP1_scRR1_MidS_096`.
  - example features: `ENSG00000290825.1, ENSG00000223972.6, ENSG00000227232.6, ENSG00000278267.1, ENSG00000243485.5`.
- `circ` shape `63 x 3209`; X storage `csr_matrix`, nnz `3674`, density `0.0181731`.
  - obs columns: `membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support, original_obs_id, gsm_id, canonical_cell_id, rna_cell_id, dna_cell_id, sample_title, molecule, original_circ_obs_id`.
  - var columns: `circ_id, chrom, start, end, strand, host_gene, host_gene_from_gtf, host_gene_id, host_genes_multi, host_gene_ids_multi, host_gene_n, host_gene_source, host_gene_from_circatlas, host_gene_from_circatlas_id`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `none`.
  - example cells: `HAP1_scRR1_MidS_104, HAP1_scRR1_MidS_103, HAP1_scRR1_MidS_102, HAP1_scRR1_MidS_101, HAP1_scRR1_MidS_096`.
  - example features: `chr9:27482842|27482998, chr9:131430090|131432794, chr7:5596478|5596886, chr7:72879795|72891083, chr7:100358690|100373720`.
- `rt` shape `56 x 56881`; X storage `ndarray`, nnz `3165828`, density `0.993876`.
  - obs columns: `cell_id, canonical_cell_id, dna_cell_id, rna_cell_id, source_cell_id, source_molecule, pairing_strategy`.
  - var columns: `chr, start, end, seqname, bin_size, feature_id, feature_type`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `circyto`.
  - example cells: `HAP1_scRR1_MidS_031, HAP1_scRR1_MidS_051, HAP1_scRR1_MidS_024, HAP1_scRR1_MidS_023, HAP1_scRR1_MidS_053`.
  - example features: `chr2:0-40000, chr2:40000-80000, chr2:80000-120000, chr2:120000-160000, chr2:160000-200000`.

### Compression And Stored Matrices

- HDF5 datasets: `171`; compressed datasets: `0`; algorithms: `{}`.
- Matrix datasets: `[{"chunks":"1837","compression":"","dataset":"mod/circ/X/data","dtype":"int32","shape":"3674"},{"chunks":"5165","compression":"","dataset":"mod/rna/X/data","dtype":"int32","shape":"661042"},{"chunks":"","compression":"","dataset":"mod/rt/X","dtype":"int16","shape":"56x56881"}]`.
- Manuscript necessity assessment: RNA `X`, circRNA `X`, and RT/CNV `X` are directly needed for manuscript validation. Any CNV `mappabilitynorm` layer is useful provenance for processed public CNV summaries. No object rewrite was performed.

### Local/Private Metadata Scan

| Location | Category | Value |
| --- | --- | --- |
| mod/rt/uns/circyto/avg_rt_bedgraph/source_avg_rt_bedgraph | server_path | [local path redacted] |
| mod/rt/uns/circyto/avg_rt_bedgraph/source_avg_rt_bedgraph | known_local_username | [local username redacted] |
| mod/rt/uns/circyto/avg_rt_bedgraph/source_avg_rt_bedgraph | known_local_username | [local username redacted] |
| mod/rt/uns/circyto/source_avg_rt_bedgraph | server_path | [local path redacted] |
| mod/rt/uns/circyto/source_avg_rt_bedgraph | known_local_username | [local username redacted] |
| mod/rt/uns/circyto/source_avg_rt_bedgraph | known_local_username | [local username redacted] |
| mod/rt/uns/circyto/source_rt_table | server_path | [local path redacted] |
| mod/rt/uns/circyto/source_rt_table | known_local_username | [local username redacted] |
| mod/rt/uns/circyto/source_rt_table | known_local_username | [local username redacted] |
| uns/circyto/rna_import_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_import_summary_json | known_local_username | [local username redacted] |
| uns/circyto/source_workdir | server_path | [local path redacted] |
| uns/circyto/source_workdir | known_local_username | [local username redacted] |
| uns/circyto/source_workdir | known_local_username | [local username redacted] |
| uns/circyto/workflow_summary_json | known_local_username | [local username redacted] |
| uns/circyto/workflow_summary_json | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/cell_map | server_path | [local path redacted] |
| uns/circyto_scrr_remap/cell_map | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/cell_map | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/source_h5mu | server_path | [local path redacted] |
| uns/circyto_scrr_remap/source_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/source_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_rt_trimodal/source_rna_circ_h5mu | server_path | [local path redacted] |
| uns/circyto_scrr_rt_trimodal/source_rna_circ_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_rt_trimodal/source_rna_circ_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_rt_trimodal/source_rt_h5ad | server_path | [local path redacted] |
| uns/circyto_scrr_rt_trimodal/source_rt_h5ad | known_local_username | [local username redacted] |
| uns/circyto_scrr_rt_trimodal/source_rt_h5ad | known_local_username | [local username redacted] |

## IMR90

- Accession/source: `GSE278958`.
- Role: Manuscript-scale RNA+circRNA+CNV reproducibility object.
- Public-data assessment: `yes_public_accession_or_canonical_cell_ids`.
- Pairwise shared cells: `{"circ&cnv":23,"rna&circ":23,"rna&cnv":23}`.
- Median detected circRNAs per cell: `119.0`.
- Median circRNA support per cell: `317.0`.
- Host-gene source counts: `{"gtf":2429,"missing":14}`.
- Missing host-gene provenance fields: `none`.
- Top-level obs columns: `rna:membership, rna:total_rna_counts, rna:detected_genes, rna:mitochondrial_fraction, rna:ribosomal_fraction, rna:circRNA_count, rna:circRNA_total_support, rna:original_obs_id, rna:gsm_id, rna:canonical_cell_id, rna:rna_cell_id, rna:dna_cell_id, rna:sample_title, rna:molecule, rna:original_rna_obs_id, circ:membership, circ:total_rna_counts, circ:detected_genes, circ:mitochondrial_fraction, circ:ribosomal_fraction, circ:circRNA_count, circ:circRNA_total_support, circ:original_obs_id, circ:gsm_id, circ:canonical_cell_id, circ:rna_cell_id, circ:dna_cell_id, circ:sample_title, circ:molecule, circ:original_circ_obs_id, cnv:cell_id, cnv:canonical_cell_id, cnv:dna_cell_id, cnv:rna_cell_id, cnv:pairing_strategy, membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support, original_obs_id, gsm_id, canonical_cell_id, rna_cell_id, dna_cell_id, sample_title, molecule, cell_id, pairing_strategy`.
- Top-level obsm keys: `circ, cnv, rna`.
- Top-level varm keys: `circ, cnv, rna`.
- Top-level uns keys: `circyto, circyto_scrr_remap, circyto_scrr_trimodal`.

### Modalities

- `rna` shape `23 x 63187`; X storage `csr_matrix`, nnz `376158`, density `0.25883`.
  - obs columns: `membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support, original_obs_id, gsm_id, canonical_cell_id, rna_cell_id, dna_cell_id, sample_title, molecule, original_rna_obs_id`.
  - var columns: `gene_id, gene_name, chrom, start, end, strand, gene_biotype`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `none`.
  - example cells: `IMR90_A_126, IMR90_A_125, IMR90_A_122, IMR90_A_119, IMR90_A_117`.
  - example features: `ENSG00000290825.1, ENSG00000223972.6, ENSG00000227232.6, ENSG00000278267.1, ENSG00000243485.5`.
- `circ` shape `23 x 2443`; X storage `csr_matrix`, nnz `2914`, density `0.0518607`.
  - obs columns: `membership, total_rna_counts, detected_genes, mitochondrial_fraction, ribosomal_fraction, circRNA_count, circRNA_total_support, original_obs_id, gsm_id, canonical_cell_id, rna_cell_id, dna_cell_id, sample_title, molecule, original_circ_obs_id`.
  - var columns: `circ_id, chrom, start, end, strand, host_gene, host_gene_from_gtf, host_gene_id, host_genes_multi, host_gene_ids_multi, host_gene_n, host_gene_source, host_gene_from_circatlas, host_gene_from_circatlas_id`.
  - obsm keys: `none`; varm keys: `none`; layers: `None`; uns keys: `none`.
  - example cells: `IMR90_A_126, IMR90_A_125, IMR90_A_122, IMR90_A_119, IMR90_A_117`.
  - example features: `chr9:110941421|111011690, chr9:110972073|111011690, chr9:110972073|110973558, chr9:131643239|131650949, chr7:849373|849578`.
- `cnv` shape `23 x 60607`; X storage `ndarray`, nnz `1278868`, density `0.917435`.
  - obs columns: `cell_id, canonical_cell_id, dna_cell_id, rna_cell_id, pairing_strategy`.
  - var columns: `bin_id, seqname, start, end, bin_size`.
  - obsm keys: `none`; varm keys: `none`; layers: `mappabilitynorm, None`; uns keys: `circyto`.
  - example cells: `IMR90_A_126, IMR90_A_125, IMR90_A_122, IMR90_A_119, IMR90_A_117`.
  - example features: `chr1:0-50000, chr1:50000-100000, chr1:100000-150000, chr1:150000-200000, chr1:200000-250000`.

### Compression And Stored Matrices

- HDF5 datasets: `177`; compressed datasets: `0`; algorithms: `{}`.
- Matrix datasets: `[{"chunks":"2914","compression":"","dataset":"mod/circ/X/data","dtype":"int32","shape":"2914"},{"chunks":"","compression":"","dataset":"mod/cnv/X","dtype":"int16","shape":"23x60607"},{"chunks":"","compression":"","dataset":"mod/cnv/layers/mappabilitynorm","dtype":"float32","shape":"23x60607"},{"chunks":"","compression":"","dataset":"mod/cnv/uns/circyto/layers/mappabilitynorm","dtype":"object","shape":""},{"chunks":"5878","compression":"","dataset":"mod/rna/X/data","dtype":"int32","shape":"376158"}]`.
- Manuscript necessity assessment: RNA `X`, circRNA `X`, and RT/CNV `X` are directly needed for manuscript validation. Any CNV `mappabilitynorm` layer is useful provenance for processed public CNV summaries. No object rewrite was performed.

### Local/Private Metadata Scan

| Location | Category | Value |
| --- | --- | --- |
| mod/cnv/uns/circyto/layers/mappabilitynorm | server_path | [local path redacted] |
| mod/cnv/uns/circyto/layers/mappabilitynorm | known_local_username | [local username redacted] |
| mod/cnv/uns/circyto/layers/mappabilitynorm | known_local_username | [local username redacted] |
| mod/cnv/uns/circyto/source_cnv_mappabilitynorm | server_path | [local path redacted] |
| mod/cnv/uns/circyto/source_cnv_mappabilitynorm | known_local_username | [local username redacted] |
| mod/cnv/uns/circyto/source_cnv_mappabilitynorm | known_local_username | [local username redacted] |
| mod/cnv/uns/circyto/source_cnv_states | server_path | [local path redacted] |
| mod/cnv/uns/circyto/source_cnv_states | known_local_username | [local username redacted] |
| mod/cnv/uns/circyto/source_cnv_states | known_local_username | [local username redacted] |
| uns/circyto/rna_circ_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_circ_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_import_summary_json | known_local_username | [local username redacted] |
| uns/circyto/rna_import_summary_json | known_local_username | [local username redacted] |
| uns/circyto/source_workdir | server_path | [local path redacted] |
| uns/circyto/source_workdir | known_local_username | [local username redacted] |
| uns/circyto/source_workdir | known_local_username | [local username redacted] |
| uns/circyto/workflow_summary_json | known_local_username | [local username redacted] |
| uns/circyto/workflow_summary_json | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/cell_map | server_path | [local path redacted] |
| uns/circyto_scrr_remap/cell_map | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/cell_map | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/source_h5mu | server_path | [local path redacted] |
| uns/circyto_scrr_remap/source_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_remap/source_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_trimodal/source_cnv_h5ad | server_path | [local path redacted] |
| uns/circyto_scrr_trimodal/source_cnv_h5ad | known_local_username | [local username redacted] |
| uns/circyto_scrr_trimodal/source_cnv_h5ad | known_local_username | [local username redacted] |
| uns/circyto_scrr_trimodal/source_rna_circ_h5mu | server_path | [local path redacted] |
| uns/circyto_scrr_trimodal/source_rna_circ_h5mu | known_local_username | [local username redacted] |
| uns/circyto_scrr_trimodal/source_rna_circ_h5mu | known_local_username | [local username redacted] |
