# Manuscript Result Tables

`load_work/` contains local/private analysis inputs and is excluded from git. The committed files in this directory are lightweight summaries derived from those inputs; the source `.h5mu` objects remain private.

Public Smart-seq3, HAP1 scRR-seq, and IMR90 scRR-seq data are proof-of-concept multimodal analysis inputs. They are not radiation-exposure datasets.

## Generated Tables

| Table | Script | Rows x columns | Notes |
| --- | --- | ---: | --- |
| `mudata_inventory.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 3 x 26 | Expanded modality inventory, cell overlap, obs/var columns, circRNA feature counts, and host-gene recovery. |
| `host_gene_annotation_summary.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 3 x 13 | Host-gene annotation recovery and provenance summary. |
| `top_host_genes_by_dataset.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 300 x 19 | Top 100 host-gene-level circRNA features per dataset. |
| `smartseq3_umap_cells.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 192 x 12 | Figure-ready Smart-seq3 cell coordinates and per-cell circRNA burden metrics. |
| `smartseq3_top_circRNA_candidates.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 100 x 20 | Ranked Smart-seq3 individual circRNA candidates. |
| `smartseq3_top_hostgene_features.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 100 x 18 | Ranked Smart-seq3 host-gene-level circRNA features with UMAP spread metrics. |
| `smartseq3_selected_feature_candidates.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 3 x 10 | Candidate Figure 1 panel C features selected from recurrence, support, host annotation, and UMAP spread. |
| `supplement_table_1_dataset_inventory.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 3 x 26 | Supplement table S1 copy of the expanded dataset inventory. |
| `supplement_table_2_detector_status.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 5 x 9 | Detector runtime readiness in the current environment. |
| `supplement_table_3_top_host_genes.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 300 x 19 | Supplement table S3 copy of top host genes by dataset. |
| `supplement_table_4_selected_smartseq3_features.tsv` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | 3 x 10 | Supplement table S4 copy of selected Smart-seq3 feature candidates. |
| `manuscript_consistency_report.md` | `scripts/manuscript/prepare_manuscript_consistency_outputs.py` | markdown | Consistency report comparing current manuscript numbers with inspected objects. |
| `hap1_rt_circ_correlations.tsv` | `scripts/manuscript/hap1_rt_circ_analysis.py` | 5 x 8 | HAP1 RT fractions and detected genes versus circRNA burden. |
| `hap1_rt_circ_ols.tsv` | `scripts/manuscript/hap1_rt_circ_analysis.py` | 3 x 8 | OLS for `circRNA_count ~ detected_genes + frac_rt_pos`. |
| `hap1_rt_circ_cell_metrics.tsv` | `scripts/manuscript/hap1_rt_circ_analysis.py` | 56 x 6 | Per-cell HAP1 metrics for plotting. |
| `imr90_cnv_circ_correlations.tsv` | `scripts/manuscript/imr90_cnv_circ_analysis.py` | 4 x 8 | Global CNV burden versus circRNA burden. Exploratory. |
| `imr90_cnv_high_host_enrichment.tsv` | `scripts/manuscript/imr90_cnv_circ_analysis.py` | 2005 x 13 | CNV-high versus CNV-low host-gene support using `frac_non_diploid`. Exploratory. |
| `imr90_cnv_circ_cell_metrics.tsv` | `scripts/manuscript/imr90_cnv_circ_analysis.py` | 23 x 6 | Per-cell IMR90 metrics for plotting. Exploratory. |
| `imr90_local_cnv_circ.tsv` | `scripts/manuscript/imr90_cnv_circ_analysis.py` | 2442 x 13 | Local CNV state at circRNA loci. Exploratory. |
| `cross_dataset_host_overlap.tsv` | `scripts/manuscript/cross_dataset_host_overlap.py` | 3 x 7 | Pairwise host-gene overlap across Smart-seq3, HAP1, and IMR90. |
| `cross_dataset_threeway_hosts.tsv` | `scripts/manuscript/cross_dataset_host_overlap.py` | 264 x 2 | Host genes observed in all three datasets. |

`cross_dataset_shared_positive_hosts.tsv` was not generated because matched positive host-gene enrichment inputs for both HAP1 and IMR90 were not available. `known_novel_circ_summary.tsv` was skipped because the circRNA feature annotations did not include known/novel fields such as `is_known`, `known_status`, `known_novel`, `annotation_id`, or matched database IDs.

## Key Results

The MuData inventory found complete RNA/circRNA cell overlap for Smart-seq3 (192 cells), HAP1 RNA/circRNA overlap of 63 cells with 56 cells shared with RT, and complete IMR90 RNA/circRNA/CNV overlap across 23 cells. CircRNA feature counts were 2503 for Smart-seq3, 3209 for HAP1, and 2443 for IMR90. Host-gene recovery was 95.0%, 97.1%, and 99.4%, respectively. Smart-seq3 median per-cell circRNA count was 12.0 and median per-cell circRNA support was 22.5 in the inspected object.

The Smart-seq3 object did not contain an RNA UMAP in `rna.obsm`; `smartseq3_umap_cells.tsv` therefore contains deterministic RNA-derived coordinates computed for notebook plotting. These are figure-ready data only, not a final styled figure.

HAP1 showed positive marginal correlations between `frac_rt_pos` and circRNA burden (`circRNA_count` Pearson r = 0.378, p = 0.00410; `circRNA_total_support` Pearson r = 0.305, p = 0.0222). Detected genes were more strongly associated with `circRNA_count` (Pearson r = 0.736, p = 1.05e-10). In OLS, `frac_rt_pos` did not remain significantly associated after controlling for detected genes: beta = 42.99, p = 0.244, n = 56. This does not support an independent RT association in this model.

IMR90 global CNV burden correlations were weak and not significant in 23 cells: `frac_non_diploid` versus `circRNA_count` Pearson r = -0.151, p = 0.490; `frac_non_diploid` versus `circRNA_total_support` Pearson r = 0.0656, p = 0.766; `frac_loss` versus `circRNA_count` Pearson r = -0.146, p = 0.508; `frac_gain` versus `circRNA_count` Pearson r = -0.0870, p = 0.693.

The IMR90 CNV-high versus CNV-low host-gene analysis is exploratory. With 6 CNV-high and 6 CNV-low cells, the largest positive differences included `Y_RNA`, `TUBA1A`, `TUBA1B`, `LDHA`, `LDHC`, and `COL1A1`. The local CNV-at-circRNA-locus analysis was possible because circRNA features included `chrom`, `start`, and `end`, and CNV bins included `seqname`, `start`, and `end`. It produced 2442 locus rows, 1018 rows with non-missing Pearson correlations, and 73 nominal p < 0.05 rows; these local results should be treated as exploratory because of the small cell number and sparse circRNA support.

Pairwise host-gene overlaps were 645 for Smart-seq3/HAP1, 507 for Smart-seq3/IMR90, and 720 for HAP1/IMR90. The three-way host-gene overlap contained 264 genes.
