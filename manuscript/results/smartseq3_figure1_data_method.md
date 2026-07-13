# Smart-seq3 Figure 1 Data Export

Input object: `load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu`

This script exports deterministic, notebook-ready data only. It does not render the final manuscript figure.

RNA preprocessing for `smartseq3_umap_cells.tsv`:

- subset to cells shared by RNA and circRNA modalities;
- `scanpy.pp.normalize_total(target_sum=10000)`;
- `scanpy.pp.log1p`;
- highly variable genes with `flavor="seurat"` and `n_top_genes=2000`;
- PCA with at most 50 components and fixed `random_state=17`;
- neighbors with at most 15 neighbors;
- UMAP with fixed `random_state=17`.

Candidate tables are ranked by prevalence first, then total support. The selected-candidate export reports prevalence, total support, annotation source, sparsity, and per-cell support values so final notebook choices are not forced by total support alone.
