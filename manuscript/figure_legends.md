# Application Note figure legends

**Figure 1. circyto workflow and annotated single-cell circRNA outputs.**
(A) circyto orchestrates protocol-aware alignment and established circRNA
detectors, constructs cell-by-circRNA candidate matrices with QC and annotation,
and exports AnnData/MuData for scverse analysis. circyto is not a new back-splice
detection algorithm. Optional processed CNV, RT, and exploratory candidate-variant
integrations have the evidence boundaries shown. (B) RNA-derived UMAP of 192
Smart-seq3 cells from E-MTAB-8735, colored by the number of circRNA candidates
with nonzero support. RNA preprocessing used normalization to 10,000 counts,
log1p, 2,000 variable genes, 50 PCs, 15 neighbors, and seed 17. CircRNA features
did not contribute to the embedding. The object contains 63,187 RNA features,
2,503 circRNA candidates, and 2,659 nonzero circRNA matrix entries. Host genes
are annotated for 2,379 candidates (95.0%). Median detected circRNA candidates
and total support per cell are 12 and 22.5, respectively. (C) The same RNA UMAP
coordinates overlaid with detection of `chr1:117402186|117420649`, annotated to
MAN1A2: 6/192 cells, total support 15. Gray denotes no observed support under
this assay; it does not establish biological absence. This candidate illustrates
the output interface and does not establish MAN1A2 regulation, function, or
independent circRNA validation. Full methods and input/output checksums are in
`results_v2/provenance.json`.

**Supplementary Figure S1. Protocol/workflow validation.**
The Smart-seq3 paired-end/indexed route, IMR90 single-end scRR route, and HAP1
paired-end scRR route use the indicated aligners and CIRI3 modes. Smart-seq3
and IMR90 object quantities are independently regenerated from checksum-matched
processed objects. HAP1 refers to the committed real-data 10-cell RNA/circ route;
its unresolved historical full-object/RT results are excluded. Execution and
interoperability validation do not establish detector accuracy.

**Supplementary Figure S2. Multimodal object organization.**
Each AnnData modality retains its own cell and feature axes. The Smart-seq3
RNA+circ object shares 192 cells. The IMR90 object contains RNA (23 × 63,187),
circRNA (23 × 2,443), and processed CNV (23 × 60,607) with 23 shared cells;
2,429/2,443 circRNA candidates have host-gene annotations (99.4%). CNV is
imported from processed GEO summaries and is used only to demonstrate
interoperability. The dashed RT schema is an implemented contract validated
with synthetic fixtures; reconciled real processed-file validation remains
pending. No RT or CNV biological relationship is inferred.
