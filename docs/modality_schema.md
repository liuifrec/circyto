# Modality Schema

This page describes the expected `circyto` modality schema for manuscript-scale AnnData and MuData objects.

## RNA Modality

Expected modality name: `rna`.

Shape:

- observations: cells
- variables: genes

Expected matrix:

- `X`: gene expression counts or imported gene-count values, cells x genes.

Common `obs` fields:

- `cell_id` or the object index as cell ID
- `total_counts`
- `detected_genes`
- `mitochondrial_fraction`
- `ribosomal_fraction`
- optional workflow membership fields when RNA and circRNA axes differ

Common `var` fields:

- `gene_id`
- `gene_name`
- `gene_biotype`

## circRNA Modality

Expected modality name: `circ`.

Shape:

- observations: cells
- variables: circRNAs

Expected matrix:

- `X`: back-splice junction support/count matrix, cells x circRNAs.

Expected circRNA identifiers:

- `var_names` are circRNA IDs.
- IDs may be detector-native or normalized coordinate IDs.

Expected coordinate fields when available:

- `chrom`, `start`, `end`, `strand`
- accepted synonyms in manuscript scripts include `seqname`, `seqnames`, `chr`, `chromosome`, `donor`, and `acceptor`

Expected host-gene fields:

- `host_gene`
- `host_gene_source`
- `host_gene_from_gtf`
- `host_gene_from_circatlas`
- `host_gene_from_circatlas_id`
- `host_gene_id`
- `host_genes_multi`
- `host_gene_ids_multi`
- `host_gene_n`

Known circRNA annotation fields are database-dependent, but may include:

- `annotation_status`
- `known_status`
- `circatlas_id`
- database-specific matched IDs or labels

Expected provenance:

- `uns["circyto"]` records workflow, command, version, source paths, and related summaries where available.

## RT Modality

Expected modality name: `rt`.

Shape:

- observations: cells
- variables: genomic bins or gene-intersect genomic features

Expected matrix:

- `X`: processed replication timing/state values.

Common `obs` fields:

- `cell_id`
- `canonical_cell_id`
- `dna_cell_id`
- `rna_cell_id`
- `frac_rt_pos`
- `frac_rt_neg`

Common `var` fields:

- `feature_id`
- `feature_type`
- `seqname`
- `start`
- `end`
- optional `avg_rt`

## CNV Modality

Expected modality name: `cnv`.

Shape:

- observations: cells
- variables: genomic bins

Expected matrix:

- `X`: integer copy-number state or processed copy-number state value.
- optional `layers["mappabilitynorm"]`: mappability-normalized values from processed scRR summaries.

Common `obs` fields:

- `cell_id`
- `canonical_cell_id`
- `dna_cell_id`
- `rna_cell_id`
- `frac_non_diploid`
- `frac_loss`
- `frac_gain`

Common `var` fields:

- `feature_id`
- `seqname`
- `start`
- `end`
- `bin_size`

## candidate_snv Modality

Expected modality name: `candidate_snv`.

This layer is exploratory and represents RNA-derived candidate variant signals imported from SComatic-style outputs. It should not be interpreted as orthogonally confirmed DNA variation.

Expected matrix:

- `X`: candidate signal counts or normalized candidate counts, depending on import mode.

Common metadata:

- candidate locus fields
- candidate class fields
- import provenance
- filtering thresholds used by the external tool

Use conservative language when reporting this layer. The manuscript should present it as interoperability and hypothesis generation unless independent validation is added.
