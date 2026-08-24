# circyto-biogenesis schema foundation

`circyto-biogenesis` schema v1.0 defines dependency-light, model-ready tables
without adding a modelling framework to the main `circyto` package. The
foundation consists of three related UTF-8 TSV tables and one JSON provenance
sidecar.

This contract does not alter existing AnnData or MuData outputs. It is an
explicit export layer over data derived from those objects and their workflow
artifacts.

## Versioning

Every record has `schema_version=1.0`. A breaking field or semantic change
requires a new schema version. Additional modelling columns may be added
without breaking v1 validation when the required v1 fields retain their
meaning.

Canonical v1 filenames are:

- `circRNA_candidates.v1.tsv`
- `cell_contexts.v1.tsv`
- `cell_circ_observations.v1.tsv`
- `biogenesis_provenance.v1.json`

## Candidate records

One row represents one candidate circRNA.

Required fields:

| Field | Meaning |
| --- | --- |
| `schema_version` | Schema version; exactly `1.0`. |
| `circ_id` | Unique candidate identifier within the bundle. |
| `dataset_id` | Dataset or study identifier. |
| `chrom`, `start`, `end`, `strand` | Genomic back-splice interval and strand. |
| `coordinate_system` | `0-based-half-open` or `1-based-closed`. |
| `detector` | Detector or candidate generator. |
| `reference_genome` | Genome assembly/build identifier. |
| `reference_annotation` | Gene/transcript annotation release. |
| `workflow_uuid` | Stable originating workflow identifier. |
| `label_status` | `positive` or `unlabelled`; never an inferred negative. |

Optional fields:

- host annotation: `host_gene`, `host_gene_id`
- sequence windows: `donor_window`, `acceptor_window`
- geometry: `exon_count`, `circ_length`, `upstream_exon_length`,
  `downstream_exon_length`, `upstream_intron_length`,
  `downstream_intron_length`
- repeats: `upstream_repeat_family`, `downstream_repeat_family`,
  `repeat_pairing_score`
- splice sites: `donor_splice_site`, `acceptor_splice_site`,
  `donor_splice_score`, `acceptor_splice_score`
- database annotation: `known_status`, `known_database`,
  `known_database_id`

`known_status` may be `known`, `novel`, or `unknown`. Database novelty is an
annotation and is independent of the modelling label.

## Cell-context records

One row represents the technical and biological context of one cell.

Required fields:

| Field | Meaning |
| --- | --- |
| `schema_version`, `cell_id` | Schema version and unique cell identifier. |
| `dataset_id`, `donor_id` | Dataset and biological donor provenance. |
| `protocol` | Library protocol, such as `ramda` or `smartseq3`. |
| `detector` | Detector associated with this exported bundle. |
| `reference_genome`, `reference_annotation` | Reference provenance. |
| `workflow_uuid` | Originating workflow identifier. |
| `total_rna_counts`, `detected_genes` | RNA depth/detectability covariates. |
| `circRNA_count`, `circRNA_total_support` | Cell-level circRNA detection summaries. |

Optional fixed fields include `condition`, `cell_type`, `batch_id`,
`library_id`, `read_layout`, `strandedness`, `assigned_reads`,
`mitochondrial_fraction`, `ribosomal_fraction`, `alignment_status`, and
`detector_status`.

Wide numeric feature namespaces are supported:

- `host_expr__GENE` for candidate host-gene expression features
- `rbp_program__NAME` for a constrained, documented RBP program
- `splicing_program__NAME` for a constrained splicing program
- `latent__NAME` for optional latent cell state
- `multimodal__NAME` for optional RT, CNV, or other multimodal covariates

Program columns should come from a predeclared feature set. They should not be
used to add an unconstrained, data-dependent vocabulary silently.

## Cell-by-circRNA observations

One row represents one cell/candidate pair.

Required fields:

| Field | Meaning |
| --- | --- |
| `schema_version`, `cell_id`, `circ_id` | Schema and foreign keys. |
| `dataset_id`, `donor_id`, `protocol` | Dataset and cell provenance. |
| `detector`, `reference_genome`, `reference_annotation`, `workflow_uuid` | Detection, reference, and workflow provenance. |
| `count` | Exported pair-level count. |
| `detected` | Whether `count` or `bsj_support` is positive. |
| `bsj_support` | Back-splice junction support. |
| `candidate_label_status` | Copy of candidate `positive`/`unlabelled` status. |
| `detectability_offset` | Precomputed technical detectability term, which may be centered and therefore negative. |

Optional fields are `host_gene_expression`, `label_confidence`,
`count_uncertainty`, `detection_probability`, and `observation_weight`.
These fields reserve stable locations for future count, hurdle, or
positive-unlabelled models without selecting a modelling library now.

## Positive/unlabelled semantics

The v1 schema deliberately has no `negative` label.

- `positive` means the candidate has the positive evidence defined by the
  dataset construction policy.
- `unlabelled` means no positive label has been assigned.
- `detected=false` means no signal was observed for that cell/candidate pair
  under the recorded technical context.
- `novel` means not matched to the selected known-circRNA database.

None of `unlabelled`, `detected=false`, `novel`, or `unknown` is a biological
negative. Future modelling must use positive-unlabelled, censored,
zero-inflated, or other explicitly justified semantics.

## Detectability and technical covariates

Detection opportunity varies with RNA depth, protocol, alignment, detector,
reference, and cell QC. The cell-context table preserves core depth metrics,
while the observation table carries a pair-level `detectability_offset`.

The offset's construction must be recorded by the producing workflow. Suitable
future choices include a log-depth offset, a calibrated detection propensity,
or an exposure term. It must not be interpreted as biological repression
without a separate model and supporting evidence.

## Provenance requirements

The three tables duplicate the critical join provenance intentionally. Bundle
validation requires exact agreement for:

- dataset
- donor
- protocol
- detector
- reference genome
- reference annotation
- workflow UUID
- candidate label status

The JSON sidecar records the `circyto` version, schema versions, source and
output paths, record counts, unique provenance values, and the explicit
`non_detection_is_negative=false` policy.

## Validation and export CLI

Validate and canonicalize three source TSVs:

```bash
circyto export-biogenesis \
  --candidates candidates.tsv \
  --cell-contexts cell_contexts.tsv \
  --observations observations.tsv \
  --outdir biogenesis_bundle
```

Use `--overwrite` only to replace the four canonical bundle files. The command
uses pandas and the Python standard library only. It does not run a detector,
load a genomic language model, or require external executables.

## Intended future modelling use

The three-table split supports future candidate encoders, cell-context
encoders, count/detection heads, donor-aware evaluation, and calibrated
positive-unlabelled learning. Sequence windows and genomic geometry stay on
the candidate axis; expression and technical state stay on the cell axis;
observed support, detectability, and uncertainty stay on the pair axis.

PyTorch, transformers, genomic language models, and two-tower implementations
are intentionally outside this foundation.
