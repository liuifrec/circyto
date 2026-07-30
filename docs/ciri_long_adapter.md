# Optional CIRI-long adapter

## Scientific boundary

This adapter is for Oxford Nanopore circRNA libraries generated with
circRNA enrichment and rolling-circle reverse transcription (RCRT). It is
separate from the generic `circyto nanopore` interoperability workflow.

It does not accept:

- ordinary poly(A)-selected or oligo(dT)-primed ONT cDNA;
- direct RNA;
- unknown library chemistry;
- BAM, SAM, CRAM, or generic minimap2 alignments;
- libraries without circRNA enrichment;
- inputs without explicit reference identity.

Mandalorion `SRR4048177` is not a CIRI-long input. The official CIRI-long demo
and CRR194209 are bulk libraries and do not establish single-cell validation.

## Manifest

The TSV schema version is `circyto.ciri_long_input.v1`.

| Field | Requirement |
| --- | --- |
| `schema_version` | Exactly `circyto.ciri_long_input.v1`. |
| `sample_id` | Unique, path-safe identifier. A sample is not assumed to be a cell. |
| `reads_path` | Raw FASTA/FASTQ, optionally gzip-compressed. |
| `sequencing_platform` | Exactly `OXFORD_NANOPORE`. |
| `molecule_type` | Exactly `cDNA`. |
| `library_preparation` | Exactly `rolling_circle_reverse_transcription`. |
| `circRNA_enrichment` | Must be `true`. |
| `source_accession` | Source run/sample accession; may be empty only when `dataset_id` is populated. |
| `dataset_id` | Dataset identity; may be empty only when `source_accession` is populated. |
| `reference_id` | Explicit reference identity; never inferred from a filename. |
| `reference_build` | Explicit build or assembly release. |
| `biological_interpretation_boundary` | Exactly `full_length_circrna_detection_from_rcrt_library`. |

Additional columns are retained. If `library_selection`, `rna_selection`,
`selection_method`, `priming_method`, or `library_preparation_summary`
documents ordinary poly(A)/oligo(dT) selection, validation fails even when the
fixed chemistry fields were copied incorrectly.

```bash
circyto ciri-long validate-manifest \
  --manifest ciri_long_manifest.tsv \
  --strict
```

## Readiness

The adapter checks:

- the CIRI-long executable and `CIRI-long --version`;
- the BWA executable and version;
- the reference FASTA and SHA-256;
- classic BWA index files `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`;
- optional GTF and known-circRNA annotation paths and SHA-256 values.

```bash
circyto ciri-long doctor \
  --ciri-long CIRI-long \
  --bwa bwa \
  --reference reference.fa \
  --gtf annotation.gtf
```

CIRI-long and its Python/native dependencies remain external and optional.
They are not installed by circyto and are not required by default tests.

## Stages and CLI

`call` and `collapse` are intentionally distinct. Both commands are dry-run by
default and require `--execute` to invoke the external executable.

```bash
# Validate tools and write per-sample argument-array plans
circyto ciri-long plan \
  --manifest ciri_long_manifest.tsv \
  --reference reference.fa \
  --gtf annotation.gtf \
  --outdir work

# Raw reads -> per-sample candidate files
circyto ciri-long call \
  --manifest ciri_long_manifest.tsv \
  --reference reference.fa \
  --gtf annotation.gtf \
  --outdir work \
  --execute

# Candidate FASTA list -> cohort collapse outputs
circyto ciri-long collapse \
  --call-manifest work/ciri_long/call/call_manifest.tsv \
  --reference reference.fa \
  --gtf annotation.gtf \
  --prefix cohort \
  --outdir work \
  --execute

# Official collapse outputs -> separate normalized tables
circyto ciri-long import \
  --collapse-dir work/ciri_long/collapse \
  --prefix cohort \
  --outdir work
```

All external commands use subprocess argument arrays with `shell=False`.
Provenance includes the exact argv and an escaped display-only rendering.

## Outputs

```text
work/ciri_long/
  call/
    call_manifest.tsv
    call_summary.json
    SAMPLE_ID/
      SAMPLE_ID.cand_circ.fa
      SAMPLE_ID.low_confidence.fa
      SAMPLE_ID.json
      SAMPLE_ID.log
      wrapper.stdout.log
      wrapper.stderr.log
      provenance.json
      tmp/
  collapse/
    calls.list
    cohort.info
    cohort.expression
    cohort.isoforms
    cohort.reads
    cohort.log
    wrapper.stdout.log
    wrapper.stderr.log
    provenance.json
  normalized/
    circRNA_bsj.tsv
    circRNA_isoforms.tsv
    circRNA_expression.tsv
    circRNA_isoform_usage.tsv
    circRNA_read_assignments.tsv
    import_summary.json
```

The adapter does not write AnnData/MuData and does not modify circyto's
validated `circ` modality.

## Official-output mapping

| Official source | Normalized output | Preserved information |
| --- | --- | --- |
| `cohort.info` | `circRNA_bsj.tsv` | Official circRNA ID, contig, BSJ coordinates, strand, support score, circRNA type, splice-site/equivalent-sequence attributes, host-gene fields, full original attributes and source line. |
| `cohort.info` `isoform` attribute | `circRNA_isoforms.tsv` | Parent circRNA, official-style isoform ID, ordered exon blocks, calculated transcript length, reported-major ordering, host-gene fields and original structure. |
| `cohort.expression` | `circRNA_expression.tsv` | Long-form circRNA, sample and support values plus the original value. |
| `cohort.isoforms` | `circRNA_isoform_usage.tsv` | Long-form isoform, parent circRNA, sample and usage. Usage is not reinterpreted as an absolute count. |
| `cohort.reads` | `circRNA_read_assignments.tsv` | Every official column, source cohort and an exact JSON representation of the original row. |

BSJ-level and isoform-level records remain separate. Bulk sample identifiers
are never relabelled as cell identifiers.

## Coordinate convention

CIRI-long `.info` positions and exon-block strings are treated as
**CIRI-series 1-based closed** coordinates. Normalized tables use
**0-based half-open** coordinates:

```text
normalized_start = original_start - 1
normalized_end   = original_end
```

For BSJs, normalized values are written in `start` and `end`; the official
values remain in `original_start` and `original_end`. Each row records:

- `original_coordinate_system`;
- `normalized_coordinate_system`;
- `coordinate_conversion_rule`.

Isoform `exon_block_structure` is JSON containing normalized `[start, end]`
blocks. `original_exon_block_structure` retains the official text.
`transcript_length` is calculated from the 1-based closed blocks as
`sum(end - start + 1)`.

The first official isoform structure is marked `reported_major`; this records
CIRI-long's reported order and major-length association, not an independently
re-estimated abundance rank.

## Provenance

Call and collapse sidecars record:

- CIRI-long path, version, exact argv and return code;
- BWA path and version;
- reference path, ID, build, SHA-256 and all BWA index assets;
- GTF and optional known-circRNA annotation paths and SHA-256 values;
- input read/candidate paths and checksums;
- the accepted chemistry gate and rejected assumptions;
- source accessions/dataset IDs;
- official output paths, wrapper logs, timestamps and warnings.

`import_summary.json` records all official and normalized paths, checksums,
record counts, coordinate policy, entity separation and scientific warnings.

## Official demo

The opt-in script queries the official GitHub release before transfer, verifies
the current asset identity and compressed size, and verifies an official
SHA-256 if the release API supplies one. Otherwise it records a local SHA-256
without describing it as an archive-provided checksum.

```bash
# No network or filesystem writes
python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_demo \
  --dry-run

# Query release metadata only
python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_demo \
  --metadata-only

# Download and extract, but do not run external tools
python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_demo

# Opt-in end-to-end detector smoke
python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_demo \
  --execute
```

The archive is downloaded on demand and is never redistributed by circyto.
CRR194209 is not downloaded by this script.
