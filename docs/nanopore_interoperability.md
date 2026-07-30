# Path A: experimental single-cell Nanopore interoperability

## Scientific boundary

This implementation is an experimental long-read ingestion and alignment
demonstration. It is not biological validation of circRNA detection.

> **EXPERIMENTAL LONG-READ INTEROPERABILITY ONLY.** This workflow does not
> validate circRNA detection. SRR4048177 was generated using a poly(A)-oriented
> cDNA protocol, so conventional non-polyadenylated circRNAs are expected to be
> depleted. Non-detection must not be interpreted as biological absence.

No CIRI3 or other short-read circRNA detector is invoked. The workflow does not
write AnnData or MuData, does not modify the validated `circ` modality, and
does not promote exploratory long-read evidence to circRNA calls.

This is the generic long-read path for Mandalorion, SiCeLoRe, ScNapBar, and
ordinary ONT cDNA. Its minimap2 BAM is not input to the separate future
CIRI-long detector path. CIRI-long requires raw reads from a compatible
rolling-circle reverse-transcription circRNA library; see
[the CIRI-long compatibility assessment](ciri_long_feasibility.md).

## Archive identity

The opt-in smoke run is fixed to Mandalorion run `SRR4048177` from
`PRJNA339767`. Before any download, circyto queries the ENA run-report API and
compares the result with
`circyto/resources/nanopore/srr4048177_expected.json`.

Hard failures cover:

- run and project accessions;
- Oxford Nanopore platform and instrument;
- archive library strategy, source, and selection;
- the expected compressed FASTQ filename;
- the archive MD5.

URL, compressed bytes, read count, and base count are warning-level
comparisons because archives can revise these representations without changing
run identity. The current URL and size are always taken from the fresh ENA
response. The checked-in compressed-byte value is only a dated comparison
baseline and is never used as a download constant.

ENA currently reports `library_selection=RACE`. This value is retained as
`archive_library_selection`; it is not reinterpreted as poly(A). The separate
`library_preparation_summary` records the paper-derived description of the
modified Smart-seq2-style, oligo(dT)-primed cDNA preparation.

Query metadata without downloading reads:

```bash
circyto nanopore query-run \
  --accession SRR4048177 \
  --expected srr4048177 \
  --output archive_metadata.json
```

## Long-read manifest

The TSV schema version is `circyto.long_read_single_cell.v1`. Every row
requires:

- `schema_version`
- `cell_id`
- `long_read_fastq`
- `protocol`
- `sequencing_platform`
- `archive_library_selection`
- `library_preparation_summary`
- `molecule_type`
- `barcode_status`
- `source_accession`
- `dataset_id`
- `biological_interpretation_boundary`

`molecule_type` is either `cdna` or `direct_rna`. The interpretation boundary
must be:

```text
experimental_long_read_interoperability_not_circrna_validation
```

Cell identifiers are restricted to path-safe characters so two logical cells
cannot collapse onto the same output directory.

```bash
circyto nanopore validate-manifest \
  --manifest long_read_manifest.tsv \
  --strict
```

## Alignment

Reference identity is never inferred from a filename. The caller must supply a
reference ID, reference build, and the expected SHA-256 of the exact FASTA.
circyto recomputes and verifies the digest before starting minimap2.

For cDNA, the alignment process uses the official splice-aware form:

```text
minimap2 -ax splice -t THREADS reference.fa reads.fastq
```

For an input explicitly declared as `direct_rna`, circyto adds `-uf -k14`.
Those molecule-specific options cannot be injected through generic extra
arguments.

minimap2 stdout is connected directly to `samtools sort` using subprocess
argument arrays with `shell=False`. A full SAM is not retained by default.
`--keep-sam` is available only for debugging and writes the SAM obtained from
the sorted BAM during QC inspection.

```bash
circyto nanopore align \
  --manifest long_read_manifest.tsv \
  --reference reference.fa \
  --reference-id EXPLICIT_REFERENCE_ID \
  --reference-build EXPLICIT_REFERENCE_BUILD \
  --reference-sha256 EXPECTED_FASTA_SHA256 \
  --archive-metadata archive_metadata.json \
  --outdir work/srr4048177 \
  --threads 8
```

Artifacts are cell-scoped:

```text
work/srr4048177/nanopore_alignment/
  alignment_manifest.tsv
  nanopore_alignment_summary.json
  WT_B1a_Cell1/
    alignment.bam
    alignment.bam.bai
    alignment_qc.json
    exploratory_bsj_evidence.tsv
    provenance.json
    minimap2.stderr.log
    samtools_sort.stderr.log
    samtools_index.log
```

Each alignment-manifest row points to its own BAM, index, QC, evidence, and
provenance files.

## QC definitions

QC reports record-level and query-level quantities separately:

- `input_query_count`: complete records in the input FASTQ;
- `mapped_primary_query_count`: input query names with a mapped primary
  alignment;
- `unmapped_query_count`: input queries without a mapped primary alignment;
- `secondary_alignment_record_count`: SAM records carrying flag `0x100`;
- `supplementary_alignment_record_count`: SAM records carrying flag `0x800`;
- `queries_with_supplementary_alignments`: query names represented by at least
  one supplementary record;
- `queries_with_sa_tag`: query names with at least one `SA` tag;
- `spliced_primary_query_count`: mapped primary query names whose primary CIGAR
  contains `N`.

`mapped_primary_query_fraction` uses `input_query_count` as its denominator.
`spliced_primary_query_fraction` uses `mapped_primary_query_count`. An ordinary
`N` operation is a linear splice and is never sufficient exploratory
back-splice evidence.

## Exploratory alignment patterns

The diagnostic requires a mapped primary segment and a mapped supplementary
segment on the same contig and strand. Secondary-only combinations are
discarded. Query intervals are normalized to the original query orientation,
reference order is normalized by strand, and the two orders must be inverted.
Overlapping query segments are rejected.

Every retained row includes the read name, flags, mapping qualities, CIGARs,
query and reference intervals, strand, all supporting segments as JSON, and
reason code
`query_reference_order_inversion_same_contig_strand`.
`circRNA_call` is always `false`.

> `exploratory_bsj_evidence` contains alignment patterns compatible with a
> possible back-splice. Entries are not circRNA calls. PCR and template-switch
> artifacts, chimeric cDNA, basecalling errors, and long-read alignment errors
> can produce similar patterns.

## Opt-in real-data smoke run

The script below is never invoked by default tests or CI and has no accession
option: it can download only `SRR4048177`.

```bash
python scripts/run_srr4048177_nanopore_smoke.py \
  --outdir work/srr4048177 \
  --reference reference.fa \
  --reference-id EXPLICIT_REFERENCE_ID \
  --reference-build EXPLICIT_REFERENCE_BUILD \
  --reference-sha256 EXPECTED_FASTA_SHA256 \
  --threads 8
```

Use `--metadata-only` to stop before download or `--download-only` to stop
before alignment. Downloads use a `.part` file, timeout-bounded requests,
remote identity sidecar data, HTTP Range when safe, and MD5 verification
before final promotion. If a server ignores Range, the temporary file is
rewritten from byte zero. Oversized or stale partial files are discarded;
failed MD5 content is retained with a `.md5_failed` suffix and is never
promoted to the final FASTQ path.
