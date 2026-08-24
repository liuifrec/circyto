# SRR4048177 Nanopore interoperability validation

## Decision and scientific boundary

**Validation date:** 2026-07-31 (Asia/Tokyo; archive and workflow timestamps are
recorded in UTC)

**Decision:** PASS on the validation worktree after correcting a provenance
sidecar omission described below. The archive identity and checksum gates,
manifest validation, cDNA splice-aware alignment, BAM sorting and indexing,
query-level QC, cell-scoped output isolation, and conservative exploratory
alignment diagnostics all passed. The code correction and this report require
review before commit.

This validation demonstrates provenance-aware ingestion and splice-aware
alignment of a physically isolated single-cell Oxford Nanopore cDNA dataset. It
validates circyto's experimental long-read interoperability workflow, including
archive verification, manifest validation, cell-scoped output generation,
alignment quality control, and conservative split-alignment diagnostics.
Because the library was generated using a poly(A)-oriented cDNA protocol, this
analysis does not validate circRNA detection, and circRNA non-detection cannot
be interpreted as biological absence.

No CIRI-long, CIRI3, or other circRNA detector was invoked. No validated circ
modality or AnnData/MuData object was created or modified.

## Repository and execution environment

The validation began from a clean `feature/nanopore-feasibility` branch at
`cf02922`, above the generic interoperability commit `a9ab739`, the optional
CIRI-long adapter commit `12f5478`, and the committed official CIRI-long
validation report. The branch matched
`origin/feature/nanopore-feasibility`.

| Component | Value |
| --- | --- |
| Host | validation workstation |
| OS | Ubuntu 22.04.5 LTS; Linux 6.8.0-136-generic x86_64 |
| CPU | Intel Core Ultra 9 185H; 22 online logical CPUs |
| circyto | 0.10.0 |
| circyto Python | 3.10.20 |
| circyto executable | `$CIRCYTO_ENV/bin/circyto` |
| minimap2 | 2.31-r1302 |
| minimap2 executable | `$NANOPORE_ENV/bin/minimap2` |
| samtools/HTSlib | 1.24/1.24 |
| samtools executable | `$NANOPORE_ENV/bin/samtools` |
| External-tool environment | isolated micromamba environment `nanopore-interop` (`$NANOPORE_ENV`) |

minimap2 and samtools were kept outside the circyto environment. All sequence,
reference, alignment, index, log, and temporary files were placed under the
Git-ignored `work/srr4048177_nanopore/` directory.

## Live ENA metadata and download integrity

Metadata were retrieved at `2026-07-31T00:54:23+00:00` through the ENA Portal
API URL recorded verbatim in `archive_metadata.json`. The hard identity check
against `circyto/resources/nanopore/srr4048177_expected.json` passed. There were
no warning-level URL, byte-count, read-count, or base-count differences.

| Field | Live value |
| --- | --- |
| Run | SRR4048177 |
| Study/project | PRJNA339767 |
| Sample | SAMN05607017 |
| Experiment | SRX2039007 |
| Organism | *Mus musculus* |
| Platform/instrument | OXFORD_NANOPORE / MinION |
| Library strategy | RNA-Seq |
| Library source | TRANSCRIPTOMIC |
| Archive library selection | RACE |
| Bases | 70,823,600 |
| Reads | 52,696 |
| FASTQ object | `SRR4048177_1.fastq.gz` |
| Compressed bytes | 50,513,678 |
| FASTQ URL | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR404/007/SRR4048177/SRR4048177_1.fastq.gz` |
| Archive MD5 | `7e841c27e8afd6c2303f06e7a470d4a6` |

The remote object advertised byte ranges, a content length of 50,513,678, a
last-modified date of 2019-05-11, and ETag
`"302c70e-58895bc0c9880"`. Before transfer, the destination was confirmed to
match `.gitignore`'s `work*/` rule and 715+ GB was available.

The opt-in downloader used its `.part`/remote-identity mechanism and promoted
the object only after MD5 verification. This cold transfer resumed from zero
bytes and the server did not ignore a range request. The final file was
50,513,678 bytes, its MD5 matched ENA, and its local SHA-256 was
`09b35d9511b756b1dd5fc0d6597b7edd95999fcaceee5db7f190d4230e774525`.
Streaming decompression counted exactly 52,696 four-line FASTQ records; the
source archive was not modified or recompressed.

Archive `library_selection=RACE` is preserved as archive metadata and is not
reinterpreted as poly(A). Separately, the publication-derived preparation
summary is: "Modified Smart-seq2-style, oligo(dT)-primed single-cell cDNA
sequenced with Oxford Nanopore 2D reads."

## Reference

The Mandalorion publication reports GRCm38 for mouse alignment and GENCODE vM10
for annotation. This run used the official GENCODE mouse release M10 GRCm38
primary assembly, with identity supplied explicitly rather than inferred from
the filename.

| Field | Value |
| --- | --- |
| Reference ID | `GENCODE_GRCm38_primary_assembly_release_M10` |
| Reference build | `GRCm38` |
| Source | GENCODE mouse release M10 |
| Source URL | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/GRCm38.primary_assembly.genome.fa.gz` |
| Compressed bytes | 773,507,376 |
| Upstream compressed MD5 | `7e386bbdc10ebfb644328668fb413719` |
| Compressed SHA-256 | `09de3c56747c69948d5ca39125caa5bebeac41e0de2e4a67bac947bfc80864f0` |
| Uncompressed bytes | 2,776,387,546 |
| Uncompressed FASTA SHA-256 | `767f48e8ae5bda117ad519fb694cbf7e24f161bf094c274bbbcd16e52aedf802` |
| Contigs | 66; `chr1`-`chr19`, `chrX`, `chrY`, `chrM`, and GRCm38 unplaced/unlocalized sequences |
| Annotation in this workflow | None |

The generic alignment interface does not consume an annotation. No GTF was
downloaded or used, so an annotation SHA-256 is not applicable. GENCODE vM10 is
recorded as study context, not as an unverified alignment input.

## Manifest and pre-execution gates

The generated one-row
`circyto.long_read_single_cell.v1` manifest passed strict validation. Its
identity and interpretation fields were:

| Field | Value |
| --- | --- |
| `cell_id` | `WT_B1a_Cell1` |
| Physical context | one physically isolated mouse B1a cell |
| `protocol` | `mandalorion_modified_smartseq2` |
| `sequencing_platform` | `OXFORD_NANOPORE` |
| `archive_library_selection` | `RACE` |
| `molecule_type` | `cdna` |
| `barcode_status` | `not_applicable_physical_single_cell` |
| `source_accession` | `SRR4048177` |
| `dataset_id` | `PRJNA339767` |
| Interpretation boundary | `experimental_long_read_interoperability_not_circrna_validation` |

The manifest SHA-256 was
`baaccafc2b4b0f629f3a68e89262deef2f1eaef95585415d1bbcf83f973f334a`.
The reference ID/build are required alignment arguments and are captured in
alignment provenance rather than guessed from the manifest's FASTA name.

The implemented smoke script provides `--metadata-only` and `--download-only`;
it does not provide separate download-plan, manifest-only dry-run, or alignment
plan flags. Therefore no nonexistent mode was claimed. Before download,
metadata-only execution, a remote HEAD request, Git-ignore verification, and a
disk-space check supplied the non-mutating gates. Before alignment, the actual
argument builder was invoked without execution to inspect the minimap2 and
samtools argv arrays. `--download-only` then performed the checked transfer and
manifest generation without alignment.

## Commands

The substantive commands were:

```text
circyto nanopore query-run \
  --accession SRR4048177 \
  --expected srr4048177 \
  --output work/srr4048177_nanopore/archive_metadata_query.json \
  --timeout-seconds 60

python scripts/run_srr4048177_nanopore_smoke.py \
  --outdir work/srr4048177_nanopore \
  --metadata-only \
  --timeout-seconds 60

python scripts/run_srr4048177_nanopore_smoke.py \
  --outdir work/srr4048177_nanopore \
  --download-only \
  --timeout-seconds 120

circyto nanopore validate-manifest \
  --manifest work/srr4048177_nanopore/long_read_manifest.tsv \
  --strict

circyto nanopore align \
  --manifest work/srr4048177_nanopore/long_read_manifest.tsv \
  --reference work/srr4048177_nanopore/reference/GRCm38.primary_assembly.genome.fa \
  --reference-id GENCODE_GRCm38_primary_assembly_release_M10 \
  --reference-build GRCm38 \
  --reference-sha256 767f48e8ae5bda117ad519fb694cbf7e24f161bf094c274bbbcd16e52aedf802 \
  --outdir work/srr4048177_nanopore/validated_run \
  --threads 8 \
  --minimap2 "$NANOPORE_ENV/bin/minimap2" \
  --samtools "$NANOPORE_ENV/bin/samtools" \
  --archive-metadata work/srr4048177_nanopore/archive_metadata.json

circyto nanopore inspect-bsj \
  --alignment work/srr4048177_nanopore/nanopore_alignment/WT_B1a_Cell1/alignment.bam \
  --cell-id WT_B1a_Cell1 \
  --input-query-count 52696 \
  --output work/srr4048177_nanopore/diagnostic_review/exploratory_bsj_evidence.tsv \
  --qc-output work/srr4048177_nanopore/diagnostic_review/exploratory_bsj_qc.json \
  --samtools "$NANOPORE_ENV/bin/samtools"
```

The argv captured in provenance, shown with workstation paths normalized, were:

```text
$NANOPORE_ENV/bin/minimap2 -ax splice -t 8 work/srr4048177_nanopore/reference/GRCm38.primary_assembly.genome.fa work/srr4048177_nanopore/downloads/SRR4048177_1.fastq.gz

$NANOPORE_ENV/bin/samtools sort -@ 8 -o work/srr4048177_nanopore/validated_run/nanopore_alignment/WT_B1a_Cell1/alignment.partial.bam -

$NANOPORE_ENV/bin/samtools index work/srr4048177_nanopore/validated_run/nanopore_alignment/WT_B1a_Cell1/alignment.bam
```

All subprocesses used argument arrays and `shell=False`. minimap2 stdout was
streamed into `samtools sort`; no full SAM was retained. The cDNA command used
`-ax splice` and did not use the direct-RNA-only `-uf -k14` options.

## Alignment results

The confirmation run exited 0 in 1 minute 13.34 seconds wall time. minimap2
reported mapping all 52,696 sequences in 63.051 seconds, using 16.999 GB peak
RSS. The coordinate-sorted BAM and BAI were created under the cell-scoped
`WT_B1a_Cell1/` directory and passed `samtools quickcheck`.

| Metric | Value |
| --- | ---: |
| `input_query_count` | 52,696 |
| `alignment_record_count` | 80,102 |
| `mapped_primary_query_count` | 39,537 |
| `unmapped_query_count` | 13,159 |
| `secondary_alignment_record_count` | 26,867 |
| `supplementary_alignment_record_count` | 539 |
| `queries_with_supplementary_alignments` | 537 |
| `queries_with_sa_tag` | 537 |
| `spliced_primary_query_count` | 33,270 |
| mapped primary / input | 0.750285 |
| unmapped / input | 0.249715 |
| spliced primary / mapped primary | 0.841490 |

These are explicit query-level and alignment-record-level metrics. No
MAPQ-positive category is described as "unique."

The first and post-fix confirmation runs produced identical QC JSON and
byte-identical exploratory evidence. Small BAM byte-size variation reflects
run-specific BAM header metadata rather than a result difference.

## Exploratory alignment diagnostics

`inspect-bsj` reported five diagnostic patterns across five reads:

| Property | Result |
| --- | --- |
| Reason code | 5 × `query_reference_order_inversion_same_contig_strand` |
| Contigs | 4 × chr6; 1 × chr18 |
| Normalized strand | 3 × `+`; 2 × `-` |
| Segment MAPQ pairs | 5 × (60, 60) |
| `circRNA_call` | `false` for all five |

The diagnostic requires a mapped primary plus a non-secondary supplementary
segment on the same contig and normalized strand, non-overlapping query order,
and inverted reference order. Consequently, an ordinary N-CIGAR junction alone
is insufficient and secondary-only combinations are rejected. The TSV retains
read names, flags, MAPQ values, CIGARs, query/reference intervals, normalized
strand, all supporting segment details, and a reason code.

Four patterns share a narrow chr6 region and similar segment endpoints; the
fifth connects chr18 segments approximately 6.59 Mb apart. This recurrence and
the extensive clipping/indel structure are compatible with PCR amplification,
template switching, chimeric cDNA, basecalling error, or long-read alignment
error. The five patterns are not circRNA calls, do not measure circRNA
abundance, and do not establish biological accuracy.

## Provenance audit and defect found

The initial real run exposed a provenance-only defect: the adapter wrote a
cell provenance sidecar but omitted the three workflow log paths and its own
path, and it did not emit the required run-level provenance sidecar. Alignment,
QC, and exploratory evidence were unaffected.

The focused correction:

- records BAM, BAI, QC, exploratory evidence, retained-SAM status, all three
  tool logs, and cell provenance in each cell's `paths`;
- emits `nanopore_alignment/provenance.json` using schema
  `circyto.nanopore_alignment_run_provenance.v1`;
- records the root alignment manifest, summary, run provenance, every cell's
  output paths, exact command arrays, QC, archive snapshot, reference, tool
  identity, stage status, and scientific disclaimers;
- keeps `detector_invoked=false`, `detector_backend=null`, and
  `circRNA_validation_status=false`.

The post-fix real run confirmed that every non-null path named by both
provenance levels exists. It also confirmed:

| Audit item | Status |
| --- | --- |
| ENA snapshot, retrieval URL/time | PASS |
| Source accession/project | PASS |
| FASTQ archive MD5 and local SHA-256 | PASS |
| Explicit reference ID/build/path/SHA-256 | PASS |
| minimap2/samtools paths and versions | PASS |
| Exact mapper/sort/index/view argv | PASS |
| Manifest path and SHA-256 | PASS |
| Complete QC metrics | PASS |
| All workflow-generated output paths | PASS after focused fix |
| `detector_invoked=false` | PASS |
| `circRNA_validation_status=false` | PASS |
| Poly(A)-oriented depletion warning | PASS |
| Non-detection-is-not-absence warning | PASS |
| CIRI-long execution | Not invoked |
| CIRI3/short-read detector execution | Not invoked |
| Validated circ modality | Not created |
| AnnData/MuData behavior | Not changed |
| Generated data tracked by Git | None |

## Limitations

- SRR4048177 represents one physically isolated cell and cannot establish
  performance across cells, tissues, library batches, or ONT chemistries.
- This legacy 2D cDNA dataset is poly(A)-oriented. Conventional
  non-polyadenylated circRNAs are expected to be depleted.
- The alignment uses a current minimap2 release rather than the publication's
  historical BLAT workflow; this is intentionally an interoperability test,
  not a reproduction of published transcript models.
- GENCODE vM10 annotation was not consumed by this annotation-free alignment
  path.
- Split/supplementary alignment patterns are vulnerable to PCR,
  template-switch, chimeric-cDNA, basecalling, and alignment artifacts.
- The smoke script has metadata-only and download-only modes but no dedicated
  plan-only CLI. The missing planning modes were handled transparently with
  metadata/HEAD checks and direct argv inspection rather than simulated.
- No statement about circRNA detection sensitivity, specificity, abundance, or
  biological absence is supported by this validation.
