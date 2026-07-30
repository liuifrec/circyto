# CIRI-long compatibility assessment

Status: architecture/source review completed 2026-07-29; optional adapter
implemented for issue #14 on 2026-07-30. The implementation is covered by
hermetic tests, but the official demo has not yet been downloaded or run.
Operational details are in [ciri_long_adapter.md](ciri_long_adapter.md).

## Decision

CIRI-long belongs in a separate, chemistry-gated circRNA detection path. It is
not an alternative aligner for ordinary Oxford Nanopore cDNA and must not
consume the BAMs produced by circyto's generic minimap2 interoperability path.

The optional adapter is appropriate only when source metadata establishes that
the library used rolling-circle reverse transcription (RCRT) of enriched
circular RNA and therefore contains concatemeric reads with repeated copies of
a circular template. A FASTA/FASTQ filename alone is not evidence of this
compatibility.

The official CIRI-long workflow is bulk rather than single-cell. Neither the
official mouse-brain demonstration nor its protocol establishes single-cell
support. circyto must preserve `single_cell_status` and may describe a result
as single-cell only when the source experiment genuinely resolves cells.

Primary sources:

- [CIRI-long software and MIT license](https://github.com/bioinfo-biols/CIRI-long)
- [official usage and output documentation](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-long_2_usage.html)
- [official installation requirements](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-long_1_installation.html)
- [CIRI-long method publication](https://www.nature.com/articles/s41587-021-00842-6)
- [CIRI-long experimental protocol](https://www.nature.com/articles/s41596-023-00815-w)

## Required input and read structure

`CIRI-long call` consumes raw long reads, not a generic genomic alignment:

```text
CIRI-long call \
  -i reads.fastq.gz \
  -o output_directory \
  -r reference.fa \
  -a annotation.gtf \
  -p sample_prefix
```

The implementation accepts FASTA or FASTQ, optionally gzip-compressed. It
also requires an indexed reference genome; the official instructions build a
BWA index for the FASTA. An optional known-circRNA file can be supplied with
`-c`.

The expected molecule is an RCRT product in which reverse transcriptase has
traversed a circular template repeatedly. The resulting long cDNA read
contains tandem copies of the circle. CIRI-long detects repeated segments,
constructs a consensus sequence, and corrects the candidate against the
reference and splice annotation. Its intermediate `*.ccs.fa` and `*.raw.fa`
files reflect this repeat/consensus step.

The protocol uses RNase R treatment, enzymatic poly(A) tailing, RCRT, PCR, and
size selection to enrich these molecules. Enzymatic tailing during a
circRNA-specific protocol is not conventional poly(A)-selection of linear
transcripts and must be represented separately in provenance.

Ordinary ONT cDNA, including Mandalorion, SiCeLoRe, ScNapBar, and other 10x- or
Smart-seq-derived cDNA, generally lacks the required tandem-repeat structure.
Although the parser can read its FASTQ, that does not make the library a
scientifically valid CIRI-long input. Those datasets remain in the generic
interoperability path.

## External dependencies

The official installation documentation and package metadata require:

- Python 3.7 or newer;
- GCC 4.8 or newer, or Clang 3.4 or newer;
- CMake 3.2 or newer;
- samtools 1.9 or newer;
- a BWA-indexed reference;
- Python packages including Cython, mappy, numpy, pandas, pysam,
  python-Levenshtein, scipy, scikit-learn, Biopython, pyspoa, bwapy, edlib,
  and pyccs;
- native compression headers needed to build pysam, where binary wheels are
  not available;
- the bundled striped Smith-Waterman native component built during package
  installation.

This dependency stack must remain optional. It must not be added to circyto's
default installation, default tests, or generic Nanopore path. A future
adapter should preflight the external `CIRI-long` executable and record its
reported version rather than importing its dependency stack into circyto.

## Official outputs

### `call`

The documented `call` outputs are:

| Output | Meaning and adapter treatment |
| --- | --- |
| `PREFIX.cand_circ.fa` | Candidate circular consensus sequences. This is the input to `collapse`, not the final normalized call table. |
| `PREFIX.json` | Per-run machine-readable call metadata. Preserve as a detector artifact. |
| `PREFIX.log` | Execution log. Preserve in provenance. |
| `PREFIX.low_confidence.fa` | Low-confidence sequences. Preserve separately; do not promote to validated calls. |
| `tmp/PREFIX.ccs.fa` | Reconstructed circular consensus sequences. Intermediate evidence. |
| `tmp/PREFIX.raw.fa` | Associated raw read sequences. Intermediate evidence. |
| `tmp/ss.idx` | Internal splice-site index. Implementation artifact. |

### `collapse`

`CIRI-long collapse` takes a two-column list of sample names and paths to
`*.cand_circ.fa` files. Its final outputs include:

| Output | Documented fields and future normalization |
| --- | --- |
| `PREFIX.info` | GTF-like rows containing chromosome, source, feature type, BSJ start/end, support score, strand, and attributes including `circ_id`, splice sites, equivalent sequence, circRNA type, length, isoform exon structures, and gene identifiers. Normalize to circyto circRNA candidate and full-length isoform records. |
| `PREFIX.expression` | BSJ-level sample expression table indexed by `circ_ID`. Preserve the counting unit and sample mapping. |
| `PREFIX.isoforms` | Isoform-usage table indexed by `isoform_ID`; usage is relative to reads supporting the same BSJ. Do not reinterpret usage as an absolute count. |
| `PREFIX.reads` | Per-read assignments containing read, circRNA, temporary isoform, strand, circular-exon, signal, alignment, segment, sample, and type fields. Preserve as detector evidence. |
| `PREFIX.log` | Collapse log. Preserve in provenance. |
| `tmp/ss.idx`, `tmp/corrected.pkl` | Internal splice-site and correction artifacts. |

The future normalizer should emit two related, versioned tables:

1. circRNA records: canonical BSJ coordinates and strand, `circ_id`,
   splice-site and equivalent-sequence annotations, circRNA type and length,
   gene identifiers, sample-level support, detector confidence class, and
   complete reference/protocol/detector provenance;
2. full-length isoform records: `circ_id`, stable `isoform_id`, ordered
   circular-exon blocks, reported usage and its denominator, derived support
   count only when justified by the official evidence, sample identity, and
   the same provenance keys.

No CIRI-long table should be normalized into the established single-cell
`circ` modality merely because it contains a sample column. Bulk samples and
cells are different observation units.

## Public datasets and practical sizes

| Candidate | Role | Transfer size | Eligibility and limitations |
| --- | --- | ---: | --- |
| Official `v0.6-alpha` demonstration archive | First adapter contract test | 36,945,437 compressed bytes (about 35.2 MiB), per the GitHub release API on 2026-07-29 | Small, official, and self-contained with test reads, a chromosome-12 reference, and annotation. It validates invocation and parsing, not independent biological performance. It is tied to an old release and has no published asset digest or separately stated data license. |
| CRR194209, CRA003317 | First real-library follow-up | 1,293,153,736 compressed bytes (about 1.20 GiB), from the official GSA object on 2026-07-29 | Adult mouse-brain RCRT/CIRI-long data associated with the method and protocol. Bulk, not single-cell. Requires a separately obtained compatible mouse reference, annotation, and BWA index. |
| GSE197872 / PRJNA812557, for example SRR18213999 | Conditional, not currently recommended | SRR18213999 is approximately 0.71 GB in ENA metadata | The study reports CIRI-long processing of circRNA-enriched ONT libraries, but the archive description reviewed here does not by itself establish the required RCRT concatemer chemistry. It fails the strict Path B admission gate until the exact library protocol is verified. It is bulk cell-line data, not single-cell. |

### Recommendation

Use the **official 36,945,437-byte demonstration archive** as the first future
CIRI-long adapter test. It is the smallest practical way to verify executable
discovery, raw-read invocation, output collection, schema parsing, and
provenance without a large transfer. Fetch it on demand and record the
observed checksum; do not vendor it.

After that contract test passes, use **CRR194209** as the first real
RCRT-library validation exercise. Its size should be queried again before
download and its archive checksum verified. It supports detector-level
reproducibility, not a single-cell claim.

Mandalorion `SRR4048177` is explicitly excluded from both tests: it remains
the small generic long-read interoperability smoke dataset.

## Licensing and redistribution

The CIRI-long source repository declares the MIT license, which covers the
software. It does not automatically license sequencing reads, bundled
reference sequences, annotations, or third-party databases.

The official demo release does not publish a cryptographic digest or a
separate redistribution license for its data assets. circyto should therefore
store only a download plan, expected filename and byte-size warning, source
URL, retrieval timestamp, and locally calculated checksum. The archive itself
must not be committed or repackaged until redistribution rights are
established.

CRR194209 and any other public reads should be retrieved from their archive
under the archive's data-use terms. Reference FASTA and GTF files retain their
own provider licenses and likewise should not be bundled by default.

## Adapter contract

The CIRI-long path is a detector adapter, not an extension of
the minimap2 alignment command:

```text
RCRT/enriched raw FASTA or FASTQ
        |
        v
chemistry-gated circRNA-library manifest
        |
        v
CIRI-long call
        |
        v
sample cand_circ.fa files
        |
        v
CIRI-long collapse
        |
        +--> normalized circRNA records
        +--> normalized full-length isoform records
        `--> retained per-read detector evidence
```

Minimum manifest/provenance fields are:

- dataset, sample, run, and source accessions;
- biological observation unit (`bulk_sample`, `single_cell`, or another
  explicitly defined unit);
- organism and tissue/cell line;
- circRNA enrichment steps;
- `rolling_circle_rt_confirmed` plus a citable protocol source;
- reverse transcriptase, PCR, size-selection, and sequencing protocol when
  known;
- raw-read path, checksum, and read count;
- reference ID, build, FASTA checksum, annotation release and checksum;
- CIRI-long version, exact `call` and `collapse` commands, optional known-circ
  source, exit status, and output checksums;
- detector confidence class and sample/read counting units.

Admission must fail clearly when RCRT is false, unknown, or unsupported; when
raw FASTA/FASTQ is absent; or when reference identity is incomplete. An
override must not silently convert an ineligible ordinary cDNA library into a
supported detector input.

Default tests should use synthetic tandem-repeat reads and a fake CIRI-long
executable. The official demo and CRR194209 must remain opt-in tests requiring
network access and an independently installed CIRI-long environment.

## Scientific wording

Before detector and artifact validation:

> circyto provides a planned adapter boundary for CIRI-long-compatible
> rolling-circle reverse-transcription libraries. This path is separate from
> generic Nanopore cDNA alignment and is not exercised by the single-cell
> Mandalorion interoperability example.

After reproducible execution on the official demo and CRR194209:

> We demonstrated experimental interoperability with CIRI-long using public
> bulk Oxford Nanopore libraries generated with circRNA enrichment and
> rolling-circle reverse transcription. Detector outputs were normalized with
> protocol, reference, and software provenance. This result does not establish
> single-cell circRNA support.

Neither statement is biological validation of individual circRNA calls.
Claims about sensitivity, specificity, or single-cell biology require
additional truth sets and orthogonal validation.
