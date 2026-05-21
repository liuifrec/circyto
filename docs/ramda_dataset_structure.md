# RamDA/scRR Dataset Structure

This note clarifies whether the current public human scRR / RamDA-like datasets should be processed as direct per-cell runs or as pooled multiplexed libraries that require demultiplexing.

## Short answer

Both currently targeted public human datasets are best treated as direct one-library-per-cell runs.

Classification:

- `GSE278958` IMR90: direct single-cell-per-run
- `GSE278952` HAP1: direct single-cell-per-run

Demultiplexing is not required for either dataset in the current `circyto` workflow design.

## Evidence from repo-tracked public metadata

The current repo snapshots resolve the following RNA-side SRRs:

### `GSE278958` IMR90

- `SRR30918126`
- `SRR30918117`

Properties recorded in repo metadata:

- RNA-side only
- human
- single-end
- one SRR per RNA library / cell candidate
- recommended direct execution with `circyto workflow full-length-circrna`

Estimated RNA cell count from the current validated public set captured in the repo:

- 2 RNA-side cell libraries

### `GSE278952` HAP1

- `SRR30911454`
- `SRR30911453`
- `SRR30911559`

Properties recorded in repo metadata:

- RNA-side only
- human
- paired-end
- one SRR per RNA library / cell candidate
- recommended direct execution with `circyto workflow full-length-circrna --star-index ... --allow-paired-ramda`

Estimated RNA cell count from the current validated public set captured in the repo:

- 3 RNA-side cell libraries

## Why demux is not required here

The practical unit exposed by the public metadata is already the per-cell RNA run:

- each candidate RNA-side cell is represented by its own SRR accession
- IMR90 uses one FASTQ object per SRR because the runs are single-end
- HAP1 uses one R1/R2 FASTQ pair per SRR because the runs are paired-end
- the current manifest abstraction maps naturally to one row per SRR-backed cell library

That is fundamentally different from a pooled full-length design where one sequencing object contains multiple cells that must be separated by barcode parsing before per-cell alignment and circRNA detection.

## Expected processing route

### IMR90

Processing route:

`single-end FASTQ -> BWA-MEM -> direct SAM -> CIRI3 -> matrix -> h5ad`

Manifest expectation:

- one row per IMR90 SRR
- no demultiplexing stage

### HAP1

Processing route:

`paired-end FASTQ -> STAR -> BWA rescue -> CIRI3 STAR tuple mode -> matrix -> h5ad`

Manifest expectation:

- one row per HAP1 SRR
- no demultiplexing stage

## Why Smart-seq3 differs

The validated `E-MTAB-8735` Smart-seq3 workflow in `circyto` is structurally different:

- sequencing starts from pooled FASTQs
- per-cell libraries are recovered by a demultiplexing step
- only after demux does `circyto` build a per-cell manifest and continue with alignment, circRNA detection, matrix export, and `h5ad`

So Smart-seq3 needs:

`pooled FASTQ -> demux -> per-cell manifest -> alignment -> CIRI3 -> matrix`

The current public RamDA/scRR datasets instead start effectively at:

`per-cell FASTQ(s) -> manifest -> alignment -> CIRI3 -> matrix`

## Implications for circyto workflow design

Current implication:

- the present `full-length-circrna` architecture is already aligned correctly for the public human RamDA/scRR datasets

That means:

- no new demultiplexing logic is needed for `GSE278958` IMR90
- no new demultiplexing logic is needed for `GSE278952` HAP1
- one-manifest-row-per-SRR remains the right abstraction

Future implication:

- if later RamDA-like datasets appear as pooled barcode-multiplexed libraries, they should be treated as a different structural class and should not be forced into the current direct-per-run assumption without metadata proof
