# Current Project Status

This page is the short operational snapshot for `circyto` as of the current manuscript-preparation phase.

## Validated

- `SMART-Seq3 all192`
  - public pooled full-length workflow validated
  - demultiplexing, circRNA matrix generation, and downstream integration layer established
- `IMR90 full23`
  - human scRR / RamDA-like single-end full 23-cell run validated on hg38
  - route: `BWA-MEM -> direct SAM -> CIRI3 -> matrix -> h5ad -> RNA profile -> RNA+circ MuData`
  - processed GEO CNV states and mappability-normalized signals imported at 50 kb
  - GSM-to-biological-cell remapping validated
  - tri-modal RNA+circ+CNV MuData validated:
    - `rna`: 23 x 63187
    - `circ`: 23 x 2443
    - `cnv`: 23 x 60607
    - trimodal overlap: 23
- `HAP1 batch10`
  - human scRR / RamDA-like paired-end 10-cell batch validated on hg38
  - route: `STAR -> BWA rescue -> CIRI3 -> matrix -> h5ad -> RNA profile -> MuData`
  - cleanup-workflow reduced the workdir from about 132G to about 9.2G while preserving final artifacts
- `HAP1 processed RT/state integration`
  - `circyto import-scrr-rt` and `circyto merge-scrr-rt` implemented with synthetic tests
  - real GSE278952 processed-file validation is pending local availability of the named HAP1 DNA RT/state files
- `SComatic interoperability`
  - RamDA/scRR adapter validated on HAP1 batch10
  - CB tags injected, NH/nM tags preserved when present
  - BaseCellCounter, Step1, Step2, and real output normalization validated as a technical smoke
  - terminology guardrail preserved:
    - `RNA-derived candidate variant signals`
    - not validated somatic mutations

## Deferred

- `HAP1 full`
  - full public HAP1 set is pending full FASTQ download and full workflow run
- `IMR90 SComatic whole-genome candidate calling`
  - exploratory / in progress
  - candidate outputs must remain RNA-derived candidate variant signals unless orthogonal DNA validation exists
- `Raw DNA FASTQ reprocessing`
  - not implemented inside circyto; current CNV and RT support imports processed GEO-style summaries

## Current manuscript-ready layers

- `Layer 1`
  - core circRNA detection workflows
  - public dataset reconstruction infrastructure
  - workflow summaries, QC, cleanup, integrity, provenance
- `Layer 2`
  - RNA profile import / simple-overlap counting
  - RNA QC refresh
  - RNA+circ joined summaries
  - MuData export and downstream scverse inspection
  - benchmark/status/report scaffolding
- `Layer 3`
  - processed scRR CNV import
  - processed scRR replication timing/state import
  - GSM-to-biological-cell mapping
  - tri-modal RNA+circ+CNV and RNA+circ+RT MuData merge
  - SComatic RNA-derived candidate variant signals interoperability as an optional sidecar

## Caution

Validated outputs and pending outputs should be described separately.

- validated:
  - IMR90 full23 RNA+circ+CNV tri-modal MuData
  - HAP1 batch10 RNA+circ and SComatic technical smoke
- pending:
  - HAP1 full workflow
  - real HAP1 GSE278952 RT/state processed-file import validation
  - IMR90 whole-genome SComatic candidate calling
  - raw DNA FASTQ reprocessing
