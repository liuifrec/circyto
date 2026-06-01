# Current Project Status

This page is the short operational snapshot for `circyto` as of the current manuscript-preparation phase.

## Validated

- `SMART-Seq3 all192`
  - public pooled full-length workflow validated
  - demultiplexing, circRNA matrix generation, and downstream integration layer established
- `IMR90 2-cell`
  - human scRR / RamDA-like single-end pilot validated on hg38
  - route: `BWA-MEM -> direct SAM -> CIRI3 -> matrix -> h5ad -> RNA profile -> MuData`
- `HAP1 3-cell`
  - human scRR / RamDA-like paired-end pilot validated on hg38
  - route: `STAR -> BWA rescue -> CIRI3 -> matrix -> h5ad -> RNA profile -> MuData`
- `SComatic synthetic integration`
  - local chr21 synthetic candidate import and DNA/RNA/circ summary validated
  - terminology guardrail preserved:
    - `RNA-derived candidate variant signals`
    - not validated somatic mutations

## Running

- `IMR90 23-cell`
  - full public scRR reconstruction running on the server
  - intended output stack:
    - `full-length-circrna`
    - `simple-overlap` RNA profiling
    - RNA QC refresh
    - RNA+circ summary
    - MuData export
- `HAP1 batch10`
  - paired-end public scRR reconstruction running on the server
  - intended output stack is the same, with paired RamDA route enabled

## Deferred

- `Real SComatic local WSL smoke mode`
  - deferred
  - repeated native `Bus error (core dumped)` during minimal conda package installation
- `Real SComatic backend`
  - future target on:
    - HPC / server conda
    - or container / `mamba` / `micromamba`

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

## Caution

Validated outputs and running outputs should be described separately.

- validated:
  - completed pilots and tested local proof-of-concept paths
- running:
  - currently executing server jobs that still need post-run inspection before any result claims
