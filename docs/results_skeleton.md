# Results Skeleton

This page is a cautious writing scaffold for the current `circyto` project.

## 1. circyto supports validated circRNA workflows across multiple full-length single-cell protocols

Points to cover:

- `SMART-Seq3` pooled workflow validated
- single-end human scRR / RamDA-like workflow validated
- paired-end human scRR / RamDA-like workflow validated
- workflow outputs are standardized into matrices, `h5ad`, and optional `h5mu`

Language guardrail:

- describe these as validated workflow and output layers
- do not overstate biological discovery from pilot-size runs

## 2. circyto reconstructs RNA+circ multimodal outputs from completed workflows

Points to cover:

- post-hoc RNA profiling
- RNA QC refresh after cleanup
- RNA+circ joined summaries
- MuData export for `scverse` interoperability

Language guardrail:

- focus on interoperability and reproducibility
- keep exploratory downstream analyses clearly labeled

## 3. Public scRR pilots establish executable human RamDA/scRR routes

Validated completed outputs:

- IMR90 2-cell pilot
- HAP1 3-cell pilot

Running outputs:

- IMR90 23-cell server run
- HAP1 batch10 server run

Language guardrail:

- separate completed validated outputs from currently running outputs
- avoid treating in-progress server runs as final results until post-run inspection is complete

## 4. circyto provides a cautious scaffold for future DNA/RNA/circ integration

Points to cover:

- schema contracts for DNA summaries
- SComatic synthetic integration validated
- future multimodal roadmap defined

Terminology guardrail:

- SComatic outputs should be called:
  - `RNA-derived candidate variant signals`
- they should not be described as validated somatic mutations without orthogonal DNA evidence

## 5. Real local SComatic execution remains deferred

Points to cover:

- synthetic integration path validated
- real WSL local setup blocked by native conda `Bus error`
- future real backend should move to:
  - HPC/server conda
  - or container / `mamba` / `micromamba`

Language guardrail:

- classify this as a deferred infrastructure step
- do not imply real SComatic calling has already been validated locally
