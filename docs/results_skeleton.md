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

## 3. Public scRR runs establish executable human RamDA/scRR routes

Validated completed outputs:

- IMR90 full23 RNA+circ workflow
- IMR90 full23 processed CNV import
- IMR90 full23 tri-modal RNA+circ+CNV MuData
- HAP1 batch10 paired-end RNA+circ workflow
- HAP1 batch10 SComatic technical smoke
- HAP1 processed RT/state importer implemented with synthetic tests

Pending outputs:

- HAP1 full public set
- real GSE278952 HAP1 processed RT/state file import validation

Language guardrail:

- separate completed validated outputs from pending full-data outputs
- avoid treating HAP1 full or exploratory SComatic candidate calling as final biological results

## 4. circyto integrates processed scRR DNA state as dataset-specific modalities

Points to cover:

- processed GEO CNV state import for IMR90-style CNV summaries
- mappability-normalized CNV signal layer where supplied
- processed replication timing/state import for HAP1-style RT summaries
- GSM-to-biological-cell mapping
- tri-modal RNA+circ+CNV or RNA+circ+RT MuData
- SComatic interoperability as an optional sidecar

Terminology guardrail:

- SComatic outputs should be called:
  - `RNA-derived candidate variant signals`
- they should not be described as validated somatic mutations without orthogonal DNA evidence

## 5. SComatic interoperability is technical and exploratory

Points to cover:

- RamDA/scRR adapter validated
- HAP1 batch10 BaseCellCounter, Step1, and Step2 executed
- Step2 produced no candidate calls in the single-cell-type design
- local WSL setup remains an environment caveat; server/container execution is preferred for real runs

Language guardrail:

- classify this as exploratory RNA-derived candidate variant signals interoperability
- do not imply validated somatic mutation calling
