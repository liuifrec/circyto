# scRR DNA Roadmap

## Recommendation

The biologically justified near-term DNA branch for circyto should be:

```text
RNA + circ + CNV
```

not:

```text
RNA + circ + SComatic
```

SComatic interoperability remains useful as an exploratory RNA-derived candidate-signal branch, but it should not define the main scRR DNA architecture.

## Rationale

The scRR DNA modality is derived from scRepli-seq. The publication emphasizes genome-wide copy number, replication state, replication timing, S-phase progression, CNV detection, haplotype-specific analysis, and genome instability.

GEO exposes processed CNV tables by genomic bin size. These files map naturally into a `cnv` AnnData modality and can be paired with RNA by sample suffix.

By contrast, the SComatic path uses RNA alignments and produces candidate variant signals under assumptions that are not central to scRR-seq. The completed HAP1 SComatic validation is a technical integration success, but Step2 produced no candidate calls in the single-cell-type test design. That outcome is consistent with treating SComatic as optional exploratory interoperability rather than the core DNA branch.

## Low-Risk Path

1. Import processed GEO CNV state tables with `circyto import-scrr-cnv`.
2. Preserve original DNA IDs, inferred RNA IDs, canonical IDs, bin coordinates, bin size, and provenance.
3. Add `cnv.h5ad` as a first-class DNA output under `WORKDIR/dna/`.
4. Extend MuData export to include `cnv` when `dna/cnv.h5ad` exists.
5. Add lightweight per-cell CNV summaries only after the raw imported matrix is stable.

This path uses public processed outputs and avoids rerunning DNA workflows.

## High-Impact Path

After the low-risk import path is stable:

- correlate circRNA abundance with CNV burden and chromosome-arm copy-number state
- compare circRNA expression across replication-state or cell-cycle groups
- use paired RNA and DNA to analyze gene expression dynamics during S-phase progression
- integrate aphidicoline and genome-instability conditions where available
- preserve haplotype-specific hooks for future allele-aware analysis

These analyses align with the paper's biological framing more directly than RNA-derived SNV calling.

## SComatic Branch

Keep SComatic under a conservative namespace:

```text
candidate_snv
```

Required terminology:

- `RNA-derived candidate variant signals`
- not `validated DNA somatic mutations`

Use cases:

- technical interoperability with external SComatic outputs
- exploratory burden summaries when cell-type design is appropriate
- optional association with circRNA host genes

Do not use this branch to make primary scRR DNA claims unless orthogonal DNA validation exists.

## Manuscript Opportunities

The stronger circyto manuscript angle is multimodal infrastructure for full-length single-cell RNA/circRNA plus scRR DNA copy-number state:

- reusable `rna + circ + cnv` MuData schema
- public scRR GEO CNV importer
- direct DNA/RNA cell-pairing validation
- genome-bin CNV tracks linked to circRNA and RNA summaries
- explicit separation between copy-number DNA evidence and RNA-derived candidate variant signals

## Implementation Milestones

Current:

- real scRR RNA/circ path validated
- SComatic technical path validated, with no candidate calls in the single-cell-type HAP1 design
- `circyto import-scrr-cnv` parser skeleton added for processed CNV summaries

Next:

- validate importer on one small processed GEO CNV table
- extend `export-mudata` to include `cnv`
- add CNV-specific inspection summaries
- document exact genome build and accession provenance per imported dataset

Deferred:

- rerunning scRepliseq DNA processing
- CNV segmentation inside circyto
- SNV claims from RNA-only SComatic output
