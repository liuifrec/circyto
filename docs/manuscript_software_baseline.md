# Manuscript software baseline

## Identity

This document is the authoritative circyto software snapshot for the
Bioinformatics Application Note.

| Field | Frozen value |
| --- | --- |
| Project | circyto |
| Software version | 0.10.0 |
| Candidate manuscript baseline commit | `3eb4acdf1c6f55c6edd88e8aa4e21f8f15376d06` |
| Validation date | 2026-08-24 |
| Full source and clean-wheel validation | Python 3.10.20 |
| Targeted MuData 0.4 compatibility check | Python 3.12.14; not a full release-stack validation |

The primary final-validation environment used NumPy 2.2.6, SciPy 1.15.3,
Pandas 2.3.3, AnnData 0.11.4, MuData 0.3.10, h5py 3.16.0, Typer 0.27.0,
and pytest 9.1.1. The clean wheel resolved the same numerical/scverse stack
and Typer 0.27.1. The targeted Python 3.12 check used MuData 0.4.1 and AnnData
0.13.2 only for the synchronization compatibility paths described below.

The candidate tree is the exact tip of the PR #13 to PR #15 stack. It contains,
in order, the issue #12 portability fixes, biogenesis schema/export foundation,
generic Nanopore interoperability, CIRI-long adapter and validations,
release-readiness cleanup, and explicit MuData synchronization compatibility.

## Core manuscript-safe capabilities

1. **Short-read full-length single-cell circRNA workflow.** circyto provides
   detector orchestration, cell-level collection, circRNA annotation,
   circRNA-by-cell matrix construction, `h5ad` export, and RNA+circ MuData
   integration. Advanced detector interfaces and high-level workflows retain
   their documented experimental status.

2. **Full-length single-cell protocol support.** The implemented paths cover
   pooled Smart-seq3 and manifest-driven RamDA, Shin-RamDA, and scRR inputs,
   including their documented single- and paired-end routes. This is a
   software capability statement, not biological validation beyond the
   evidence already committed for each route.

3. **Multimodal interoperability.** The software constructs RNA+circ objects
   and can add exploratory RNA-derived `candidate_snv`, processed replication
   timing/state (`rt`), and processed copy-number (`cnv`) modalities. MuData
   construction, update, read, and write operations explicitly preserve
   circyto's validated `pull_on_update=True` synchronization policy.

4. **Generic Nanopore interoperability.** The experimental Nanopore path
   provides provenance-aware ingestion, archive identity/integrity checks,
   minimap2/samtools alignment, query-level QC, and conservative alignment
   diagnostics. It does not perform or validate circRNA detection.

5. **CIRI-long interoperability.** The optional adapter chemistry-gates inputs
   to RCRT/circRNA-enriched libraries and supports separate call, collapse, and
   import stages. Imports preserve BSJs, full-length isoforms, expression,
   isoform usage, and read assignments.

6. **Biogenesis export.** Versioned dependency-light schemas and validation
   cover candidate, cell-context, and cell-by-candidate observation tables.
   Labels use positive/unlabelled semantics; no predictive model or true-negative
   interpretation is implemented or claimed.

## Validation evidence

The [official CIRI-long demo validation](validation/ciri_long_official_demo.md)
recorded 149 BSJs, 149 isoforms, 149 expression rows, 149 isoform-usage rows,
and 1,133 read assignments. All five normalized artifacts matched their
official inputs with zero unmatched identifiers and zero value discrepancies.
This is detector interoperability evidence from a bulk official demonstration,
not single-cell biological performance evidence.

The [SRR4048177 Nanopore interoperability validation](validation/srr4048177_nanopore_interoperability.md)
processed 52,696 reads. It observed 39,537 mapped primary queries (75.03%), of
which 33,270 were spliced primaries (84.15% of mapped primaries). These results
validate generic long-read ingestion, alignment, QC, and provenance only.

The final package gate produced **345 passed, 8 skipped, 5 warnings**. The
remaining warnings are four intentional overlapping-feature `var_names`
warnings and one AnnData string-index conversion warning. The sdist and wheel
built successfully; diff and compile checks passed. A clean Python 3.10.20
environment installed only the wheel plus declared dependencies, resolved
imports from `site-packages`, passed `pip check`, exercised the requested CLI
and diagnostic surfaces, found the packaged CIRI3, CIRI-full, CIRI2, and
Nanopore resources, and completed a warning-free RNA+circ H5MU round trip.
Default tests required neither network access nor detector executables. See the
[release-readiness audit](release_readiness_audit.md) for the preceding package
audit.

The [MuData compatibility audit](mudata_compatibility.md) reduced 126
behavior-change warnings to zero by making the already-validated synchronization
policy explicit. Regression tests verified logical equivalence for cell and
feature identity/order, metadata and missingness, modality inventory, maps,
matrix shapes/values, and representative H5MU round trips. Broader AnnData
0.13 dtype migration is not claimed.

## Explicit non-claims

- circyto is not a new circRNA calling algorithm.
- Generic ONT cDNA is not validated for circRNA detection.
- The official bulk CIRI-long demo does not prove single-cell CIRI-long
  biological performance or biological accuracy.
- Independent RCRT validation on CRR194209 is deferred.
- Biogenesis modelling is a future extension; this baseline includes schemas
  and export only.
- Unobserved or unlabelled circRNA candidates are not assumed to be true
  negatives.
- Advanced detector interfaces and high-level workflows remain experimental
  where documented; this snapshot is not a frozen v1 API promise.

## Frozen or deferred engineering

The following work does not block manuscript writing:

- CRR194209 independent RCRT validation;
- broad MuData 0.4 / AnnData 0.13 dtype migration;
- biogenesis predictive modelling;
- additional detector integrations;
- CLI redesign; and
- repository history slimming.

## Freeze decision

**PASS.** PR #13 independently passed its tests, build, compile, diff, and CLI
checks. The exact simulated PR #13 to PR #15 final tree passed the full source,
package, clean-wheel, semantic, and repository-hygiene gates. This software
tree can serve as the stable manuscript baseline without adding biological
scope. PR #13 and PR #15 remain unmerged at the time of this record.
