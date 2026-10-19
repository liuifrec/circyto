# CIRCYTO / CMPB - Codex task A05 (execute now)

## Objective and decision boundary

Execute ONLY the registered A05 descriptive analysis for the primary journal Computer Methods and Programs in Biomedicine (CMPB), subscription/non-open-access route. A02 genuine matched-cell acceptance and A03/A04 reference-conversion/robustness studies are already completed. Do not rerun them. The purpose is to produce source-linked, reproducible within-dataset descriptive evidence for an IMR90 RNA/circRNA/CNV use-case figure, including a mandatory CNV missingness/callability gate.

Manuscript interpretation and publication-ready figure selection will be reviewed in ChatGPT after your factual handoff. This task does not approve any manuscript rewrite, package release, additional workflow feature, detector comparison, or other experiment.

## Frozen scientific sources and branches

- Repository: liuifrec/circyto.
- Accepted A03/A04 head: 02ee3d4acfed63c5b2f02c4baef6d22c822d4dac (analysis/cmpb-a03-a04). Verify it and base a NEW worktree and branch `analysis/cmpb-a05-descriptive` on that exact commit. Never reset, rewrite or merge main, manuscript, PR #19, or prior analysis branches.
- Repaired evaluated production implementation: e6843a38469d49f96d44758a4e37893ebbe2734b. No production change is authorized.
- Original processed object archive: c99cddaae2314d3af7d5232150c5b86bb6b52e25.
- Smart-seq3 object: `load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu`; expected SHA256 `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339`.
- IMR90 object: `load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu`; expected SHA256 `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.
- Source registration: `manuscript/cmpb/study_design.md`, `manuscript/cmpb/benchmark_plan.yaml`, `manuscript/cmpb/figure_plan.md`, and A01 data/claim inventories. Re-read the actual current files; do not assume these summary instructions supersede a verified conflict in the original registration. Record any amendment before execution.
- Existing A02/A03/A04 evidence and derived objects remain read-only, including the source four-cell records, the accepted historical failures, and the complete prevalence dataset. Do not clone large source trees or upload data.

Verify the exact source H5MU hashes and modality axes before any analysis. If a required original object is inaccessible, stop and specify the exact missing source. Do not manufacture a replacement from A02 or a different processed cohort.

## Design: retain the source meaning

1. Analyze Smart-seq3 and IMR90 separately; they differ in assay/material, so no cross-dataset superiority/biological-yield inference.
2. Confirm shape/orientation/nonzero/support totals against the audited original records, documenting any mismatch rather than replacing old results:
   - Smart-seq3: 192 stored identities, 2503 circ candidate features, 2659 circ nonzeros, 5711 circ support, 2379 recorded host-annotated candidates.
   - IMR90: 23 matched identities, 2443 circ candidate features, 2914 circ nonzeros, 7422 circ support, 2429 recorded host-annotated candidates. Imported CNV representation: 23 x 60607 bins.
3. Smart-seq3 `DIYHEK_192` is RNA-only with a padded circ zero row; explicitly separate physical object rows, modality membership and original circ assay completion. Unknown historical completion MUST remain unknown. Do not treat 192 as 192 verified completed circ assays.
4. The historical RNA is a lightweight simple-overlap profile; its original reference digest has not been recovered. Do not describe it as independently quantified transcriptomics or a newly verified reference.
5. All per-cell and per-candidate metrics are descriptive of stored measurements, not true biological prevalence, donor-independent replication, or circRNA detection sensitivity.

## A05-A: scientific summary tables

Compute from each original sparse matrix without modifying the H5MU:

- Per cell: exact biological cell ID, modality membership, stored circ row presence, historical processing status if recorded (else unknown), positive circ candidate count, sum of support, total RNA overlap-profile signal, RNA positive-feature count; optional zero support-row flag distinct from completed-zero.
- Per circ candidate: original ID, coordinates/strand, host-gene annotation fields and provenance, number of stored cells with positive support, total junction-read support, maximum per-cell support. Distinguish stored-cell fraction from verified completed-cohort prevalence; do not calculate the latter when statuses are unknown.
- Characterize within-dataset distributions with clear denominators and medians/quantiles. Carry all zero and missing members; do not drop low/zero support by aesthetic preference.
- Prespecified threshold sensitivity, separately: (A) per-cell support >= 1, 2, 3, 5 on individual entries and (B) pooled candidate support >= 1, 2, 3, 5 across the original matrix. Start EVERY threshold from the original X, not from sequential filtering. Do not recommend a threshold based on these descriptive plots.
- A cautiously labelled RNA-versus-circ burden view, within each dataset only and only for actually matched identities. Use descriptive scatter/aggregates without causal interpretation or cross-assay claims. Any correlation coefficient must state N, choice of metric and limitations; no selected positive example.

Deliver source data, concise figure drafts, and a factual Methods/Results data note. Keep every value traceable to the two original object hashes. All source data for plotted observations should be retained in UTF-8 TSV/CSV plus provenance metadata.

## A05-B: IMR90 CNV eligibility gate (DO NOT bypass)

Before aggregating or plotting CNV, inspect exact original /mod/cnv fields, source metadata, associated GEO state and mappabilitynorm origins, bin coordinates, numerical encodings, and any documented missing/uncallable/callability mask. Check explicit NaN/missing values and layers, but DO NOT infer that numerical zero means missing, callable or diploid. Successful H5MU import does NOT certify original bin callability.

Decision:
- PASS only if a source-supported missing/uncallable definition and compatible numeric semantics are recovered and recorded. If PASS, create a display-only fixed 1-Mb median of valid underlying 50-kb bins; require >= 50% valid covered bases per window and record contributing bin IDs, number of bins, valid-base fraction, source range, aggregated value, display order. Keep all 23 cells in source/stably defined order. Mark invalid windows missing. Never overwrite native 50-kb states or claim the plot diagnoses CNVs.
- NOT VERIFIED if missingness/callability cannot be established. In this case do NOT draw a misleading CNV heatmap or invent a mask. Provide `cnv_eligibility_decision.md` with exact missing evidence and a proposal for a descriptive RNA/circ-only substitute for Figure 2C. The remaining A05 metrics MUST still be completed.
- FAIL if source conflicts are found; preserve them, stop ONLY the CNV panel, and report the contradiction.

No radiation/clonality/circRNA-causation assertion, chromosome-based selection by attractive appearance, unregistered association test, or post-hoc change to the 1-Mb rule.

## Inputs, reproducibility and execution controls

Use a NEW results root outside the object archive and all historical workdirs. Reuse the existing accepted Python/AnnData/MuData environment; avoid dependency installs, large downloads, raw read acquisition, BWA/STAR/CIRI3 execution, genome indexing and all A03/A04 grid runs.

Use the registered A05 resource proposal as the initial maximum: 2 CPU cores, 4 GiB controlled memory, 2 GiB ADDITIONAL persistent outputs, 1 hour active job window. Report exact actual resource accounting and enforcement; if the cap cannot be implemented as previously verified, stop for owner review. Do not silently extend budgets. Small input inventories, scripts, TSVs and vector figures are sufficient; no dense conversion of full genomic matrices.

At minimum, verify source matrix shape, exact sums/nnz by independent sparse operations, cell ID uniqueness and matching, missing handling, deterministic figure-source ordering and all plotted values. Keep original H5MU unmodified and record before/after SHA256. Annotate when a technical field is not available rather than inventing it. Existing datasets and generated summaries must be distinguishable from simulated A04 records and the new four-cell A02 rerun.

## Outputs and handoff

In `manuscript/cmpb/a05_descriptive/` commit only small factual records and source-linked plot drafts (SVG and PNG as useful), NOT raw H5MU or large per-bin arrays. Suggested files:

- `a05_results_for_review.md` - factual findings, actual counts and limitations; no polished claims.
- `descriptive_cell_metrics.tsv`
- `descriptive_candidate_metrics.tsv`
- `threshold_sensitivity.tsv`
- `modality_membership_audit.tsv`
- `cnv_eligibility_decision.md` and `cnv_source_inventory.tsv`
- `cnv_display_bins.tsv` ONLY if eligibility PASS (else not created, not a fake empty PASS table).
- `figure2c_draft.svg` and `.png` or a clearly named RNA/circ-only alternative.
- `inputs_outputs_sha256.json`, `execution_record.json`, `figure_data_dictionary.md`.

Any plot must have source table, source object hash, code/helper commit SHA, precise axis label/units and denominator; no opaque pre-rendered images without source. Do not edit or replace existing Figure 1, circulated manuscript, or original Biorender materials. Commit/push ONLY the new A05 analysis branch after validation; report exact head and whether push succeeded. Do not open a PR, merge, tag, release or submit.

CRITICAL USER HANDOFF: Create a verified ZIP in a NEW, non-overwriting directory under `~/Downloads/CIRCYTO_CMPB_A05_20261019/` with the small committed source tables, figure drafts, report, scripts needed to recreate plots, SHA256SUMS and README. Do NOT include raw reads, 23-cell original H5MU, CNV full raw arrays, environments, executables, JARs, genome indices or secrets. Report the exact ZIP path, size and SHA256 so the user can upload it to ChatGPT. If `~/Downloads` is unavailable, use a clearly reported user-owned alternative, not a transient /tmp path.

STOP after A05. Do not start A06 until its separately prepared scope is approved following review of this evidence. If CNV fails eligibility, that is a valid A05 outcome, not an invitation to change the source rules.
