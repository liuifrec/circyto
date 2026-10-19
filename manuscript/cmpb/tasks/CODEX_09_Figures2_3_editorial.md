# CODEX 09 — CMPB Figures 2 and 3: science-led editorial design

**Staged task only.** Execute only when owner asks Codex to run this file. No analysis reruns, production edits, new detector workflow, A07, manuscript submission, release, push or merge.

## Provenance and worktree

Accepted A02 \`307758f57a39c009393d03fc0f69e4638b2498c6\`; A03/A04 \`02ee3d4acfed63c5b2f02c4baef6d22c822d4dac\`; A05 \`90381340cc0365bc66c815d06ff01b2978dc3d34\`. A06 **local** \`f858edff019b295c8e92e1c8645884ab9d2b721e\`; accepted Figure 4 **local** \`81ffd7509fcaa4cbb44b7f9e82064edfd40e9837\` (branch \`manuscript/cmpb-figure4-review\`). **Verify** these local refs and ancestry, do not assume A06 or Figure 4 has been pushed. This planning branch starts from A05 only to expose this task file. Create a separate \`manuscript/cmpb-figures2-3-review\` worktree from the **local Figure 4 review commit**, or stop if absent/ambiguous. Existing originals, prior outputs, BioRender, manuscripts, worktrees and failure records remain immutable.

Read \`manuscript/cmpb/figure_plan.md\`, A02 \`cohort_acceptance.md\` and \`pilot_cell_trace.tsv\`, A03/A04 \`results_for_review.md\` and \`robustness_outcomes.tsv\`, A05 \`a05_results_for_review.md\` and \`cnv_eligibility_decision.md\`, and Figure 4 \`captions.md\` / SVG to align style. Resolve accepted datasets using their hashes and manifest paths. Do not substitute a review extract for a complete source table.

Optional full A04 completed-cohort candidate table: \`a04_combined/completed_cohort_prevalence.tsv.gz\` under the accepted evaluation root, expected 71,757,173 compressed bytes, 5,159,619 data rows, SHA256 \`8163178301f702fb4217fc0e48d4390ed82aa39244a8d1efd2cc681ad9152b8c\`. Find its root from the A04 manifest, not by guessing paths. Hash and stream if available; if missing, choose a truthful lower-resolution alternative panel. **Never rerun 715 perturbations.**

## Figure design objective

Produce visually sophisticated, *journal-readable* two-column-width **Figures 2 and 3**: clean hierarchy, real identity connections, exact-data mini-matrices, state logic, source-anchored quantitative plots, consistent typography and restrained colorblind-safe palette. SVG, 600-dpi PNG, 180-dpi preview and grayscale proofs; readable around 180 mm print width. Keep information dense but not crowded; no needless icons, artistic data manipulation, generic stock figures or decorative science. Existing BioRender Figure 1A remains unchanged.

### Figure 2: genuine identity-preserving integration and biological-data context

**2A** – True four-cell A02 pipeline/evidence (A100–A103; accession/run labels as space permits). A100 142 candidates/474 support; A101 95/317; A102 112/419; A103 123/489. All four together have **435 coordinate/strand-compatible unique candidates, 472 positive cell-candidate entries, 1,699 total CIRI3 junction-read support**; distinguish all four quantities. 4 x 435 circ output. A100's original run and repaired status replay had mixed lineage with normal repaired A101–A103 execution. Never label all four as one raw-read workflow.

**2B** – **Actual exact-ID matching** from historical RNA input order A103,A102,A101,A100 to circ metadata order A100,A101,A102,A103. Draw restrained crossing connectors/ribbons, with source cells and real canonical ID labels, not invented merges. RNA historical slice 4 x 63,187 retains original RNA values and annotations; reference digest unknown, RNA simple-overlap not independently quantified. The distinct new same-run three-cell RNA remains separate. State matching by ID, not positional order.

**2C** – A05 *23-cell historical* IMR90 RNA/circ-only descriptive figure, redesigned as source-anchored aligned small multiples or a compact heat-strip. Source four metrics per original stored cell: circ-positive candidates, CIRI3 support, RNA positive features and historical simple-overlap signal; preserve original 23-cell order, all cells, actual physical-zero versus modality member semantics. This 23-cell historical scope is **different** from the A02 four-cell new detection experiment; all original per-cell completion statuses unknown. **CNV is NOT VERIFIED**; do not render a CNV heatmap or assume zeros are uncallable.

**2D** – Compact exact technical reconciliation (native evidence→circ AnnData→matched RNA MuData→QC/HTML), with genuine observed counts and an appropriately constrained source/value matrix. Distinguish fresh workflow denominator 3 from mixed-lineage native import denominator 4, and "technical exact equality" from biological truth. Avoid transforming thousands of assertions into biological replicates or claiming 100% circRNA accuracy. Full 435 candidates must not be silently sampled. If compact exemplars are used, predefine a non-outcome-based selection rule and label the subset.

### Figure 3: tested observation-state contract and perturbations

**3A** – Clear state/denominator logic: observed positive support, genuine completed-zero, failed/unknown/unprocessed and no-completed denominator→undefined prevalence. Genuine A02 only observed positive support. Completed-zero/failure/unknown controls are **simulated/synthetic**, not observed biological states in A02. No dropout correction claim.

**3B** – Full registered robustness landscape for R01–R11, including valid and invalid controls. A04: **715 configurations per arm** (705 valid, ten invalid). Both circyto and independent competent reference accepted/passed all 705 valid and rejected all ten invalid. Deliberately unsafe positional/zero-fill control passed 172 valid and failed 533 valid value/state checks. Its invalid rejection is **inherited from the competent validator**, not an independent unsafe validation success. Use source-linked scenario x arm tiles or grouped outcomes, not only the favorable scenarios. The unsafe method represents no external package or typical analyst.

**3C** – Graphical example of same-*completed*-cohort prevalence and undefined denominator. R10 has 2,443 undefined candidate results correctly retained in both competent implementations; the unsafe zero-fill creates misleading defined numbers. Same-completed-cohort error is distinct from selection bias relative to a full cohort. No false missingness-bias-correction claim.

**3D** – Error magnitude/types from the *complete* accepted A04 candidate table if present. Stream through it, reconcile arm/config totals and preserve exact counts. Prior A04 reports **781,440 nonzero candidate/configuration errors for unsafe** and zero in both competent arms, plus R10 undefined-to-defined outcomes. Rows are correlated candidate/configuration observations, not independent cells or error probabilities. If the full gzip is inaccessible, use \`robustness_outcomes.tsv\` and state the narrower endpoint. No cherry-picked striking genes/candidates; use registered scenarios and fixed selection rules.

Figure 3 should show *why a contract matters*, while transparently demonstrating that a competent ID-aware custom reference performs just as correctly as circyto in these tests.

## Accuracy, deliverables and handoff

- Plot from **accepted results only**. Every displayed numerical mark, source ID, state label, denominator and annotation must link to a machine-readable figure-source row and a recorded source hash. Build scripts may compute plotting aggregates; they must not mutate or rerun accepted scientific analyses.
- Output \`figure2_main.svg/png\`, \`figure3_main.svg/png\`, grayscale and compact previews, \`figure2_sources.tsv\`, \`figure3_sources.tsv\` (plus source-linked A04 aggregates if needed), source hashes, a high-quality caption per figure, a short draft Results paragraph per figure, plotting code and an independent numerical/SVG check not importing plot-generator code.
- Maintain existing Figure 1A/1B/1C, Figure 4, original A02–A06 outputs and all manuscripts unchanged. Do not claim a new analytical model, biological circularization validation, radiation effect, clone/CNV association, novel missingness estimator, speed superiority, or human-usability study.
- Figure-only target 60–90 minutes, up to 2 CPUs/4 GiB process memory/2 GiB **additional** storage. Treat limits as ceilings; use the already installed environment. No reinstallation, long resource-control tests, raw FASTQ, CIRI3, source index builds, or rescanning 65,839 synthetic files. Stop and report when a required input cannot be recovered within those limits.
- Review all output figures at final physical size and grayscale; repair overcrowding without omitting unfavorable or null outcomes. Keep accepted old A05 alternative/other figures unchanged, use **new** files only under \`manuscript/cmpb/figures2_3_review/\`.
- Commit **new figure-review materials only** to a local isolated branch, with clean diff; no push/merge/tag/publication without fresh approval. Export a compact verified ZIP under \`~/Downloads/CIRCYTO_CMPB_FIGURES2_3_<timestamp>/\`, with source tables, code, figures, captions, complete small review data and hashes. The full 71 MB A04 gzip may remain outside the ZIP if the SHA and complete derived aggregates are included. Verify all archive members, report ZIP absolute path/size/SHA256, local commit, figure audit, and any gated panel. Stop for owner review.

**After acceptance:** prepare an integrated CMPB manuscript draft; do not perform that step in this task.
