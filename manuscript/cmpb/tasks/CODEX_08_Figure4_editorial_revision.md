# CODEX 08 — CMPB Figure 4 editorial revision (FIGURE ONLY)

**Execution status: not started.** Read this file in Codex, then perform only this scoped graphics/manuscript-evidence handoff if the owner instructs you to. This is NOT a request for new sequencing, scientific analysis, feature development, release, merger, or manuscript submission.

## 0. Exact source lineage and safety

The accepted A05 evidence commit is \`90381340cc0365bc66c815d06ff01b2978dc3d34\`; A03/A04 is \`02ee3d4acfed63c5b2f02c4baef6d22c822d4dac\`. The **A06 local** analysis commit is \`f858edff019b295c8e92e1c8645884ab9d2b721e\`. A06's owner-uploaded handoff ZIP has SHA-256 \`87d5c3c1879d7cc9a2bdea8bc2f801f3413b2cbe6a938efbcae897f6b4dfa542\`. Do not assume the A06 commit has been pushed: inspect local repository state and commit parents, never reset/cherry-pick/force-push prior analysis refs merely to access it. The task-planning GitHub branch was based on A05 for discoverability, **not** to replace A06 as the data-source commit.

Use an isolated worktree/branch from the A06 local commit (choose a new name, e.g. \`manuscript/cmpb-figure4-review\` only if unused). Preserve existing scientific outputs, accepted source branches, original Application Note figures, recheck trees, A100 failure records, original H5MU, prior A02–A06 files, and local execution roots. All production files remain unchanged. No \`pip install\`, reference/FASTQ/JAR download, CIRI3 invocation, raw analysis, new A03/A04/A05/A06 computation, release/version bump, or publication.

Read these existing files at their accepted commits:
- \`manuscript/cmpb/a03_a04/reference_comparison.tsv\` and \`operation_inventory.tsv\`, plus \`results_for_review.md\`.
- \`manuscript/cmpb/a06_scaling/scaling_runs.tsv\`, \`scaling_summary.tsv\`, \`scaling_inputs.tsv\`, \`scaling_figure_sources.tsv\`, \`results_for_review.md\`, and original \`figure4_synthetic_scaling.svg/png\`.
- \`manuscript/cmpb/figure_plan.md\` and \`manuscript/cmpb/a05_descriptive/cnv_eligibility_decision.md\`. The original Figure 4 layout is planning, not binding; write a small editorial decision explaining the new layout.

If any committed table, sample identity, exact measured value, or required run-level row is unavailable or conflicts with the reports, stop and report rather than reconstruct from remembered narrative.

## 1. Editorial decision — four panels in main Figure 4

Produce a **journal-legible, vector-first, 2x2 main Figure 4**. Each panel has a clear question and an explicit source table; no decorative title saying "Draft" inside the artwork. Use readable text at a typical two-column print width, sensible ticks and labels, consistent fonts, and a compact legend. Preserve monochrome distinguishability; no superfluous scientific assertions in the graphics.

- **4A (A03, genuine four-cell endpoint):** circyto versus independent bespoke reference, wall time (seconds). Show three genuine measured fresh-process repeats per arm, with median plus observed min/max; do not include warm-up as a fourth repeat. No significance tests or CI from three repeats. State these are **four real IMR90 cells**, native-to-RNA/multimodal objects/evidence/QC/offline HTML, startup included, warm-cache. Reference was faster: approximately 0.654 versus circyto 1.199 seconds median.
- **4B (A03):** the **same matched endpoint**, peak process RSS (MiB), three fresh repeats, median and observed min/max. Convert KiB to MiB with 1024 correctly. The independent reference used less process RSS (~156,000 KiB versus ~187,572 KiB). Use \`peak_rss_kib\` from A03's measurement rows, not A06 process-tree or cgroup peaks.
- **4C (A06 synthetic scaling):** registered S01/S02/S03 scenario labels plus actual measured native input bytes, biological cell count/feature/nnz descriptors and **wall time** (seconds). Plot all three measured repeats per scenario and median/range. Choose a clear axis (actual input size in MiB with explicit scenario labels, or categorical workload size), not an unjustified fitted curve. These scenarios simultaneously change cells, features, number of per-cell files and nnz: do not attribute scaling to nnz alone. S03: 10,000 synthetic cells, 50,000 features, 1,000,000 nnz, median 47.879 s.
- **4D (A06):** peak **process** RSS (MiB) for the same S01/S02/S03, all three measured repeats and medians/ranges, with consistent scenario positions. S03 median ~960.28 MiB; its max measured process RSS ~960.46 MiB. Do not equate with cgroup \`memory.peak\` or sampled tree RSS.

**Critical:** A03 is a *competent custom reference versus circyto* on four genuine matched IMR90 cells; A06 is *circyto alone* on artificial 100–10,000 cell records. They are two different experiments, and the figure/caption must not imply one shared scaling curve or a comparison to nf-core/circrna/CirComPara2. Both are post-detection conversion/report tasks, not complete alignment+caller times.

## 2. Supplementary controlled-variable figure

Create a clearly labeled supplementary pair of plots (or paired plot + table) for A06 \`S02\` (32-byte IDs, circ only), \`S04\` (128-byte IDs, otherwise identical native circ input), \`S05\` (same S02 circ inputs + synthetic RNA 1,000x2,000, 100,000 RNA nnz and MuData). Show raw repeats + medians for wall time and process RSS, and actual native/total input bytes and output bytes in a concise table or caption. S02-to-S04 median wall ratio is ~1.037; S02-to-S05 ~1.262. These are descriptive controlled contrasts on one machine, not generalized overhead effects.

Move redundant A06 plots of wall/RSS against multiple related size proxies to supplemental material only if useful. The existing six-panel A06 figure should be retained unchanged as a historical draft; generate **new file names** rather than overwriting it.

## 3. Accuracy, figure constraints and checks

- Rebuild **plots only** from the original TSV rows. Never rely solely on medians copied from reports; independently derive them from exactly three measured rows per arm/configuration, exclude controls/warm-ups, and verify against recorded summary values.
- Every visible data point/range/label must have a corresponding machine-readable figure-source row with original \`scenario\`, \`phase\`, input bytes, RSS and elapsed time. If jittering overlapping points on categorical axes, explain jitter as display-only; numerical values remain exact. Use proper log tick notation only if log axes are actually used. No normal-distribution CI, trend fit, benchmark significance claim, speed-superiority claim, synthetic biological inference or arbitrary drop of the reference's better result.
- Include vector SVG and high-resolution PNG for main and supplementary figures, exact \`figure4_main_sources.tsv\` and \`figure4_supp_sources.tsv\`, plotting script, source SHA-256s, an independent plot/value validation summary, complete captions with all benchmark limitations, and a **short** neutral draft Results paragraph referencing A03/A06 separately.
- Keep production, original figures and all accepted data unchanged. The descriptive A05 CNV gate is NOT VERIFIED: do not add any CNV panel or claim. The Figure 2C RNA/circ-only alternate remains a separate reviewed draft.
- Record revised editorial placement of A03 operation inventory (supplement/table, no quantified "human effort" claim). No full-paper rewrite.

## 4. Low-friction workflow and handoff

This is a **small figure formatting session**; aim for ~30–45 minutes, not a new scientific computation campaign. Reuse existing plotting tools/environments. Do not rerun the 20 A06 processes, regenerate 65,839 source files, run all 715 A04 scenarios, launch full production tests, download anything, or repeat the earlier full-source hash audit. Make lightweight source-table/plot/preservation checks appropriate to the changed files.

Create one small local review commit on the new figure branch after checking the diff; **do not push that analysis branch without a separate owner instruction**, merge, tag, publish, or submit. Put the new materials under \`manuscript/cmpb/figure4_review/\` or similarly isolated path. Export a **verified compact ZIP under \`~/Downloads/CIRCYTO_CMPB_FIGURE4_<timestamp>/\`**, containing figure SVG/PNG, scripts, complete source tables, captions, Results paragraph, hashes, and a README. No external JAR, H5MU, FASTQ, environments, bulk synthetic inputs/objects or private machine paths.

Return the figure preview paths, the ZIP absolute path/size/SHA-256, local commit, figure-check outcome, and *what scientific limitations the caption explicitly preserves*. Stop for owner review before any subsequent figure or manuscript tasks.
