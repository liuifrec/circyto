# CIRCYTO CMPB handoff and figure review - 2026-10-19

## Target

Primary: Computer Methods and Programs in Biomedicine, subscription/non-open-access route approved in principle by RERF Chief Scientist Dr. DeWees. TCBB is the alternate. This planning record does not replace the corrected institutional approval form or journal clearance.

## Accepted computational evidence and source references

A02 genuine four-cell native/RNA acceptance: `analysis/cmpb-a02-complete` SHA `307758f57a39c009393d03fc0f69e4638b2498c6`. Per-cell circ candidates/support: A100 142/474; A101 95/317; A102 112/419; A103 123/489. Four-cell union 435 candidates, 472 nonzero cell-candidate entries and 1699 CIRI3 junction-read support units. Mixed scientific execution lineage; A100 was replayed for postprocessing only. Historic same-cell RNA has four rows and 63187 features in reverse cell order, and original reference identity remains unverified. These are technical integrations, NOT independent validation of circularization.

A03/A04 completed: `analysis/cmpb-a03-a04` SHA `02ee3d4acfed63c5b2f02c4baef6d22c822d4dac`. User review ZIP received and its 40 entries pass ZIP CRC. Complete 5,159,619-row candidate prevalence table has SHA256 `8163178301f702fb4217fc0e48d4390ed82aa39244a8d1efd2cc681ad9152b8c`, gzip intact. Both circyto and competent bespoke ID-aware reference passed 705/705 valid simulation configurations and rejected all 10 invalid ones; unsafe illustrative positional/status-filling arm failed at least one value or state check in 533/705 valid configurations. The unsafe arm shares the reference validator and is NOT an independent comparator to nf-core/circrna or CirComPara2.

A03 four-cell exact-task fresh-process medians: circyto 1.199263 s and 187572 KiB RSS; bespoke reference 0.653512 s and 156000 KiB RSS. The reference was faster and leaner at this size. The numbers are startup-inclusive small-workload comparisons, not scaling, developer effort or human usability. Circyto provides a reusable, tested workflow; no speed superiority claim.

## Figure disposition

- Figure 1A: preserve approved BioRender diagram exactly. Figure 1B/C: preserve audited Smart-seq3 panels; 192 stored identities do not mean 192 verified completed circ assays; distinguish padded RNA-only DIYHEK_192.
- Figure 2A/B/D: A02 source evidence established; begin paper-quality assembly after confirming source formats. Figure 2C: WAIT for A05 CNV semantics gate. If unsupported, use an explicitly identified RNA/circ-only panel instead of a speculative CNV heatmap.
- Figure 3A/B: valid/invalid configuration outcome comparison and correct-vs-unsafe semantics, explicitly simulation. Both competent arms achieved perfect results on the tested grid; do not imply superiority over the reference.
- Figure 3C: same-completed-cohort prevalence correctness and undefined-denominator preservation; full candidate table available. Simulation seed repetitions are NOT independent biological samples and large candidate-case counts are NOT p-values or external replications.
- Figure 4A: task/operation inventory. No usability/user-effort measurement, no named-tool execution; rely on documentation-based Table 1 for nf-core/CirComPara2 comparison.
- Figure 4B: A03 equivalent-endpoint timing/RSS, showing the reference faster.
- Figure 4C/D: A06 synthetic scaling, still pending. All synthetic panels labelled as such; no claims of detector/runtime end-to-end scalability until measured.

## Immediate journal-quality scientific writing boundaries

Results should be evidence-first: A02 = real matched-cell acceptance; A03 = correctness equivalent to a competent bespoke reference, with measured overhead; A04 = correct state/identity handling under prespecified perturbations relative to deliberately unsafe heuristic; A05 = real stored-signal characterization; A06 = synthetic scaling. Comparisons do not establish new circRNA detection sensitivity/specificity, biological prevalence or caller accuracy. Historic run provenance remains incomplete and must be disclosed separately from new technical rerun.

Main manuscript should NOT call 523 passing software tests, 2145 case/arm runs, or 5,159,619 candidate rows independent biological replications. No hypothesis tests against invented or unsafe baselines. No claim of workflow package publication before actual PyPI/Bioconda release. Figure labels need data version, preprocessing units and separate captions for real versus simulated evidence.

## What to ask Codex for next

Start with CODEX_06_A05_descriptive_and_CNV_gate.md. Do not launch A06 before A05 review. Codex must provide a direct ZIP under ~/Downloads with checksum, rather than only /tmp pointers, and small source tables + figures. Bring the A05 ZIP back to ChatGPT; together we will decide the final Fig 2C and scientific Results prose before scheduling A06.
