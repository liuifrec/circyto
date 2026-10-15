# CMPB figure and table registration

2026-10-15. Four figures and two main tables are planning choices, not verified
CMPB limits. This file implements the supplied expansion blueprint; the M0
Application Note figure/table plans and circulated files remain unchanged.
All new assets are **PENDING**, with no plots or previews generated this session.

| Item / question | Inputs and analysis | Required source tables or assets | State and interpretation limit |
| --- | --- | --- | --- |
| Figure 1A: what does the framework connect? | Existing approved BioRender panel | Original supplied BioRender file, kept exactly | PRESERVE. The historical script-generated architecture SVG is not a replacement. No new schematic generated. |
| Figure 1B/C: what is the established Smart-seq3 output? | A01; original checksum-matched object; existing RNA-derived coordinates | M0 `manuscript/results_v2/smartseq3_umap_cells.tsv`, `smartseq3_selected_candidate.tsv`, `figure1b_smartseq3_circrna_burden.svg`, `figure1c_smartseq3_man1a2_detection.svg` and their recorded hashes | PRESERVE audited assets. Circ support did not define embedding. Clarify 192 stored identities and RNA-only padded row; do not label 192 completed circ assays. |
| Figure 2A: which real cells and inputs were used? | A02 fixed GSM8558852–GSM8558855; original recovery A or new rerun B | `pilot_cell_trace.tsv`: canonical/GSM/SRX/SRR/RNA/DNA IDs, input hashes, layout, producer/ref hashes, status/exit code, native hashes, inclusion rule | PENDING. Public mapping and same-cell RNA available; original native/status absent. Distinguish recovery from NEW validation. |
| Figure 2B: are native values preserved? | A02 independent native CSV oracle, imported CSR and matched RNA | `native_value_checks.tsv`: every cell/candidate support, coordinate, strand, gene/type annotation, axis, RNA/layer and round-trip check; totals and discrepancies derived from complete table | PENDING. Expected numerical equivalence is valid; no circRNA signal threshold for successful scientific acceptance. If all genuine outputs are zero, retain that outcome and reconsider panel scope. |
| Figure 2C: what can a matched multimodal view display? | A05 original IMR90 23-cell RNA/circ/CNV object, separately from the four-cell pilot | `descriptive_cell_metrics.tsv`; `cnv_display_bins.tsv` with source-bin IDs, source coordinates, valid fraction and aggregated value; missing-mask provenance | PENDING. Keep original 50-kb bins; display-only fixed 1-Mb median with >=50% valid covered bases. Uncallable/missing remain missing. No CNV causation, treatment effect or clone claim. Missing-mask uncertainty blocks this panel until reviewed. |
| Figure 2D: do reports show those values? | A02 authentic copied/new workflow directory versus separate native IMPORT bundle | `report_value_checks.tsv` linking native/oracle to workflow/import JSON, TSV and displayed HTML; separately labeled offline screenshots | PENDING. Two report types have different provenance and denominator contracts. A new workflow report does not recover original-workflow acceptance. |
| Figure 3A: do valid permutations preserve meaning? | A04 R01/R02/R04/R09 valid inputs; circyto and ID-aware reference | `robustness_outcomes.tsv` with analysis/config/seed/fraction/realized fraction, arm, expected validity, assignment and value errors | PENDING. Unsafe positional/zero-fill arm labeled illustration, never a proxy for nf-core or CirComPara2. |
| Figure 3B: what happens on invalid or incomplete observations? | A04 duplicate/mixed-declaration/missing-source errors and valid failed/unknown/partial-positive cases | Same complete outcome table; error type, exit code, retained values/statuses, false rejection on valid cases | PENDING. Show all controls and failures. Do not select only scenarios favorable to circyto. |
| Figure 3C: which denominator is being estimated? | A04 random and RNA-profile-depth-associated missingness; twenty fixed seeds | `completed_cohort_prevalence.tsv`: candidate, completed-cohort truth, arm estimate, denominator, absolute error or undefined | PENDING. Simulation repeats are not biological replicates. Compare same completed cohort; full-cohort difference is selection bias, not corrected by metadata. |
| Figure 4A: what researcher preparation is required? | A03 equivalent native/RNA conversion plus limited nf-core input-channel task | `reference_comparison.tsv` (operations/code/configuration requirements); `nf_input_mapping.tsv` (same four real libraries, emitted IDs/layout/lane groups); reviewed reference and harness hashes | PENDING. Metadata-only pipeline task is explicitly limited. Ordinary sample IDs can preserve cell identity. No human usability study. |
| Figure 4B: what does the same post-processing endpoint cost? | A03 native+RNA -> objects/evidence/QC/HTML in both implementations | Three fresh-process rows per arm, warm-cache policy, input/output bytes, wall time, peak RSS and all failed outcomes | PENDING. Never compare import-only cost with alignment/detection; no forced speed advantage. |
| Figure 4C/D: what constrains scaling and reproduction? | A06 five non-Cartesian synthetic configurations; A07 existing installed artifacts | `scaling_runs.tsv` with actual native bytes/file counts/read-ID lengths/shape/nnz/RNA size/threads/env; linked P0 installed evidence | Scaling PENDING; installation evidence already recorded. Synthetic workloads remain labeled. Local wheel/conda validation is not PyPI/Bioconda publication. |
| Table 1: which tasks are supported or evaluated? | A03 source-linked five-tool assessment in study_design.md; expand only with measured evidence | Pinned CIRI3, nf-core/circrna, CirComPara2, circtools 2.0.4 docs and N0 comparison; empirical task outcomes later | Documentation assessment available. Use documented / tested / requires code / not documented / not assessed / not applicable, without turning absence of documentation into inability. |
| Table 2: what supports each claim? | A01/A07 existing records and future A02–A06 | `claim_evidence.tsv`, `data_provenance.tsv`, completed output records when available | Evidence map available; new experiment rows explicitly PENDING. Historical analysis versions and eventual installation release remain distinct. |

Every new result table must retain source/implementation/helper hashes, argv,
environment, run status and error information. Integer equality and source axes
are checked before plotting. Display aggregation does not overwrite source data.
Source tables contain all eligible cells/candidates/configurations; graphical
summaries may aggregate only by the prespecified rules. No output is selected
because it strengthens a desired conclusion.

Proposed supplemental groups are descriptive working groups only: original
run/selection provenance; full perturbations, thresholds and negative outcomes;
environments/commands/hashes/acceptance and complete evidence maps. Final S1/S2/S3
numbering will be reconciled during manuscript review. Do not overwrite existing
M0 `supplement_v2/table_s1_*`, `table_s2_*`, `table_s3_*` or figures.

The pilot's immediate contribution is Figure 2A/B/D and the real-input portion
of Figure 4. Figure 2C requires the separate CNV missingness review. Figures 3
and 4 scaling/comparative measurements remain later reviewed executions. Preserve
Figure 1 throughout; no figure depends on a positive biological association.
