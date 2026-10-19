# A05 results for review — 2026-10-19

**A05 COMPLETE WITH LIMITATIONS. CNV eligibility: NOT VERIFIED.** The registered RNA/circ descriptive outputs are complete. The CNV panel is withheld; the proposed Figure 2C substitute displays the 23 original IMR90 cells with RNA/circ metrics in their stored order. No `cnv_display_bins.tsv` or CNV heatmap was produced. These are factual review drafts, not manuscript revisions.

## Source identity and reused evidence

The branch `analysis/cmpb-a05-descriptive` starts exactly at accepted A03/A04 commit `02ee3d4acfed63c5b2f02c4baef6d22c822d4dac`. Accepted A02 is `307758f57a39c009393d03fc0f69e4638b2498c6`; production remains `e6843a38469d49f96d44758a4e37893ebbe2734b`. A02’s four genuine cells and A04’s simulations supply prior technical evidence only; none of their measured or simulated values replaced the historical originals here. The 5,159,619-row accepted A04 prevalence gzip was rehashed read-only: `8163178301f702fb4217fc0e48d4390ed82aa39244a8d1efd2cc681ad9152b8c`; no grid was rerun.

Both H5MU byte hashes match archive commit `c99cddaae2314d3af7d5232150c5b86bb6b52e25`, including a fresh hash of each Git archive blob. They are unchanged after analysis.

| Dataset | Original object SHA-256 |
| --- | --- |
| Smart-seq3 | `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339` |
| IMR90 | `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8` |

## Stored measurements, analyzed separately

| Dataset | Stored rows / recorded circ members | Candidates | Positive entries | Support | Host-annotated | Candidate/support medians per stored row |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Smart-seq3 | 192 / 191 | 2503 | 2659 | 5711 | 2379 / 2503 | 12.0 / 22.5 |
| IMR90 | 23 / 23 | 2443 | 2914 | 7422 | 2429 / 2443 | 119.0 / 317.0 |

Smart-seq3 `DIYHEK_192` is physically stored in both RNA and circ axes, labelled `rna_only`, and has a padded circ zero row. The full stored-row distribution includes it (N=192); the recorded circ-member and matched-scatter denominator is 191. Its original circ completion remains unknown. The support median among those 191 members is 23, while the 192-row median is 22.5. IMR90 has 23 shared RNA/circ/CNV identities. All 215 historical per-cell circ statuses remain unknown; none is certified completed-zero.

All RNA matrices have 63,187 stored features. Smart-seq3 RNA total signal is 42,535,540 with 1,695,807 positive entries; IMR90 has 67,525,881 with 376,158 positive entries. These are the historical lightweight simple-overlap profile counts, not independent transcript quantification, TPM, or verified sequencing depth. Original reference digests remain unrecovered.

| Dataset | Metric | N | Q25 | Median | Q75 | Range |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Smart-seq3 | circ_positive_candidates | 192 | 7.0 | 12.0 | 19.0 | 0.0–53.0 |
| Smart-seq3 | circ_support | 192 | 11.0 | 22.5 | 43.25 | 0.0–197.0 |
| Smart-seq3 | rna_overlap_signal | 192 | 112611.0 | 205045.5 | 322025.25 | 13725.0–736965.0 |
| Smart-seq3 | rna_positive_features | 192 | 8355.0 | 9046.0 | 9507.75 | 2393.0–10933.0 |
| IMR90 | circ_positive_candidates | 23 | 108.0 | 119.0 | 134.5 | 87.0–244.0 |
| IMR90 | circ_support | 23 | 271.0 | 317.0 | 344.5 | 188.0–509.0 |
| IMR90 | rna_overlap_signal | 23 | 2799210.5 | 2913778.0 | 3112403.0 | 2438226.0–3470571.0 |
| IMR90 | rna_positive_features | 23 | 15941.5 | 16171.0 | 16687.0 | 15207.0–18352.0 |

Quantiles use linear interpolation with all zero rows retained. All per-candidate distributions and original annotation fields, including multiple host genes and blanks, are in `descriptive_candidate_metrics.tsv` and `distribution_summary.tsv`. No gene or candidate was selected for a positive example.

## Prespecified threshold sensitivity

A applies the threshold independently to each original cell–candidate entry. B applies it to the original pooled candidate sum, retaining all original entries of eligible candidates. Each threshold starts from original X. All physical rows remain in denominators, including threshold-induced zero rows.

| Dataset | Mode | Threshold | Retained candidates | Positive entries | Retained support | Zero stored rows |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Smart-seq3 | A: entry | 1 | 2503 | 2659 | 5711 | 1 |
| Smart-seq3 | A: entry | 2 | 1321 | 1383 | 4435 | 12 |
| Smart-seq3 | A: entry | 3 | 706 | 738 | 3145 | 38 |
| Smart-seq3 | A: entry | 5 | 219 | 224 | 1430 | 99 |
| Smart-seq3 | B: pooled | 1 | 2503 | 2659 | 5711 | 1 |
| Smart-seq3 | B: pooled | 2 | 1332 | 1488 | 4540 | 10 |
| Smart-seq3 | B: pooled | 3 | 722 | 867 | 3320 | 28 |
| Smart-seq3 | B: pooled | 5 | 243 | 370 | 1716 | 65 |
| IMR90 | A: entry | 1 | 2443 | 2914 | 7422 | 0 |
| IMR90 | A: entry | 2 | 1116 | 1344 | 5852 | 0 |
| IMR90 | A: entry | 3 | 509 | 645 | 4454 | 0 |
| IMR90 | A: entry | 5 | 118 | 195 | 2965 | 0 |
| IMR90 | B: pooled | 1 | 2443 | 2914 | 7422 | 0 |
| IMR90 | B: pooled | 2 | 1164 | 1635 | 6143 | 0 |
| IMR90 | B: pooled | 3 | 588 | 1018 | 4991 | 0 |
| IMR90 | B: pooled | 5 | 184 | 509 | 3640 | 0 |

These threshold summaries do not select or validate a biological cutoff. Stored-cell and recorded-member fractions are explicitly separate from undefined completed-cohort prevalence.

## CNV gate and proposed replacement

All 60,607 bin coordinates and 1,393,961 cell–bin values match each retrieved GEO 50-kb table after exact ID alignment and the documented float32 conversion for mappabilitynorm. No numerical or axis contradiction was found. The original state array contains 115,093 zeros; its mappabilitynorm layer contains 115,834 zeros. Neither has explicit NaN/Inf, and neither public table contains the tested missing-value tokens. The object has no explicit missing/uncallable mask. None of these facts certifies callability.

The missing evidence is the exact cell/bin callability mask or authoritative export rule distinguishing uncallable/missing bins from true numeric zero, tied to these coordinates and the original modified blacklist. See `cnv_eligibility_decision.md`. The fixed 1-Mb/50%-valid-base rule was not attempted or altered.

Review `figure2c_rna_circ_only_alternative.svg` or its PNG: four columns show circ candidate count, circ support, RNA positive-feature count and RNA overlap-profile total for all 23 cells. Separate dataset-specific scatter, full empirical distributions and threshold drafts accompany it. No correlation coefficient, fit, p-value, treatment/clonality analysis or host-gene association test was computed.

## Verification, execution and limitations

Independent raw-CSR Python integer oracles verified every cell/candidate/threshold against scipy-based summaries. All 10,874 plotted observations in 7 SVG drafts have source-table checks; raster previews were visually inspected after fixing a renderer compatibility issue. Technical assertion counts are not biological replicates. `preservation_record.json` verifies 654 protected files, frozen refs, both originals and 87 installed production Python files.

The reused Python 3.10.22 environment contains AnnData 0.11.4, MuData 0.3.10, NumPy 2.2.6, SciPy 1.15.3, pandas 2.3.3 and h5py 3.16.0. No dependencies were installed. `execution_record.json` contains exact stage timings, CPU accounting, sampled memory peaks, disk accounting and the fixed one-hour deadline. All science ran within verified two-CPU affinity, 4-GiB cgroup memory and zero-swap limits, with automatic time/disk/process-tree stops. Kernel CPU quota was not delegated; inherited affinity plus all-thread monitoring provided the two-core bound, as in accepted A03/A04.

The initial metadata inventory exited 1 only when printing a NumPy integer to JSON after writing its inventory; the corrected stage passed and the failed record is retained. Initial PNGs silently omitted nested-title circles in ImageMagick; self-closing SVG circles corrected this, and reviewed figures were regenerated and revalidated in a new directory. Scientific tables were unchanged. Expected resource-probe failures are kept separately from these administrative issues.

Unknown original per-cell completion, original aligner/detector execution provenance and original reference/GTF digests remain unknown. Stored junction-read support is not independent validation of circularization, PCR-deduplicated molecule abundance, biological prevalence or detection sensitivity. Cells are not donor-independent replication. Assays/material differ, so no cross-dataset yield/superiority inference is supported. These limitations also apply to visually apparent scatter patterns.

Existing Figure 1, circulated manuscript, BioRender and prior analysis branches were preserved. No sequencing, alignment, detector, grid, A06, release, PR or merge was started. Final figure choice and manuscript interpretation remain for owner review.
