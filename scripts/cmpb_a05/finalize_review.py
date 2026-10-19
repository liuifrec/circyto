"""Assemble small A05 review records from verified results, without matrix analysis."""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

from inspect_sources import digest
from run_stage import footprint


def read(p):
    with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    root=Path(sys.argv[1]);repo=Path(__file__).resolve().parents[2]
    source=root/'results_review';out=repo/'manuscript/cmpb/a05_descriptive'
    assert json.loads((source/'validation_record.json').read_text())['passed']
    replacements={str(root):'${A05_ROOT}',str(repo):'${A05_WORKTREE}',
        '/tmp/circyto-import/load_work':'${ORIGINAL_OBJECT_ROOT}',
        '/home/liuyuchen/data/circyto_cmpb_a02_20261015':'${A02_ROOT}',
        '/home/liuyuchen/Github/circyto':'${ORIGINAL_CHECKOUT}',
        '/user/ifrec/liuyuchen/circyto_redo':'${HISTORICAL_PROJECT_ROOT}',
        '/sys/fs/cgroup/user.slice/user-1001.slice/user@1001.service':'${USER_CGROUP}'}
    def sanitize(s):
        for a,b in replacements.items():s=s.replace(a,b)
        return s
    for p in source.iterdir():
        assert p.suffix in ['.tsv','.json','.svg','.png'],p
        if p.suffix=='.png':shutil.copyfile(p,out/p.name)
        else:(out/p.name).write_text(sanitize(p.read_text()),encoding='utf-8')
    controls=out/'execution_details';controls.mkdir(exist_ok=True)
    for p in (root/'records').glob('*.json'):
        if p.name in ['protected_before.json','refs_before.json']:continue
        (controls/p.name).write_text(sanitize(p.read_text()))
    for p in (root/'logs').glob('*.log'):(controls/p.name).write_text(sanitize(p.read_text()))
    auth=out/'authorization';auth.mkdir(exist_ok=True)
    for p in (root/'inputs').glob('*.md'):shutil.copyfile(p,auth/p.name)
    cs=read(out/'descriptive_cell_metrics.tsv');vs=read(out/'descriptive_candidate_metrics.tsv');ts=read(out/'threshold_sensitivity.tsv')
    summary=read(out/'distribution_summary.tsv');numeric=read(out/'cnv_cell_missingness_audit.tsv')
    validation=json.loads((out/'validation_record.json').read_text())
    manifest=json.loads((out/'figure_manifest.json').read_text())
    provenance=json.loads((out/'analysis_provenance.json').read_text())
    report=['# A05 results for review — 2026-10-19','',
        '**A05 COMPLETE WITH LIMITATIONS. CNV eligibility: NOT VERIFIED.** The registered RNA/circ descriptive outputs are complete. The CNV panel is withheld; the proposed Figure 2C substitute displays the 23 original IMR90 cells with RNA/circ metrics in their stored order. No `cnv_display_bins.tsv` or CNV heatmap was produced. These are factual review drafts, not manuscript revisions.','',
        '## Source identity and reused evidence','',
        'The branch `analysis/cmpb-a05-descriptive` starts exactly at accepted A03/A04 commit `02ee3d4acfed63c5b2f02c4baef6d22c822d4dac`. Accepted A02 is `307758f57a39c009393d03fc0f69e4638b2498c6`; production remains `e6843a38469d49f96d44758a4e37893ebbe2734b`. A02’s four genuine cells and A04’s simulations supply prior technical evidence only; none of their measured or simulated values replaced the historical originals here. The 5,159,619-row accepted A04 prevalence gzip was rehashed read-only: `8163178301f702fb4217fc0e48d4390ed82aa39244a8d1efd2cc681ad9152b8c`; no grid was rerun.','',
        'Both H5MU byte hashes match archive commit `c99cddaae2314d3af7d5232150c5b86bb6b52e25`, including a fresh hash of each Git archive blob. They are unchanged after analysis.','',
        '| Dataset | Original object SHA-256 |','| --- | --- |']
    for ds,sha in provenance['source_hashes'].items():report.append(f'| {ds} | `{sha}` |')
    report+=['','## Stored measurements, analyzed separately','',
        '| Dataset | Stored rows / recorded circ members | Candidates | Positive entries | Support | Host-annotated | Candidate/support medians per stored row |','| --- | ---: | ---: | ---: | ---: | ---: | --- |']
    for ds in ['Smart-seq3','IMR90']:
        cc=[r for r in cs if r['dataset']==ds];vv=[r for r in vs if r['dataset']==ds]
        q={r['metric']:r for r in summary if r['dataset']==ds and r['subset']=='all_stored_circ_rows'}
        members=sum(r['recorded_circ_member']=='True' for r in cc);ann=sum(bool(r['host_gene'].strip()) for r in vv)
        nnz=sum(int(r['circ_positive_candidates']) for r in cc);support=sum(int(r['circ_support']) for r in cc)
        report.append(f"| {ds} | {len(cc)} / {members} | {len(vv)} | {nnz} | {support} | {ann} / {len(vv)} | {q['circ_positive_candidates']['median']} / {q['circ_support']['median']} |")
    report+=['','Smart-seq3 `DIYHEK_192` is physically stored in both RNA and circ axes, labelled `rna_only`, and has a padded circ zero row. The full stored-row distribution includes it (N=192); the recorded circ-member and matched-scatter denominator is 191. Its original circ completion remains unknown. The support median among those 191 members is 23, while the 192-row median is 22.5. IMR90 has 23 shared RNA/circ/CNV identities. All 215 historical per-cell circ statuses remain unknown; none is certified completed-zero.','',
        'All RNA matrices have 63,187 stored features. Smart-seq3 RNA total signal is 42,535,540 with 1,695,807 positive entries; IMR90 has 67,525,881 with 376,158 positive entries. These are the historical lightweight simple-overlap profile counts, not independent transcript quantification, TPM, or verified sequencing depth. Original reference digests remain unrecovered.','',
        '| Dataset | Metric | N | Q25 | Median | Q75 | Range |','| --- | --- | ---: | ---: | ---: | ---: | --- |']
    for r in summary:
        if r['subset'] in ['all_stored_circ_rows','all_stored_rna_rows']:
            report.append(f"| {r['dataset']} | {r['metric']} | {r['n_total']} | {r['q25']} | {r['median']} | {r['q75']} | {r['minimum']}–{r['maximum']} |")
    report+=['','Quantiles use linear interpolation with all zero rows retained. All per-candidate distributions and original annotation fields, including multiple host genes and blanks, are in `descriptive_candidate_metrics.tsv` and `distribution_summary.tsv`. No gene or candidate was selected for a positive example.','',
        '## Prespecified threshold sensitivity','',
        'A applies the threshold independently to each original cell–candidate entry. B applies it to the original pooled candidate sum, retaining all original entries of eligible candidates. Each threshold starts from original X. All physical rows remain in denominators, including threshold-induced zero rows.','',
        '| Dataset | Mode | Threshold | Retained candidates | Positive entries | Retained support | Zero stored rows |','| --- | --- | ---: | ---: | ---: | ---: | ---: |']
    for r in ts:report.append(f"| {r['dataset']} | {'A: entry' if r['threshold_mode']=='individual_entry' else 'B: pooled'} | {r['threshold']} | {r['retained_candidates']} | {r['positive_cell_candidate_entries']} | {r['retained_support']} | {r['stored_zero_rows']} |")
    report+=['','These threshold summaries do not select or validate a biological cutoff. Stored-cell and recorded-member fractions are explicitly separate from undefined completed-cohort prevalence.','',
        '## CNV gate and proposed replacement','',
        'All 60,607 bin coordinates and 1,393,961 cell–bin values match each retrieved GEO 50-kb table after exact ID alignment and the documented float32 conversion for mappabilitynorm. No numerical or axis contradiction was found. The original state array contains 115,093 zeros; its mappabilitynorm layer contains 115,834 zeros. Neither has explicit NaN/Inf, and neither public table contains the tested missing-value tokens. The object has no explicit missing/uncallable mask. None of these facts certifies callability.','',
        'The missing evidence is the exact cell/bin callability mask or authoritative export rule distinguishing uncallable/missing bins from true numeric zero, tied to these coordinates and the original modified blacklist. See `cnv_eligibility_decision.md`. The fixed 1-Mb/50%-valid-base rule was not attempted or altered.','',
        'Review `figure2c_rna_circ_only_alternative.svg` or its PNG: four columns show circ candidate count, circ support, RNA positive-feature count and RNA overlap-profile total for all 23 cells. Separate dataset-specific scatter, full empirical distributions and threshold drafts accompany it. No correlation coefficient, fit, p-value, treatment/clonality analysis or host-gene association test was computed.','',
        '## Verification, execution and limitations','',
        f"Independent raw-CSR Python integer oracles verified every cell/candidate/threshold against scipy-based summaries. All {validation['plotted_observations_checked']:,} plotted observations in {validation['figures_checked']} SVG drafts have source-table checks; raster previews were visually inspected after fixing a renderer compatibility issue. Technical assertion counts are not biological replicates. `preservation_record.json` verifies 654 protected files, frozen refs, both originals and 87 installed production Python files.",'',
        'The reused Python 3.10.22 environment contains AnnData 0.11.4, MuData 0.3.10, NumPy 2.2.6, SciPy 1.15.3, pandas 2.3.3 and h5py 3.16.0. No dependencies were installed. `execution_record.json` contains exact stage timings, CPU accounting, sampled memory peaks, disk accounting and the fixed one-hour deadline. All science ran within verified two-CPU affinity, 4-GiB cgroup memory and zero-swap limits, with automatic time/disk/process-tree stops. Kernel CPU quota was not delegated; inherited affinity plus all-thread monitoring provided the two-core bound, as in accepted A03/A04.','',
        'The initial metadata inventory exited 1 only when printing a NumPy integer to JSON after writing its inventory; the corrected stage passed and the failed record is retained. Initial PNGs silently omitted nested-title circles in ImageMagick; self-closing SVG circles corrected this, and reviewed figures were regenerated and revalidated in a new directory. Scientific tables were unchanged. Expected resource-probe failures are kept separately from these administrative issues.','',
        'Unknown original per-cell completion, original aligner/detector execution provenance and original reference/GTF digests remain unknown. Stored junction-read support is not independent validation of circularization, PCR-deduplicated molecule abundance, biological prevalence or detection sensitivity. Cells are not donor-independent replication. Assays/material differ, so no cross-dataset yield/superiority inference is supported. These limitations also apply to visually apparent scatter patterns.','',
        'Existing Figure 1, circulated manuscript, BioRender and prior analysis branches were preserved. No sequencing, alignment, detector, grid, A06, release, PR or merge was started. Final figure choice and manuscript interpretation remain for owner review.']
    (out/'a05_results_for_review.md').write_text('\n'.join(report)+'\n')
    (out/'cnv_eligibility_decision.md').write_text('''# IMR90 CNV eligibility: NOT VERIFIED

The mandatory missingness/callability gate does not pass. No CNV heatmap or 1-Mb aggregation is produced, and `cnv_display_bins.tsv` does not exist. No source-value conflict was found, so this is NOT VERIFIED rather than FAIL.

The read-only original `/mod/cnv` contains a 23 × 60,607 int16 X and one float32 `mappabilitynorm` layer, with cell IDs, bin IDs, coordinates and importer metadata. All fields are listed in `cnv_source_inventory.tsv`; per-cell numerical missingness and chromosome coordinate summaries are separate TSVs. There is no explicit callability/missingness mask or missing-labelled field. Source bins are 50,000 bases, with stored end-start widths and IDs. They remain unmodified.

The two public 50-kb source tables were acquired from [GEO GSE278958](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278958) on 2026-10-19 Japan time. All 23 IDs, 60,607 coordinates, and 1,393,961 values per table match the original after exact ID mapping (state integer parsing and mappabilitynorm float32 import conversion). File names, URLs, access times, sizes and SHA-256s are recorded in `cnv_source_inventory.tsv`. The original server-local files remain unavailable; matching public bytes at the value level does not recover their historical file digests.

X has 115,093 numerical zeros and mappabilitynorm has 115,834. There are no explicit NaN/Inf values or tested public missing tokens. The source `<N>-somy` tokens parse to integers, but neither that conversion nor successful H5MU import establishes original callability. Zero is not assigned to missing, callable, diploid or deletion here.

The [source paper's CNV Methods](https://www.nature.com/articles/s41467-025-64688-1) describes AneuFinder 1.2.1, dnacopy, mappability correction, an ENCODE blacklist and merged G1 reads. Its [pinned Fig. 6 code](https://github.com/mcbmieu/scRR-seq/blob/0910c9a30ce1b7e8a56cb97ac07dcb63a585c626/MainFigures/Fig6/scripts/Codes_Fig6_part2.R) points to `hg38-blacklist.v2.mod.bed`, but the inspected files do not contain that exact modified mask, the per-cell model/bin validity flags, or the export rule linking excluded/model bins to these zero tokens. This is a bounded finding in the original H5MU, importer, GEO metadata/tables, paper and pinned Fig. 6 scripts; it does not establish absence from all author storage. The study's RSEM/GENCODE RNA methods do not describe the historical circyto simple-overlap matrices used in A05.

Required evidence to reconsider: an authoritative missing/uncallable definition and compatible numerical semantics for these exact exported states, plus its bin/cell mapping (or original model validity mask), modified blacklist identity and export code/parameters sufficient to distinguish genuine zero from missing/uncallable. A fresh standard blacklist or a rule inferred from observed zeros is not a substitute.

Proposed Figure 2C replacement: `figure2c_rna_circ_only_alternative.svg`/`.png`, all 23 cells in original order, showing stored circ candidate count/support and RNA positive-feature count/overlap-profile signal. No source chromosome selection, CNV diagnosis, clone, radiation, causal or association claim is made. The registered 1-Mb median and >=50% valid-covered-base requirement remain unchanged and unexecuted.
''')
    (out/'methods_results_data_note.md').write_text('''# Factual Methods/Results data note for review

Inputs were the two byte-identified original H5MU files at archive c99cddaae2314d3af7d5232150c5b86bb6b52e25. h5py opened them in read-only mode; AnnData read_elem read their stored CSR matrices and axis/annotation tables. Cell and feature IDs were unique. Sparse scipy reductions computed per-row/per-column positive counts, integer support sums and candidate maxima, with no full genomic matrix densification. An independent verifier read the raw CSR indptr/indices/data and used Python integer dictionaries to check every result. Before/after and archive-blob SHA-256s matched.

All stored cells and candidates were retained. The RNA-only padded Smart-seq3 row was distinguished from recorded shared membership and from unknown historical circ completion. Candidate fractions explicitly use either physical stored circ rows or recorded circ members; completed-cohort prevalence is left undefined. Original coordinate/strand and all host-gene annotation fields were copied without reannotation or imputation. Historical RNA is a lightweight simple-overlap profile; original alignment/reference/GTF identities were not recovered.

Individual-entry and pooled-candidate thresholds 1, 2, 3 and 5 were independently applied to original X, with separate output tables and unchanged row denominators. Summaries use linear-interpolated quartiles. Scatter drafts show all recorded matched identities (Smart-seq3 191, IMR90 23), with no correlation, regression or hypothesis test. Distribution drafts retain every observation, including zero rows, with empirical cumulative fractions. No cross-dataset numerical comparison is interpreted biologically.

The audited circ totals were reproduced: Smart-seq3 192 stored rows × 2,503 candidates, 2,659 positive entries, 5,711 support, 2,379 host-annotated candidates; IMR90 23 × 2,443, 2,914 positive entries, 7,422 support, 2,429 annotated. Median positive candidates/support per stored row were 12/22.5 and 119/317 respectively. These are descriptions of separate stored datasets.

The CNV source-value inspection found exact coordinate/value concordance with the public GEO 50-kb tables and no explicit NaN/Inf values, but could not recover an authoritative callability/missingness definition for exported zeros. The gate is NOT VERIFIED. No CNV aggregation or heatmap was made. The RNA/circ-only Figure 2C alternative uses all 23 original cells in source order.

These results describe stored measurements and technical preservation only. They do not establish true prevalence, biological replication, circularization truth, detection sensitivity/specificity, donor-independent effects, causal relationships, radiation effects or new clones. A02 genuine rerun evidence and A04 simulated evidence were reused as prior milestones, not substituted for original A05 matrices. No A06 work was performed.
''')
    (out/'figure_data_dictionary.md').write_text('''# A05 source tables and figure dictionary

All TSVs are UTF-8 with a header, tab delimiter and LF newlines. Blank numerical fields mean unavailable/undefined, not zero. Literal IDs and blank annotation strings are preserved. Boolean fields are True/False. Indices/orders are zero based. `dataset`, `source_object_sha256` and `helper_commit` identify lineage. These are historical stored measurements, not A02 rerun or A04 simulation values.

| Table | Observation/key | Field meanings and denominators |
| --- | --- | --- |
| descriptive_cell_metrics.tsv | 215 rows; dataset + exact cell_id | `source_cell_order`: original global obs order; `circ_source_row`/`rna_source_row`: original modality positions. Physical row presence, recorded membership and historical processing status are distinct. `circ_positive_candidates` counts X>0; `circ_support` sums X in stored junction-read units. `rna_overlap_signal` sums historical simple-overlap counts; `rna_positive_features` counts positive RNA features. `matched_rna_circ` requires shared membership and exact ID intersection. `circ_zero_support_row` does not certify `completed_zero`. All completion fields stay unknown. |
| modality_membership_audit.tsv | same 215 identities/order | Physical circ/RNA/CNV row presence; recorded circ/RNA membership; matched flag, source indices and unknown completion. Smart-seq3 DIYHEK_192 is physically present, rna_only, unmatched, circ zero. |
| descriptive_candidate_metrics.tsv | 4,946 rows; dataset + candidate_id | Original var order; original circ_id/chrom/start/end/strand and every host_gene* field, including multi annotations, are retained. `positive_stored_cells`, `pooled_support`, `maximum_per_cell_support` are original X reductions. `stored_cell_fraction` uses 192/23 physical rows; `recorded_circ_member_fraction` uses 191/23 recorded members. Both denominators are explicit. Completed denominator/prevalence are blank because historical status is unknown. Host annotation provenance is original var, with original GTF digest unverified. |
| circ_nonzero_entries.tsv | 5,573 positive cell/candidate records | Exact original positive support, biological IDs and source axis positions; implicit zeros reconstructed from complete cell/candidate axes, not interpreted as completed assays. This small circ-only sparse table is not a CNV per-bin export. |
| threshold_sensitivity.tsv | 16 rows; dataset + mode + threshold | `individual_entry`: retain entries with original X>=t. `pooled_candidate`: original column sum>=t, retain all original entries in eligible columns. Counts, support, zero rows and medians use every physical row. No sequential filtering. |
| threshold_cell_metrics.tsv | 1,720 rows; dataset + mode + threshold + cell_id | Every original stored cell at each threshold; remaining positive candidates/support and membership. |
| distribution_summary.tsv | dataset + subset + metric | `n_total/n_finite/n_missing` describe the explicit subset. Minimum, Q25, median, Q75, maximum; linear interpolation. Candidate fractions are stored-cell fractions, not verified completed-cohort prevalence. |
| matrix_audit.tsv | dataset + modality | Cells-by-features shape, CSR stored entries vs numerical nonzeros, explicit zeros, integer total signal, dtype and uniqueness/missing-value checks. |
| registration_checks.tsv | dataset | Exact original audited tuple: rows, candidates, positive entries, support, host-annotated candidates. |
| cnv_source_inventory.tsv | field | Original fields, source acquisition hashes/URLs, public equality results and gate. No original full per-bin arrays. |
| cnv_numeric_encoding_audit.tsv | source + encoding | Numeric integer and public somy token frequencies only; no callability assignment. |
| cnv_cell_missingness_audit.tsv | all 23 source cells | Numerical zero/NaN/Inf and discordant-zero counts only. They are NOT valid/invalid counts or masks. |
| cnv_coordinate_audit.tsv | stored chromosome order | Bin counts, coordinate extent, width and gap/overlap/order checks; no chromosome selection or CNV aggregation. |
| figure_observations.tsv | figure + panel + display_order (within mode for threshold plots) | Every plotted x/y value, identity, source table/hash and helper/renderer commits. ECDF rows use stable value then original-order sorting and rank/N; y is fraction. Figure 2C uses raw metric x and source cell ordinal y. Threshold observation_id includes mode and threshold. |

`figure_manifest.json` maps every SVG to source tables, source object, helper and renderer. Every SVG embeds provenance and exact units/denominators in its labels/metadata. PNGs are previews of those SVGs; `png_render_record.json` binds source/render hashes, dimensions and the installed renderer. The revised SVG uses self-closing circle nodes for renderer compatibility. Source rows/IDs and all zero members remain unchanged.

| Draft | Panels / axes and population |
| --- | --- |
| figure2c_rna_circ_only_alternative | All 23 IMR90 cells in unchanged stored order; x columns: positive circ candidates, circ support units, RNA positive features, RNA overlap-profile counts. No CNV. |
| smartseq3_rna_circ_scatter_draft; imr90_rna_circ_scatter_draft | A: RNA overlap total x, circ support y. B: RNA positive features x, positive circ candidates y. N=191 recorded shared Smart-seq3 or 23 IMR90, separately. No fitted relation or coefficient. |
| smartseq3_distributions_draft; imr90_distributions_draft | ECDF y=rank/N; x: per-row positive circ candidates, per-row circ support, positive stored cells/candidate, pooled support/candidate. All 192/23 stored rows and 2,503/2,443 candidates. |
| smartseq3_threshold_sensitivity_draft; imr90_threshold_sensitivity_draft | x=registered threshold 1,2,3,5; y=retained candidate count or total support. Both threshold definitions, separate datasets. Lines only connect tested thresholds. |

`source_object_inventory.json` records the original axes, metadata and HDF5 field inventory. Published metadata replaces workstation path prefixes with documented aliases; local control logs retain exact paths. Original biological IDs, numerical data and source object hashes are unchanged. `validation_record.json`, `preservation_record.json` and `execution_record.json` define check/resource scope and limitations.
''')
    renderer=manifest[0]['renderer_commit']
    (out/'README.md').write_text(f'''# CIRCYTO CMPB A05 review packet

Start with `a05_results_for_review.md`, `cnv_eligibility_decision.md` and `figure2c_rna_circ_only_alternative.png`. All source tables and draft SVG/PNG figures are included. This packet stops after A05. The CNV gate is NOT VERIFIED; no heatmap or fake empty display table is included. Publication figure choice and manuscript prose await review.

Scientific analysis helper: `{provenance['helper_commit']}`. Figure renderer helper: `{renderer}`. The ZIP's top-level HANDOFF_README.txt records the final committed packet SHA and branch. `SHA256SUMS` in the ZIP covers every other archive member; the outer checksum file covers the ZIP. `inputs_outputs_sha256.json` binds original source hashes, outputs and scripts, explicitly excluding its own hash to avoid self-reference.

To recreate only the vector drafts from the included committed tables, unpack the ZIP, run from its root, and use a new empty destination (Python standard library only):

```bash
python3 scripts/cmpb_a05/plot_drafts.py \\
  --tables manuscript/cmpb/a05_descriptive \\
  --out /tmp/a05-recreated-figures \\
  --renderer-commit {renderer}
```

The SVG/figure-observation output should match the packet byte-for-byte with the same script and pinned renderer commit. Optional PNG conversion uses existing ImageMagick 6.9.11-60 via `python3 scripts/cmpb_a05/render_png.py /tmp/a05-recreated-figures`; raster pixels may vary with fonts/renderer version on other machines. No dependency installation is required to recreate SVGs.

Full scientific recomputation requires the two original checksum-matched H5MU files (not in this ZIP) and the accepted Python/AnnData environment. `analyze.py --original ORIGINAL_LOAD_WORK --out NEW_OUTPUT`, `cnv_gate.py A05_ROOT ORIGINAL_LOAD_WORK NEW_OUTPUT`, and `verify.py --original ORIGINAL_LOAD_WORK --tables NEW_OUTPUT` document the calculation/check interfaces. The exact executed argv, control properties, input URLs/hashes and environment versions are included in execution_details, execution_record.json and source_object_inventory.json. Running resource controls on a new machine requires verifying equivalent limits first; no uncontrolled rerun is recommended by this packet.

Only small A05 outputs/scripts are packaged. The original H5MU, downloaded full CNV source tables, raw reads, native A02 outputs, full A04 prevalence archive, environments, executables, indices and secrets are excluded. Public CNV table URLs and hashes are provided for provenance. The unchanged A02/A03/A04 evidence is linked by accepted commits and recorded preservation checks.
''')
    (out/'visual_review.json').write_text(json.dumps(dict(review='Manual inspection of PNG previews after self-closing-circle compatibility correction',
        decision='Reviewed raster markers, labels, ranges, source ordering and no clipping; numerical SVG/source checks passed.',
        limitation='Draft figures only; publication layout/selection awaits review.',corrected_initial_issue='Installed ImageMagick ignored circles with nested title nodes; first PNG set retained locally and excluded from handoff.'),indent=2)+'\n')
    print(json.dumps(dict(review_directory=str(out),files=sum(p.is_file() for p in out.rglob('*')))))


if __name__=='__main__':main()
