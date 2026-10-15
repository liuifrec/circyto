# Original detection provenance: bounded recovery

Audit date: **2026-10-15**. Comment 3 is **partially addressed**. This audit
recovers original metadata embedded in the two checksum-matched objects and
separates it from historical script reconstruction and later figure regeneration.
It does not establish a complete original FASTQ-to-object chain. No raw analysis,
detector run, figure regeneration, input repair or new biological analysis was
performed for this revision.

The field-level [provenance TSV](original_detection_provenance.tsv) records every
value, classification, exact Git/object identity, source field and missing record.
The [cell-identity TSV](original_cell_identity.tsv) gives all **192 Smart-seq3 +
23 IMR90** stored circ-modality rows in original order, with available GSM,
RNA/DNA and canonical mappings. These are final-object inclusions, not an invented
original exclusion ledger. The exact HAP1 batch10 list remains unknown.

## Evidence classes and identities

- **D — directly recorded:** a statement in a cited record, or a value actually
  stored in a checksum-matched object. This does not certify that a mutable
  summary was a pristine launch record or independently prove its contents.
- **R — reconstructed from a contemporaneous versioned script but execution
  unconfirmed:** code existed in the relevant period; no run-bound revision or
  command log proves it was the deployed code or that defaults were used.
- **U — unknown:** no adequate record recovered. Neither today's installed tools
  nor a filename/default fills an unknown historical field.

| Source | Exact identity |
| --- | --- |
| Original object archive O0 | `c99cddaae2314d3af7d5232150c5b86bb6b52e25` |
| Smart-seq3 O-SS | O0:`load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu`; SHA-256 `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339` |
| IMR90 O-IM | O0:`load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu`; SHA-256 `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8` |
| Historical May source S-May | `9e028d8dd6b054bbf34fb77debb4ec095916c359` (2026-05-20) |
| Historical June source S-June | `777dc3406e9793703aba79bac7dfa2473911d739` (2026-06-01) |
| Server wrapper lead S-Pilot | `6e08139221537d95cf0b68b53c16dc1a161939dc` (2026-05-21): `scripts/run_scrr_imr90_full_stack.sh`, `scripts/run_scrr_hap1_full_stack.sh` |
| Frozen documented implementation C0 | `44697355bcab1c525ca7ef9b130e2ad0094d9e1b`: `docs/validated_workflows_summary.md`, `docs/current_project_status.md`, `docs/ramda_full_stack_server_run.md`, `testdata/public_datasets/human_ramda_candidate_runs.tsv` |
| Existing manuscript evidence M0 | `c29831c3f6e7a31b708cce493c0be8036176077a` |

Both local originals in `/tmp/circyto-import/load_work/` matched O-SS/O-IM before
inspection. Only HDF5 `r` mode and public `anndata.io.read_elem` were used to
extract stored metadata. Source hashes were checked again during final review.
This is metadata inspection, distinct from the already completed native MuData
round-trip tests. The large objects and complete decoded metadata stay outside
this branch.

## Original-generation table

Values below are summarized for reading; the TSV retains full paths, UUIDs,
flags, field locations and classifications. A D path identifies recorded text,
not the contents of the file that occupied that path.

| Required field | Smart-seq3 / E-MTAB-8735 | IMR90 / GSE278958 | Validated HAP1 batch10 / GSE278952 |
| --- | --- | --- | --- |
| Study, input sample/run accessions | **D:** study/diySpike label and pooled filenames. **U:** sample/library/run accessions and download checksums. Repository `emtab8735_runs.tsv` contains explicitly marked planning placeholders. | **D:** GSM8558852–GSM8558874 in stored mapping. Curated metadata links GSM8558852→SRX26321174→SRR30918126 and GSM8558853→SRX26321183→SRR30918117. **U:** other 21 run links and original FASTQ execution binding. | **D:** batch10 execution reported. Three candidate/pilot mappings exist (SRR30911454, SRR30911453, SRR30911559). **U:** membership of these in batch10 and the complete ten-run list. |
| Accession→library→cell map | **D:** DIYHEK_1–DIYHEK_192 stored; annotation path/column names survive. **U:** index-pair sequences and accession/library mapping. | **D:** GSM, sample title, original RNA/DNA IDs and canonical cell IDs survive; exact 23 rows in cell TSV. **U:** original map/SOFT file hashes and complete SRR linkage. | **U:** original batch10 mapping. Pilot metadata is only a lead. |
| Inclusion/exclusion ledger | **D:** 192 selected and aligned in stored summary; final 192-row axis. **U:** original all-manifest bytes, full input annotation, exclusions/reasons. | **D:** planned/completed aggregate 23 and final 23-cell axis. **U:** manifest_all.tsv bytes, wider cohort and exclusions/reasons. | **D:** documented scale ten. **U:** exact included/excluded cells and reasons. |
| Pooled/per-cell layout and barcode handling | **D:** pooled R1/R2/I1/I2 paths, `cell_id/index1/index2` columns, `max_mismatch=0`, paired-end. **R:** historical index-pair matching and synchronized read-ID checks. **U:** upstream trimming/UMI processing. | **D:** per-cell manifest, single-end, `ramda`, unstranded; demux skipped and paired opt-in false. | **D:** reported paired-end scRR route. **R:** manifest-driven STAR route with explicit paired-RamDA opt-in. **U:** exact strand/barcode settings. |
| Top-N, min reads, read subsampling | **D:** stored `top_n=null`, `max_records=null`, 75,015,128 records, 68,649,627 assigned; 24 threads, parallel=1, chunk_size=1, resume=true. **U:** min-read rule and upstream/earlier-invocation subsampling. | **U:** top-N/min-read/subsampling commands and seeds. Two RNA template totals of 4,000,000 are observations, not proof of a sampling rule. | **U:** all original selection/subsampling settings. |
| Exact command and workflow revision | **D:** preserved option dictionary and workflow UUID `03d6099d-7800-4ac2-be15-6bb93e9debb3`. Stored command_name is `refresh-rna-qc`. **U:** launch/aligner/detector command logs and executing SHA. | **D:** UUID `6c62bf6c-8fa5-49fe-a30e-bc6db7e19232`; `refresh-rna-qc` summary and export record. **U:** exact workflow/aligner/detector commands and executing SHA. | **U:** batch10 launch record and revision. S-Pilot full-stack wrapper does not identify batch10 execution. |
| Historical circyto version | **D:** 0.9.0 in export/refreshed workflow metadata. **U:** exact revision that originally aligned/detected. | **D:** 0.9.0 export/refresh/remap/merge metadata; export timestamp 2026-06-02T04:28:14Z. **U:** exact detection revision. | **U:** original batch10 version. C0 0.10.0 is a later supported implementation snapshot. |
| STAR / BWA / CIRI3 versions and executable hashes | **U:** all three original versions and digests. | **U:** BWA/CIRI3 original versions/digests; STAR is not the documented single-end route. | **U:** all three original versions/digests. |
| Chimeric alignment, rescue and detector settings | **R:** S-May STAR chimeric segment/score/overhang 15, other flags below; BWA rescue `-T 19`; CIRI3 `-S 0 -Ma 1` template. | **R:** S-June direct BWA `-k 15 -T 15`; unsorted SAM; CIRI3 `-S 0 -Ma 0` template. | **R:** S-June paired STAR/BWA rescue/CIRI3 tuple template. **U:** actual flags for batch10. |
| Reference assembly/provider/release/FASTA hash | **D:** `/user/ifrec/liuyuchen/ref_clean/hg38/hg38.fa`, STAR-index path. **U:** authenticated assembly/accession/patch/provider/contig set, index build and FASTA digest. | **D:** same FASTA path in O-IM itself, independently of pilot guide. **U:** precise reference identity/digest. | **R:** S-Pilot hg38 path is a lead. **U:** actual batch10 reference identity/digest. |
| Annotation provider/release/GTF hash | **D:** GENCODE v45-labelled path in stored workflow/RNA summaries. **U:** original GTF bytes/provider verification, detector `-A`, and repair GTF hash. | **D:** v45-labelled path in O-IM itself. **U:** GTF digest, exact detector annotation usage and repair GTF identity. | **R:** S-Pilot v45 path. **U:** actual batch10 annotation. |
| Support units, thresholds, filters, duplicates | **R:** normalizer reads `#junction_reads`; collector support is not UMI counts; cell-total filter and duplicate-coordinate handling below. **U:** executed detector filters, raw read/UMI/PCR duplicate policy. | **R:** same historical support/collector contract. **U:** original threshold/filter/deduplication settings. | **R:** same expected collector; **U:** execution, thresholds and duplicate policy. |
| RNA profiling after detection | **D:** simple-overlap over retained alignments; QNAME grouping; count templates overlapping exactly one GTF gene, exclude ambiguous templates; automatic QC filtering false. | **D:** same stored rule, 23 cell IDs; RNA paths beneath full23_v2. | **U:** batch10 RNA input/counting provenance. |
| Annotation repair, mapping and object assembly | **D:** host_gene_source=`gtf` for 2,379/2,503; 191 shared + one RNA-only membership with 192-row stored circ axis. **U:** exact repair execution and GTF digest. | **D:** 2,429/2,443 GTF annotations; explicit GSM→canonical remap and processed CNV merge, 23 overlap, no missing map IDs. **U:** repair command and source file hashes. | **U:** batch10 final-object chain. The archived full63/RT summary refers to a different directory and cannot substitute. |
| Run→checksum-matched processed object link | **D:** O-SS embeds `emtab8735/work/diySpike_workflow_all192` and workflow UUID. **U:** checksum links for per-cell inputs/results through assembly/repair. | **D:** O-IM embeds `scrr_imr90/work_hg38_full23_v2`, source/remapped MuData paths and UUID. **U:** raw/native→object checksum chain. | **U:** no batch10 checksum-identified object recovered. |

## Historical script reconstruction, not executed Methods

S-May/S-June `circyto/pipeline/align_manifest.py` supply the following STAR
defaults when not overridden:

```text
--outSAMtype SAM --outReadsUnmapped Fastx
--outSJfilterOverhangMin 15 12 12 12
--alignSJoverhangMin 15 --alignSJDBoverhangMin 15
--outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1
--outFilterMismatchNmax 2
--chimSegmentMin 15 --chimScoreMin 15 --chimJunctionOverhangMin 15
```

The paired route feeds CIRI3 a tuple of chimeric junctions, STAR SAM and BWA
rescue SAM (`bwa mem -T 19` over both unmapped mates). Direct single-end BWA
defaults to `-k 15 -T 15`. In `circyto/detectors/ciri3.py`, the direct command
uses `-I`, `-O`, `-F`, `-T`, `-G`, optional `-A`, and defaults to stringency
`-S 0`, mapper `-Ma 1` for STAR or `-Ma 0` for BWA. Command templates and
`CIRCYTO_CIRI3_EXTRA_ARGS` can override behavior. The S-May Smart-seq3 workflow
does not pass a GTF to detector execution; its later RNA/annotation GTF path
must not be treated as proof of detector `-A` use.

The historical normalizer selects BSJ junction-read support, and the collector
constructs a sparse matrix that sums repeated candidate/cell coordinates.
`min_count_per_cell=1` filters **cell-total support**, not each candidate at two
reads. The Smart-seq3 export expands the filtered matrix back onto selected
cells. These code paths offer an explanation for 191 collector columns and a
192-row object; original execution remains unconfirmed. They do not establish
read deduplication or prove the single aggregate `empty` record belongs to
DIYHEK_192. No new thresholds are proposed.

## Conflicts retained

1. **Version:** the later validated implementation is 0.10.0; archived
   export/refresh metadata says 0.9.0. Neither supplies a raw-run source SHA.
2. **Cells/status:** Smart-seq3 stores 192 rows, while original summaries record
   191 shared cells, RNA-only DIYHEK_192, matrix `n_cols=191` but `n_cells=192`,
   and aggregate 191 success/one empty. All numerical results stay unchanged.
   Unknown per-cell completion remains unknown.
3. **Annotation:** the original objects retain a v45-labelled path; the
   manuscript command shapes and curated candidate table mention v38. Those
   commands are placeholders. Neither filename establishes an annotation digest
   or identifies the later host-gene repair input conclusively.
4. **Runtime:** Smart-seq3 stores a refreshed elapsed time of zero and resumed
   workflow timing of 1.376 seconds alongside long alignment/detector stage
   totals. These are not a defensible clean end-to-end runtime benchmark.
5. **Pilot versus final:** the server guide's two IMR90 cells and one/three HAP1
   pilot cells do not prove final23/batch10 settings. O-IM supplies an independent
   full23_v2 path. The unresolved HAP1 archived summary points to
   `work_hg38_full_v3` and 63/56 modality counts; it is excluded from batch10
   evidence and from manuscript results.
6. **Unrelated tool log:** historical `efab3322adfc88a4e12319cd59319c4697046d56:Log.out`
   records STAR 2.7.10b building a temporary test index under `/tmp/tmp.oQo6YX0e2J`.
   It has no manuscript-run link and supplies no original version for these rows.

## Figure regeneration is a separate chain

O0 archive member → exact input SHA-256 → read-only stored modality extraction
→ RNA-only preprocessing/embedding → circRNA overlays and tables → output
checksums. M0:`manuscript/results_v2/provenance.json` records the execution HEAD
`7e9c1cb9a09a4f0bfd114bb66c0c34df4f31fdbc` **with an uncommitted script**, whose
SHA-256 is `4f97474384bb975f3428d49d4effe62a4dad839af0c46c8027672a61147aec32`.
Use that digest and the committed script/requirements, not HEAD alone, to
identify the regeneration code. Scanpy 1.11.5 and seed 17 belong here; they
are not detection metadata. Existing 18/18 comparison and byte-identical
recheck results are preserved. [Reproduction guide](reproduce.md).

The independent ordinary MuData audit has its own input hashes, script digest,
stack, warnings and round-trip comparisons in
[real_object_interoperability.md](real_object_interoperability.md). Neither
chain closes the missing original native-result/command/reference records.

## Bounded search performed and stopping point

Read the existing evidence, submission/readiness/reproduction/figure/table
records; inspected C1 import/native/report acceptance and C2 distribution
validation via `git show`. Queried local/all remote branch heads and PR #19
through GitHub on 2026-10-15; none of the four supplied refs advanced.

Local search covered repository `work/` and `ramda/`, the already existing
`/tmp/circyto-{import,native,readiness,review,distribution}` evidence directories
(excluding dependency caches), the top-level `/home/liuyuchen/data` inventory,
and relevant filenames across all local Git history. Relevant historical
objects, summaries, script revisions and the unrelated root `Log.out` were
inspected. Previous broader searches are documented at
C1:`docs/reviewer_readiness.md` and `docs/external_ciri3_acceptance.md`; this
audit did not repeat a whole-home search. Current local artifacts are processed
objects, official demos, synthetic imports/reports or dry runs, not the missing
original native per-cell outputs and QC bundles.

`/user/ifrec/liuyuchen/circyto_redo` is not mounted here. No server host/session
was supplied and no remote login was attempted. Stop recovery here pending
the scoped records below. This is a statement about accessible evidence, not
that the records never existed.

## Read-only collection list for James — prepared, not sent or executed

Start with the **object-linked** roots:

```text
/user/ifrec/liuyuchen/circyto_redo/emtab8735/work/diySpike_workflow_all192
/user/ifrec/liuyuchen/circyto_redo/scrr_imr90/work_hg38_full23_v2
```

For HAP1, first identify the directory whose original summary/manifest explicitly
matches **batch10**. Inspect only the project root and immediate `work*`
directories under `/user/ifrec/liuyuchen/circyto_redo/scrr_hap1`; do not substitute
`work_hg38`, `work_hg38_fullstack`, or `work_hg38_full_v3` without that link.

Please return only the following records, with original relative paths and
SHA-256 values; indicate missing/deleted files explicitly:

| Need | Small records requested |
| --- | --- |
| Raw input and exact cells | Original selected/all/batch10 manifests; Smartseq3.diySpike.sample_annotation.normalized.tsv; original sample sheet/SDRF or run-selector rows for the included libraries; accession→filename→library→cell mapping; selection/subsampling command and inclusion/exclusion reasons. Existing FASTQ checksum records only, no read data. |
| Execution | Original invocation/job script, workflow_summary.json, align/alignment_manifest.tsv, align/alignment_prepare_summary.json, ciri3/detector_run_summary.json, per-cell alignment/CIRI3 command/version/exit-status records and provenance JSON. Original checkout SHA plus task-relevant patch if dirty; no entire shell history or environment dump. |
| Aligners/detector | Run-bound STAR Log.out/version/parameters, BWA command/version stderr, CIRI3 command/JAR identity/hash; relevant rescue/strand/stringency/mapper/annotation/duplicate settings. If only present-day binaries exist, label their hashes present-day and leave original versions unknown. |
| Reference/annotation | Original download/source/release manifest and FASTA/GTF SHA-256; STAR index genomeParameters.txt and build command. Known leads: ref_clean/hg38/hg38.fa, ref_clean/hg38/gtf/gencode.v45.annotation.gtf, ref_clean/hg38/star_index. A checksum of a currently mutable path alone is not historical identity. |
| Matrix and post-processing chain | Original matrix indices/feature table and checksums; qc/cell_qc.tsv, qc/circ_qc.tsv; rna/rna_import_summary.json; cell map; export/remap/CNV merge/host-repair command records and input/output checksums. Final O-SS/O-IM checksum confirmation. |
| Optional small native acceptance | Inventory only of unchanged per-cell native CIRI3 result files and their sample/cell binding; availability of the corresponding RNA counts/object and metadata. Select files by the rule in [the response](nikaido_response.md#small-real-cell-acceptance-plan), before inspecting circRNA signal. |

This command inventories and hashes **only named small records** in the two
known workflow roots. It prints to stdout, creates no server files, follows no
symlinks and does not read FASTQ/SAM/BAM/native result data. It does not claim
that missing names mean missing provenance; James can identify alternative
task-specific log names from the original job script.

```bash
python3 - <<'PY'
from pathlib import Path
import hashlib
roots = [
    Path('/user/ifrec/liuyuchen/circyto_redo/emtab8735/work/diySpike_workflow_all192'),
    Path('/user/ifrec/liuyuchen/circyto_redo/scrr_imr90/work_hg38_full23_v2'),
]
names = [
    'workflow_summary.json', 'manifests/all_manifest.tsv',
    'manifests/single_end_manifest.tsv', 'manifests/paired_end_manifest.tsv',
    'demux/manifest.tsv', 'demux/demux_summary.json',
    'align/alignment_manifest.tsv', 'align/alignment_prepare_summary.json',
    'ciri3/detector_run_summary.json', 'logs/run_ciri3_summary.json',
    'matrix/cell_index.txt', 'matrix/circ_index.txt',
    'matrix/circ_feature_table.tsv', 'qc/cell_qc.tsv', 'qc/circ_qc.tsv',
    'rna/rna_import_summary.json',
]
for root in roots:
    for name in names:
        path = root / name
        if any(p.is_symlink() for p in [path, *path.parents]):
            print('SYMLINK_NOT_READ', path, sep='\t')
        elif not path.is_file():
            print('MISSING', path, sep='\t')
        elif path.stat().st_size > 10 * 1024 * 1024:
            print('SIZE_REVIEW', path.stat().st_size, path, sep='\t')
        else:
            print(hashlib.sha256(path.read_bytes()).hexdigest(),
                  path.stat().st_size, path, sep='\t')
PY
```

Review the named records for unrelated information before sharing them. No
credentials, unrelated genomic data, full server mirror, raw alignments or
automatic message is requested. This collection request is the minimum next
action; a new raw-data run would be **new validation**, not provenance recovery.
