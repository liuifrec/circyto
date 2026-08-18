# CIRI-long official-demo interoperability validation

## Decision

**PASS** for detector interoperability.

On 2026-07-30, the official CIRI-long `v0.6-alpha` demonstration data ran
successfully with the official CIRI-long `v1.1.0` Python package through
circyto's separate `call`, `collapse`, and `import` stages. All official
BSJ-, isoform-, expression-, isoform-usage-, and read-assignment records were
preserved by normalization. No circyto code change was required.

This validation demonstrates interoperability between circyto and the official
CIRI-long RCRT circRNA workflow using the project's bundled demonstration data.
It validates command execution, provenance capture, output parsing, coordinate
normalization, and preservation of BSJ-, isoform-, expression-, and read-level
information. It does not establish single-cell Nanopore circRNA performance or
biological accuracy on an independent dataset.

## Scope and environment

- Repository commit: `12f5478` (`feature/nanopore-feasibility`)
- circyto version reported in provenance: `0.10.0`
- Host: Dell Precision 5690, Intel Core Ultra 9 185H, 22 online logical CPUs
- OS: Ubuntu 22.04.5 LTS, Linux `6.8.0-136-generic`, x86-64
- Isolated environment: `ciri-long-demo`
- Environment prefix: `$CIRI_LONG_ENV`
- Python: `3.10.20`
- CIRI-long: `1.1.0`
- CIRI-long executable: `$CIRI_LONG_ENV/bin/CIRI-long`
- BWA: `0.7.19-r1273`
- samtools: `1.24`
- CMake: `4.4.1`
- GNU Make: `4.4.1`

CIRI-long came from the official PyPI package `CIRI-long==1.1.0`, corresponding
to the official GitHub `v1.1.0` release. It was not installed into the circyto
environment. Commands invoking CIRI-long used `PYTHONNOUSERSITE=1`, with the
isolated environment first on `PATH`.

### Installation commands and compatibility notes

The environment was created with:

```bash
micromamba create -y -n ciri-long-demo \
  -c conda-forge -c bioconda \
  python=3.10 bwa samtools cmake make pip
```

The first official-package attempt was:

```bash
$CIRI_LONG_ENV/bin/python \
  -m pip install --no-cache-dir CIRI-long==1.1.0
```

It stopped while building `bwapy==0.1.4`: the legacy build imported
`pkg_resources`, which is absent from the setuptools 83 build-isolation
environment. Installation was completed with the following isolated
build-tool compatibility sequence:

```bash
$CIRI_LONG_ENV/bin/python \
  -m pip install --no-cache-dir setuptools==80.9.0

$CIRI_LONG_ENV/bin/python \
  -m pip install --no-cache-dir Cython

$CIRI_LONG_ENV/bin/python \
  -m pip install --no-cache-dir --no-build-isolation bwapy==0.1.4

$CIRI_LONG_ENV/bin/python \
  -m pip install --no-cache-dir CIRI-long==1.1.0

env PYTHONNOUSERSITE=1 \
  $CIRI_LONG_ENV/bin/python \
  -m pip install --no-cache-dir numpy==1.26.4 scipy==1.11.4
```

The last command prevents accidental satisfaction of CIRI-long requirements
from the user's Python site. With `PYTHONNOUSERSITE=1`,
`python -m pip check` reported `No broken requirements found`.

Relevant installed packages included CIRI-long 1.1.0, bwapy 0.1.4, Cython
3.2.9, mappy 2.31, numpy 1.26.4, scipy 1.11.4, pandas 2.3.3, pysam 0.24.0,
pyccs 1.1.0, pyspoa 0.3.2, Biopython 1.87, edlib 1.3.9.post1,
scikit-learn 1.7.2, and python-Levenshtein 0.27.3.

## Demo identity and integrity

The demo source was the official
[`v0.6-alpha` GitHub release](https://github.com/bioinfo-biols/CIRI-long/releases/tag/v0.6-alpha).
The current official usage guide still downloads this asset for its example,
while the installed software is the newer `v1.1.0` release. This version
difference was tested explicitly; the old demo completed successfully under
1.1.0.

| Property | Observed value |
| --- | --- |
| Archive | `CIRI-long_test_data.tar.gz` |
| URL | `https://github.com/bioinfo-biols/CIRI-long/releases/download/v0.6-alpha/CIRI-long_test_data.tar.gz` |
| GitHub asset ID | `26361119` |
| Content length | 36,945,437 bytes |
| Upstream digest | None supplied by the GitHub release API |
| Locally computed SHA-256 | `fb2583c0a1f2461759ad48351fee96ce32384d6703fa5ee88745156f72bcea7e` |
| Extraction directory | `work/ciri_long_official_demo/extracted/test_data/` |
| Redistribution | Not performed |

The extraction destination is ignored by Git through the repository's
`work*/` rule. Approximately 716 GB was available before download.

| Extracted file | Bytes | SHA-256 |
| --- | ---: | --- |
| `test_reads.fa` | 1,460,667 | `71cc00e69b0192e049868aa2f4f27955fdbe4d9e6ec022a2ce691875444759a3` |
| `mm10_chr12.fa` | 122,131,180 | `df14dca22a5007c6995e4320978b386bcae65e60e62b1c9c49ddfad858c3ee03` |
| `mm10_chr12.gtf` | 28,413,519 | `2e1c6e3670fbe38b8c84c60e18ab5f86602873cccc8d3308112a0cec01b8f9fe` |

The official usage documentation identifies the files as an `mm10_chr12`
demo. Validation used the explicit reference ID
`CIRI-long_v0.6-alpha_demo_mm10_chr12` and reference build label
`mm10_official_demo`. The exact upstream FASTA provider and release are not
specified by the release metadata or usage page and remain unresolved; the
validation does not infer a more specific release from the filename.

## Metadata, planning, and readiness commands

```bash
python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_official_demo \
  --metadata-only --timeout-seconds 60

python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_official_demo \
  --dry-run \
  --ciri-long "$CIRI_LONG_ENV/bin/CIRI-long" \
  --bwa "$CIRI_LONG_ENV/bin/bwa" \
  --threads 4

python scripts/run_ciri_long_official_demo.py \
  --outdir work/ciri_long_official_demo \
  --ciri-long "$CIRI_LONG_ENV/bin/CIRI-long" \
  --bwa "$CIRI_LONG_ENV/bin/bwa" \
  --threads 4 --timeout-seconds 120

$CIRI_LONG_ENV/bin/bwa \
  index -a bwtsw \
  work/ciri_long_official_demo/extracted/test_data/mm10_chr12.fa

env PYTHONNOUSERSITE=1 circyto ciri-long doctor \
  --ciri-long "$CIRI_LONG_ENV/bin/CIRI-long" \
  --bwa "$CIRI_LONG_ENV/bin/bwa" \
  --reference work/ciri_long_official_demo/extracted/test_data/mm10_chr12.fa \
  --gtf work/ciri_long_official_demo/extracted/test_data/mm10_chr12.gtf

env PYTHONNOUSERSITE=1 circyto ciri-long plan \
  --manifest work/ciri_long_official_demo/ciri_long_demo_manifest.tsv \
  --reference work/ciri_long_official_demo/extracted/test_data/mm10_chr12.fa \
  --gtf work/ciri_long_official_demo/extracted/test_data/mm10_chr12.gtf \
  --outdir work/ciri_long_official_demo \
  --threads 4 \
  --ciri-long "$CIRI_LONG_ENV/bin/CIRI-long" \
  --bwa "$CIRI_LONG_ENV/bin/bwa"
```

Doctor passed with no errors or warnings after indexing. BWA indexing took
75.74 seconds.

| BWA index asset | Bytes | SHA-256 |
| --- | ---: | --- |
| `mm10_chr12.fa.amb` | 106 | `7a6524ac1e854ddfec6a647036e77df561e1349cd16b8436c4559a35f95dd858` |
| `mm10_chr12.fa.ann` | 44 | `21de4212ed3532e214c288140d91a832b184dcf11ba858ad1bb8e838b8d8593c` |
| `mm10_chr12.fa.bwt` | 120,129,096 | `99126de615f69a3d39c03baf6580af557ac877e40c2bb61c42de174fa1f80ef7` |
| `mm10_chr12.fa.pac` | 30,032,257 | `c79784311bc7f5b2054bbf1cdfde46a5fd3844c42ff038a25809ac57243fc6d5` |
| `mm10_chr12.fa.sa` | 60,064,560 | `3cc7567c2d20941bbade955b80ae847a0d6d0cdca389df61cf694dd1004186ed` |

## CIRI-long call

The manifest chemistry gate accepted one bulk RCRT/circRNA-enriched cDNA
sample. It explicitly rejected interpretation as ordinary poly(A) cDNA,
direct RNA, generic minimap2 alignment input, or a single cell.

The detector argv recorded by circyto, shown with workstation paths normalized, was:

```text
$CIRI_LONG_ENV/bin/CIRI-long
call
-i work/ciri_long_official_demo/extracted/test_data/test_reads.fa
-o work/ciri_long_official_demo/ciri_long/call/official_demo
-r work/ciri_long_official_demo/extracted/test_data/mm10_chr12.fa
-p official_demo
-a work/ciri_long_official_demo/extracted/test_data/mm10_chr12.gtf
-t 4
```

The wrapper was executed with `--execute`, `shell=false`, and the isolated
environment first on `PATH`.

| Call metric | Result |
| --- | ---: |
| Exit status | 0 |
| circyto-recorded runtime | 38.584 s |
| External wall time | 39.24 s |
| Input reads | 1,562 |
| Cyclic consensus reads | 1,354 |
| Candidate sequences / BSJ-bearing reads | 1,204 |
| Splice-signal candidates | 1,189 |
| Low-confidence / partial sequences | 10 |

All documented outputs existed:
`official_demo.cand_circ.fa`, `official_demo.json`, `official_demo.log`,
`official_demo.low_confidence.fa`, `tmp/official_demo.ccs.fa`, and
`tmp/official_demo.raw.fa`.

CIRI-long logged `Index of reference genome not found, generating ...` during
its splice-site preparation and created `tmp/ss.idx`. This occurred despite
doctor having verified all five BWA index files; it refers to CIRI-long's
internal index generation and did not prevent successful execution.

## CIRI-long collapse

The generated `calls.list` contained the `official_demo` sample and its
candidate FASTA. The detector argv, shown with workstation paths normalized, was:

```text
$CIRI_LONG_ENV/bin/CIRI-long
collapse
-i work/ciri_long_official_demo/ciri_long/collapse/calls.list
-o work/ciri_long_official_demo/ciri_long/collapse
-p cohort
-r work/ciri_long_official_demo/extracted/test_data/mm10_chr12.fa
-a work/ciri_long_official_demo/extracted/test_data/mm10_chr12.gtf
-t 4
```

| Collapse metric | Result |
| --- | ---: |
| Exit status | 0 |
| circyto-recorded runtime | 9.984 s |
| External wall time | 10.63 s |
| Initial circular-read clusters | 210 |
| Corrected clusters | 151 |
| Final collapsed BSJs | 149 |
| Total support | 1,133 |
| Isoform structures | 149 |
| Expression matrix | 149 rows × 1 sample |
| Isoform-usage matrix | 149 rows × 1 sample |
| Read assignments | 1,133 |

The official `cohort.info`, `cohort.expression`, `cohort.isoforms`,
`cohort.reads`, and `cohort.log` files all existed.

## circyto import

```bash
env PYTHONNOUSERSITE=1 circyto ciri-long import \
  --collapse-dir work/ciri_long_official_demo/ciri_long/collapse \
  --outdir work/ciri_long_official_demo \
  --prefix cohort
```

| Normalized artifact | Records |
| --- | ---: |
| `circRNA_bsj.tsv` | 149 |
| `circRNA_isoforms.tsv` | 149 |
| `circRNA_expression.tsv` | 149 |
| `circRNA_isoform_usage.tsv` | 149 |
| `circRNA_read_assignments.tsv` | 1,133 |

`import_summary.json` and call/collapse provenance files recorded source and
normalized checksums, exact commands, tool versions, reference identity,
chemistry gating, and scientific warnings.

For real record `chr12:3617606-3773640`, CIRI-long reported 1-based closed
coordinates `3617606–3773640`. circyto retained these in `original_start` and
`original_end`, and emitted 0-based half-open coordinates
`3617605–3773640`, following:

```text
normalized_start = original_start - 1
normalized_end   = original_end
```

## Reconciliation

Numeric comparisons used relative and absolute tolerances of `1e-12`.

| Artifact | Official rows | Normalized rows | Matched identifiers | Unmatched identifiers | Value discrepancies | Status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `cohort.info` → `circRNA_bsj.tsv` | 149 | 149 | 149 | 0 | 0 | PASS |
| `cohort.info` isoforms → `circRNA_isoforms.tsv` | 149 | 149 | 149 | 0 | 0 | PASS |
| `cohort.expression` → `circRNA_expression.tsv` | 149 | 149 | 149 | 0 | 0 | PASS |
| `cohort.isoforms` → `circRNA_isoform_usage.tsv` | 149 | 149 | 149 | 0 | 0 | PASS |
| `cohort.reads` → `circRNA_read_assignments.tsv` | 1,133 | 1,133 | 1,133 | 0 | 0 | PASS |

Additional reconciliation checks passed:

- every official `circ_id` and sample column was retained;
- every official isoform structure was retained;
- expression and isoform-usage values were unchanged;
- read assignments and their complete original rows were retained;
- `.info`, expression, and normalized support totals all equalled 1,133;
- normalized BSJ, isoform, expression, usage, and read-assignment primary keys
  had zero duplicates;
- original official representations remained available in normalized fields;
- no AnnData/MuData export, generic Nanopore output, short-read detector, or
  validated `circ` modality was invoked or modified.

## Limitations and warnings

- The demo archive is attached to the historical `v0.6-alpha` prerelease,
  while this validation used CIRI-long `v1.1.0`. Compatibility succeeded for
  this asset but is not a general guarantee for all historical inputs.
- GitHub supplied no checksum for the demo asset. The recorded archive SHA-256
  is locally computed, not upstream-authenticated.
- The demo release and usage page identify `mm10_chr12`, but do not state the
  exact upstream FASTA provider or release.
- Installing current CIRI-long 1.1.0 required a legacy-build workaround for
  `bwapy==0.1.4`; this did not require a circyto code change.
- The bundled demo is a project demonstration, not an independent biological
  validation dataset. It is bulk and provides no evidence of single-cell
  performance.
- CRR194209 was not downloaded.
