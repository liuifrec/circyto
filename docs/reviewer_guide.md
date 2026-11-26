# Reviewer quick guide

**What circyto adds.** Protocol-aware, cell-resolved orchestration around
established circRNA detectors: manifests and demultiplexing, alignment reuse,
QC, host-gene provenance, sparse matrices, AnnData and optional RNA/circ MuData.
See the [documented software comparison](circrna_workflow_comparison.md).
The offline HTML reporter is now present on `main` following the November 20
merges; it should not be treated as evidence from the submitted CMPB study.
For the paper's scope and results, see [the CMPB submission record](cmpb_submission.md).

**Install and try it.** Python 3.10+ and a wheel built from the chosen commit suffice for reporting.
Given a matching wheel and source archive **built from the same source commit** in the current directory:

```bash
python3 -m venv reviewer-env
source reviewer-env/bin/activate
python -m pip install ./circyto-0.10.0-py3-none-any.whl
tar -xzf circyto-0.10.0.tar.gz circyto-0.10.0/examples/qc_report_demo.py
circyto --help
circyto doctor
circyto detectors
python circyto-0.10.0/examples/qc_report_demo.py --outdir reviewer-demo
circyto report --workdir reviewer-demo
```

Use a new demo directory. The source archive supplies the demonstration script;
it is not an installed wheel command. Open `reviewer-demo/qc/report.html`:
**synthetic** 12 cells, 6 candidates, 77 support, 2/12 zero-count cells, 4/6
recorded host genes. Viewing and regeneration require no network. Package
installation may download dependencies. To build these feature artifacts from
a checkout: `python -m pip install build` then `python -m build` (outputs in
`dist/`). Version 0.10.0 remains unchanged; retain the exact feature commit.

**Run real single-end RamDA/scRR.** Use the optional external-tool recipe
from this checkout or source archive:

```bash
tar -xzf circyto-0.10.0.tar.gz circyto-0.10.0/environment.ramda-se.yml
micromamba create -f circyto-0.10.0/environment.ramda-se.yml --strict-channel-priority -y
micromamba run -n circyto-ramda-se python -m pip install ./circyto-0.10.0-py3-none-any.whl
micromamba run -n circyto-ramda-se circyto doctor
```

Use `micromamba activate circyto-ramda-se` in an initialized micromamba shell
for the following commands. Supply your real
per-cell manifest, FASTA and matching GTF as described in
[the workflow guide](qc_reporting.md#plan-and-run-an-established-workflow):

```bash
circyto manifest validate manifest.tsv --strict
bwa index ref/genome.fa
circyto workflow full-length-circrna --manifest manifest.tsv --protocol ramda \
  --genome-fasta ref/genome.fa --gtf ref/genes.gtf \
  --outdir results/ramda --threads 8 --export-h5ad
```

The recipe was tested on Linux x86_64 with BWA 0.7.19, SAMtools 1.22.1 and
OpenJDK 17.0.18: installation, runtime help, doctor and planning only. It does
not include references or establish a new biological validation. BWA indices
must already exist for real execution; dry-run success does not check them.
The single-end direct-SAM path uses BWA/Java; SAMtools is available for BAM
handling. STAR is additionally needed for paired/pooled SMART-Seq3 routes.
Python installation remains independent of these external runtimes. Review
the [existing detector distribution gaps](../THIRD_PARTY_NOTICES.md) before
publishing a new release.

**Find outputs.** Under `results/ramda/`: `matrix/circ_counts.mtx` is candidates
× cells; use `circ_index.txt` and `cell_index.txt` in that directory for order.
`anndata/circ_counts.h5ad` is cells × candidates. `qc/report.html` is the offline
summary, and `qc/metrics.json` records definitions, exact counts and hashes.
Regenerate with `circyto report --workdir results/ramda`. Support is not a UMI
count; failed/unprocessed cells are not biological negatives; recorded
host-gene coverage is not annotation accuracy or circRNA validation.

**Submitted CMPB study:** The [CMPB submission record](cmpb_submission.md)
summarizes the 437-cell evaluation and its analysis boundaries. The historical
commit `c29831c3f6e7a31b708cce493c0be8036176077a` and its
[reproduction instructions](https://github.com/liuifrec/circyto/blob/c29831c3f6e7a31b708cce493c0be8036176077a/manuscript/reproduce.md)
refer to an earlier manuscript-generation snapshot, not a verified exact
reproduction recipe for the November 20 CMPB submission. Do not claim that
running the synthetic report demo regenerates the submitted figures.
A durable processed-data deposit remains to be verified.
