# Reproduce the Application Note evidence

The production baseline is circyto 0.10.0 at
`44697355bcab1c525ca7ef9b130e2ad0094d9e1b`. Use the manuscript branch for the
scripts and committed small outputs. Production code and package metadata are
unchanged from that baseline. circyto orchestrates established detectors; it
is not a new back-splice detection algorithm.

## Environment

Use Python 3.10.20 and an isolated environment with
[`requirements-results.txt`](requirements-results.txt). For example, after
creating and activating that environment:

```bash
python -m pip install -r manuscript/requirements-results.txt
```

Environment setup may require package downloads; the regeneration scripts
themselves are offline. No detector executable or raw sequencing input is
needed. The script imports circyto from this checkout and reports
`0.10.0 (source tree)`, including when circyto is not installed in the environment.
The numerical/plotting stack is pinned; the complete installed-distribution
snapshot, platform, Python build, script digest, time, Git commit, and dirty
worktree state are in [`results_v2/provenance.json`](results_v2/provenance.json).
This is a scientific-stack pin file, not a cross-platform binary lockfile.

## Obtain the two processed objects

| Dataset | Historical archive member | SHA-256 |
| --- | --- | --- |
| E-MTAB-8735 | `load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu` | `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339` |
| GSE278958 | `load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu` | `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8` |

Both objects were recovered from historical Git commit
`c99cddaae2314d3af7d5232150c5b86bb6b52e25`. Sizes are 30,911,142 and
41,036,354 bytes, respectively. They are deliberately absent from the current
tree. In a full clone that contains that historical commit, this read-only
recovery route creates a separate temporary data directory:

```bash
git cat-file -e c99cddaae2314d3af7d5232150c5b86bb6b52e25^{commit}
MANUSCRIPT_INPUTS=$(mktemp -d)
git archive --format=tar --output="$MANUSCRIPT_INPUTS/objects.tar" \
  c99cddaae2314d3af7d5232150c5b86bb6b52e25 \
  load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu \
  load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu
tar -xf "$MANUSCRIPT_INPUTS/objects.tar" -C "$MANUSCRIPT_INPUTS"
```

If the commit is unavailable in a shallow clone, obtain the full history or the
two checksum-matched objects from the authors. A dedicated public processed-data
deposit and final availability statement still require confirmation before
journal submission; no DOI or download endpoint is invented here. Public raw-data
accessions alone are not substitutes for these exact processed bytes.

## One command for results and Figure 1B/C

Run from the repository root after obtaining the objects:

```bash
python scripts/manuscript/regenerate_application_note_results.py \
  --smartseq3 "$MANUSCRIPT_INPUTS/load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu" \
  --imr90 "$MANUSCRIPT_INPUTS/load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu" \
  --output manuscript/results_v2_reproduced
```

Alternatively substitute the two local paths directly. `--outdir` is an alias
for `--output`. The output must be a new directory beneath this checkout's
`manuscript/`; existing directories and symlink escapes are rejected. Use
`manuscript/results_v2` only if it does not already exist. Failed runs can leave
an incomplete directory; only a completed `provenance.json` with matching
output digests identifies a finished run.

Both full input byte streams are SHA-256 checked before object reads or output
creation, and checked again before publishing numerical results. A mismatch
stops the run. Numerical disagreements with historical expectations are
published in `evidence_comparison.tsv` and reported as warnings; the regenerated
numbers remain authoritative. No source H5MU is written. Runtime plotting/JIT
caches stay under the output's ignored `.cache/` directory.

The committed run used the exact relative archive-member paths above from a
temporary materialization directory; provenance omits that machine-specific
root. The recorded command is portable to a repository-root invocation after
substituting the local input paths. Git archive commit + member path + full
checksum define the source identity independently of its local mount point.

## Reader and embedding decisions

The [real-object interoperability audit](real_object_interoperability.md)
establishes that **both checksum-identified originals** open through ordinary
MuData 0.3.10 in the pinned stack and preserve their semantics through native
read/write/read. In a fresh Python process, use:

```python
import mudata

mdata = mudata.read_h5mu("<ORIGINAL_H5MU>")
mdata.write_h5mu("<NEW_ROUNDTRIP_H5MU>")
restored = mudata.read_h5mu("<NEW_ROUNDTRIP_H5MU>")
```

Choose a new output path; preserve the original. No global option override,
warning suppression, or compatibility conversion is needed. The default
native API emits upstream FutureWarnings, documented in the audit.
`circyto.multimodal.sync.read_h5mu` is a separate interface: its explicit-pull
policy fails on Smart-seq3's named modality / unnamed global observation
indices, but reads IMR90. The audit gives the confirmed diagnosis and an
offline verification command; it does not claim a production helper fix.

The figure-regeneration reader still opens HDF5 in mode `r` and uses public
`anndata.io.read_elem` on each `mod/*` group for numerical reconstruction.
That extraction uses only modality matrices, annotations and cell IDs and is
distinct from the successful native MuData tests. No figure method was changed.

Neither source RNA modality has an embedding. Figure 1B/C therefore regenerate
the method in `c99cdda:scripts/manuscript/export_smartseq3_figure1_data.py`:
RNA/circ shared cells in RNA order; RNA `X` normalized to 10,000 counts per
cell; log1p; 2,000 Seurat-flavor highly variable genes; ARPACK PCA with 50 PCs;
15 Euclidean neighbors; UMAP with spectral initialization, min_dist 0.5,
spread 1.0, and seed 17. PCA, neighbors, and UMAP all receive seed 17. The
full parameter record includes effective UMAP epochs and fitted a/b values.
Numerical threads are fixed to one; Scanpy is 1.11.5. The historical
`use_highly_variable=True` argument is retained and emits one Scanpy deprecation
warning. It is unrelated to MuData synchronization warnings.

CircRNA values never enter embedding construction. Metrics are recalculated
from `circ.X`, aligned by cell ID, and overlaid on exactly the same exported
coordinates in both panels. Detection means nonzero support (stored sparse
zeros are excluded); candidate support is unnormalized detector support.
Host annotation fraction means nonempty `host_gene` values, excluding missing
sentinels. It is not biological validation. Floating-point embeddings may vary
across platforms/dependency stacks; the committed coordinates define the final
figure, and the recorded environment enables numerical reproduction.

## Schematics and minimal supplement

```bash
python scripts/manuscript/build_application_note_supplement.py \
  --results manuscript/results_v2_reproduced \
  --output manuscript/supplement_v2_reproduced
```

This verifies the result file digests, produces Figure 1A from the existing
documented architecture scaffold (no prior rendered panel was present), renders
S1/S2, and writes Tables S1-S3 with source references. HAP1 uses only the
committed 10-cell RNA/circ workflow report; optional RT is explicitly synthetic
contract validation. The script does not read a HAP1 object or reinterpret CNV.
The supplied supplement omits optional long-read rows to keep the scope compact.

Use [`figure_legends.md`](figure_legends.md) with the generated SVG originals
and 300-dpi PNG previews. TSVs are the editable source for supplementary tables.

## Focused checks

```bash
python -m compileall -q scripts/manuscript
python -m pytest -q tests/test_manuscript_scripts.py
git diff --check
```

The frozen production result remains 345 passed, 8 skipped, 5 warnings. Focused
manuscript tests are a separate gate, not a replacement for that baseline.
