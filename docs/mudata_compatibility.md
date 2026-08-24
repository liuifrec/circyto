# MuData synchronization compatibility

MuData 0.4 changes `MuData.update()` so it updates global indices/maps without
automatically pulling modality `obs`/`var` columns. MuData 0.3 retains the old
default and emits a `FutureWarning`; upstream provides `pull_on_update=True` to
preserve it. MuData 0.4 also replaces `set_options(...)` with
`settings.override(...)` for scoped settings.

circyto's validated policy is to preserve the historical pulls. Global `obs`
uses the deterministic MuData union/order and global `var` retains the modality
feature union, while modality names, axes, metadata, sparse matrices, and maps
remain unchanged. Missing modality values remain missing. This policy matters
to RNA+circ, candidate-SNV, RT, CNV, partial-overlap, and zero-circ objects.

`circyto.multimodal.sync` now applies `pull_on_update=True` only around each
MuData construct/read/update/write operation and restores the surrounding
setting afterward. It supports both settings APIs and does not install a
persistent global option or filter warnings. Existing multimodal provenance
records `mudata_sync_policy = "pull_on_update=True"`.

Host-gene repair needs an additional explicit step: its updated circ `var`
columns replace the stale `circ:` global columns through `pull_var(mods=["circ"])`
before writing. Without that step, the newer pull implementation can fail when
a categorical global column does not contain a newly repaired host-gene value.

Regression tests cover cell/feature identity and order, metadata values and
dtypes, missingness, modality inventory, maps, sparse shapes/values, and H5MU
round trips for RNA+circ, RNA+circ+candidate-SNV, RNA+circ+RT, RNA+circ+CNV,
zero-circ, and partially overlapping cell indices. Focused CLI tests cover
export, inspection/QC, SComatic, scRR RT/CNV, host-gene repair, and manuscript
readers.

The release environment is fully validated with MuData 0.3.10 and AnnData
0.11.4. A targeted Python 3.12 check with MuData 0.4.1 verified the explicit
construction and host-gene refresh paths, but the associated AnnData 0.13.2 /
Pandas 3 stack changes string round-trip dtypes. Broad MuData 0.4 support is
therefore not claimed until that dependency-stack migration receives its own
semantic validation.

The full release suite changed from 337 passed, 8 skipped, 131 warnings (126
MuData behavior-change warnings) to 345 passed, 8 skipped, 5 warnings (zero
MuData behavior-change warnings). The remaining warnings are four duplicate
feature-name warnings from intentionally overlapping manuscript fixtures and
one AnnData string-index conversion warning.

Upstream references: [MuData changelog](https://github.com/scverse/mudata/blob/main/CHANGELOG.md)
and [scoped settings override](https://mudata.readthedocs.io/latest/generated/mudata.settings.override.html).
