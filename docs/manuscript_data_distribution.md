# Manuscript Data Distribution Recommendation

This repository should keep source code, tests, tiny fixtures, summary tables, checksums, and reproducibility command records in ordinary Git.

Frozen manuscript-scale `.h5mu` objects should be distributed through a versioned GitHub Release and/or Zenodo archive, with SHA-256 checksums committed in Git. This gives reviewers stable, citable binary artifacts without adding more large binary updates to ordinary Git history.

## Options Considered

| Option | Fit | Recommendation |
| --- | --- | --- |
| Ordinary Git | Good for source, docs, tests, tiny fixtures, manifests, checksums, and small TSV/JSON/MD outputs. Poor for repeated large binary object updates. | Keep lightweight reproducibility assets in Git. Do not add further routine manuscript-scale object revisions to Git history. |
| Git LFS | Handles large files but changes clone and hosting behavior, requires LFS quota/account management, and would need a deliberate migration plan. | Do not add Git LFS configuration in this hardening pass. |
| GitHub Release | Good for frozen reviewer artifacts tied to a specific tag, with checksums in the repo. | Recommended for release snapshots of frozen manuscript-scale objects. |
| Zenodo | Good for DOI-backed archival distribution and manuscript citation. | Recommended for final archived manuscript objects, ideally linked from the matching GitHub Release. |
| Tiny repo demo plus archived manuscript objects | Gives reviewers a fast smoke test from Git and separate access to frozen scale objects. | Preferred path: `circyto demo mini` in Git, summary/checksum tables in Git, manuscript-scale `.h5mu` files in Release/Zenodo. |

## Recommended Policy

- Keep code, tests, tiny fixtures, summary tables, validation reports, checksums, and command tables in Git.
- Keep the tiny reviewer demo in the package so `circyto demo mini --out ... --overwrite` requires no external detector tools.
- Archive frozen manuscript-scale `.h5mu` objects in a GitHub Release and Zenodo record.
- Do not destructively migrate the current repository history during this pass.
- Avoid further repeated binary updates in ordinary Git history.
- Preserve existing tracked manuscript objects until a deliberate release/archive migration is reviewed separately.
