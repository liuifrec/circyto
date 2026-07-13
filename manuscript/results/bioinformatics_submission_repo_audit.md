# Bioinformatics Submission Repository Audit

This audit records repository-readiness risks for a Bioinformatics Application Note style submission. It does not change the manuscript version, upload data, rewrite history, or remove existing tracked objects.

## P0 Blockers

- No current P0 blocker is known for code/test publication after the hardened validation outputs are regenerated.
- Before public release, confirm the final archived manuscript-scale `.h5mu` files do not expose private local paths in object metadata, or document that the committed audit redacts those values and distribute sanitized release copies.

## P1 Strongly Recommended

- Publish frozen manuscript-scale `.h5mu` objects through a GitHub Release and/or Zenodo record with SHA-256 checksums matching `manuscript/results/manuscript_object_manifest.tsv`.
- Keep ordinary CI independent of manuscript-scale `.h5mu` files and external detector binaries.
- Keep external-tool and manuscript-object checks gated, documented, and reproducible from explicit commands.
- Maintain the tiny `circyto demo mini` smoke test as the reviewer-facing no-detector check.
- Keep public supplementary and reproducibility outputs free of usernames, workstation names, and private server paths.
- Avoid adding further large binary revisions to normal Git history.

## P2 Nice To Have

- Add a release checklist entry linking the Git tag, GitHub Release assets, Zenodo DOI, SHA-256 manifest, and exact validation command outputs.
- Add optional CI jobs for detector-backed smoke tests in an environment with STAR, BWA, Java, CIRI3, CIRI-full, and find-circ3 installed.
- Add a small generated HTML or notebook index for the manuscript result tables after the final submission package is fixed.
- Consider a future, reviewed Git history/data-layout cleanup only after release artifacts are archived and manuscript references are stable.
