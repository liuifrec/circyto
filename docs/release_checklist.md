# Manuscript Release Checklist

Use this checklist before creating the Bioinformatics Application Note release/Zenodo archive.

- Update `pyproject.toml` and package version metadata.
- Run `python -m pip install -e ".[dev,mudata]"`.
- Run `pytest -q`.
- Run `python -m pip check`.
- Run `circyto --version`, `circyto doctor`, `circyto detectors`, and `circyto detectors --json`.
- Run `circyto demo mini --out demo_out --overwrite` and confirm no large files are produced.
- Confirm external detector tests are either explicitly run in a prepared environment or skipped with clear dependency messages.
- Confirm CIRI3, Smart-seq3, RamDA/scRR, IMR90 CNV, HAP1 RT, and SComatic wording matches `docs/validated_workflows_summary.md`.
- Confirm SComatic outputs are described only as RNA-derived candidate variant signals.
- Confirm no raw FASTQs, reference genomes, completed `.h5mu` objects, private paths, or large intermediate files are included.
- Create a GitHub release from the intended source tag.
- Archive source, lightweight test fixtures, and reproducibility scripts on Zenodo.
- Archive any large public-data-derived result bundles separately if needed, with clear provenance and licensing.
- Update the Zenodo DOI in the manuscript and repository documentation.

Do not publish a GitHub release or Zenodo record from this checklist unless explicitly instructed.
