# Results Checklist

## Repository and Reproducibility

- [ ] `python -m pip install -e ".[dev]"` or documented equivalent succeeds in the manuscript environment.
- [ ] `pytest -q` passes or expected optional skips are documented.
- [ ] `circyto --help`, `circyto workflow --help`, and `circyto repair-host-genes --help` run.
- [ ] SComatic-related text uses RNA-derived candidate variant signal language.
- [ ] Public HAP1/IMR90 analyses are not described as exposure-specific results.

## Processed Objects

- [ ] Smart-seq3 MuData exists at the expected local path.
- [ ] HAP1 RNA+circRNA+RT MuData exists at the expected local path.
- [ ] IMR90 RNA+circRNA+CNV MuData exists at the expected local path.
- [ ] Each object has host-gene repaired circRNA feature metadata.

## Tables

- [ ] Dataset inventory table generated.
- [ ] Host-gene annotation recovery table generated.
- [ ] HAP1 RT/circRNA correlation table generated.
- [ ] HAP1 OLS table generated.
- [ ] IMR90 CNV/circRNA global correlation table generated.
- [ ] IMR90 CNV-high host-gene table generated.
- [ ] Cross-dataset host-gene overlap table generated.
- [ ] Known/novel circRNA summary generated if annotation fields are present.

## Biological Checks

- [ ] HAP1 RT-positive fraction has the expected positive relationship with circRNA burden before covariate adjustment.
- [ ] HAP1 OLS result is interpreted with detected gene number / transcriptional complexity in mind.
- [ ] IMR90 global CNV burden is not overinterpreted as a strong global circRNA-burden driver.
- [ ] IMR90 CNV-high host-gene programs include fibroblast/ECM/stress-associated genes where supported by the data.
- [ ] Local CNV-at-circRNA-locus analysis is reported when coordinates and copy-number state are available.
