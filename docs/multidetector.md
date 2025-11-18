# Multi-detector circRNA runs and comparison

This page documents the **multi-detector pipeline** in `circyto`:

- Running multiple circRNA detectors on the **same manifest**
- Organizing per-detector outputs in a standard layout
- Comparing detector calls via union / intersection and circ×detector matrices

> **Status (v0.6.0)**  
> - A generic detector API and multi-detector orchestrator are implemented.  
> - `CIRI-full` is fully integrated and validated on a chr21 Smart-seq2 subset.  
> - Additional detectors (e.g. `CIRI-long`, `find_circ`) will be plugged into the same framework in v0.6.x–0.7.0.

---

## 1. Concepts

### 1.1. Detector API

Each detector is wrapped in a shared interface:

- `DetectorRunInputs`: standardized inputs for one cell/sample
- `DetectorResult`: standardized outputs pointing to a final circRNA TSV
- `DetectorBase`: the protocol implemented by each detector engine

See:

- `circyto/detectors/base.py`
- `circyto/detectors/ciri_full.py`
- `circyto/detectors/__init__.py`

The current registry is exposed via:

```python
from circyto.detectors import available_detectors

engines = available_detectors()  # e.g. {"ciri-full": CiriFullDetector(...)}
