# circyto/detectors/ciri2.py
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any

from .base import DetectorBase, DetectorRunInputs, DetectorResult

# NOTE: this matches your current codespace layout:
# tools/CIRI-full_v2.0/bin/ciri2_adapter.sh
CIRI2_ADAPTER_DEFAULT = Path("tools/CIRI-full_v2.0/bin/ciri2_adapter.sh")


@dataclass
class Ciri2Detector(DetectorBase):
    """
    Detector engine for CIRI2 (CIRI2.pl-based circRNA caller).

    This wraps our ciri2_adapter.sh script, which:
      - runs BWA on the provided FASTQ(s)
      - runs CIRI2.pl with low stringency (-0)
      - normalizes the CIRI2 output into a 6-column TSV:

          circ_id, chr, start, end, strand, support

        at the path given by OUT_TSV.
    """

    name: str = "ciri2"
    input_type: str = "fastq"
    supports_paired_end: bool = True

    adapter_script: Path = CIRI2_ADAPTER_DEFAULT

    def is_available(self) -> bool:
        """Check that the adapter script, perl, and bwa exist."""
        if not self.adapter_script.exists():
            return False
        if shutil.which("perl") is None:
            return False
        if shutil.which("bwa") is None:
            return False
        return True

    def version(self) -> Optional[str]:
        """
        Try to read a version string from the bundled manual if present.
        Fall back to a generic label on failure.
        """
        # We know from your layout:
        # tools/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2_manual.txt
        manual = (
            self.adapter_script.parent
            / "CIRI_v2.0.6"
            / "CIRI2_manual.txt"
        )
        try:
            if manual.exists():
                for line in manual.read_text().splitlines():
                    if "CIRI2" in line:
                        return line.strip()
        except Exception:
            pass
        return "CIRI2 (unknown version)"

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Run CIRI2 on a single cell via ciri2_adapter.sh.

        Environment variables passed to the adapter:

          R1, R2, REF_FA, GTF, OUT_TSV, THREADS
        """
        outdir = inputs.outdir
        outdir.mkdir(parents=True, exist_ok=True)

        if inputs.ref_fa is None or inputs.gtf is None:
            raise ValueError("Ciri2Detector requires both ref_fa and gtf")

        cell_id = inputs.cell_id
        threads = inputs.threads

        out_tsv = outdir / f"{cell_id}.tsv"
        run_dir = outdir / f"{cell_id}.ciri2_run"
        log_path = outdir / f"{cell_id}.ciri2.log"

        env: Dict[str, Any] = {
            "R1": str(inputs.r1),
            "R2": str(inputs.r2 or ""),
            "REF_FA": str(inputs.ref_fa),
            "GTF": str(inputs.gtf),
            "OUT_TSV": str(out_tsv),
            "THREADS": str(threads),
        }

        real_env = os.environ.copy()
        real_env.update({k: str(v) for k, v in env.items()})

        # Let the adapter handle all the heavy lifting; we just capture logs.
        cmd = [
            "bash",
            "-lc",
            f'"{self.adapter_script}" > "{log_path}" 2>&1',
        ]

        result = subprocess.run(
            " ".join(cmd),
            shell=True,
            env=real_env,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"[ciri2] adapter exited with code {result.returncode} "
                f"for {cell_id}; see {log_path}"
            )

        return DetectorResult(
            detector=self.name,
            cell_id=cell_id,
            outdir=outdir,
            tsv_path=out_tsv,
            run_dir=run_dir,
            log_path=log_path,
            meta={"threads": threads},
        )
