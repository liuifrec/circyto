# circyto/detectors/ciri_full.py
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .base import DetectorBase, DetectorRunInputs, DetectorResult


CIRI_FULL_JAR_DEFAULT = Path("tools/CIRI-full_v2.0/CIRI-full.jar")
CIRI_FULL_ADAPTER_DEFAULT = Path("tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh")


@dataclass
class CiriFullDetector:
    """
    Detector engine for the CIRI-full Pipeline.

    This is a thin wrapper around the existing ciri_full_adapter.sh and
    the CIRI-full jar. For now we keep it simple and single-cell oriented;
    later we can add manifest-based orchestration.
    """

    name: str = "ciri-full"
    input_type: str = "fastq"
    supports_paired_end: bool = True

    ciri_full_jar: Path = CIRI_FULL_JAR_DEFAULT
    adapter_script: Path = CIRI_FULL_ADAPTER_DEFAULT

    def is_available(self) -> bool:
        # Check that the adapter script and jar exist, and that java is on PATH
        if shutil.which("java") is None:
            return False
        if not self.ciri_full_jar.exists():
            return False
        if not self.adapter_script.exists():
            return False
        return True

    def version(self) -> Optional[str]:
        """
        Try to get a version string by grepping the manual or jar banner.
        Keep this lightweight; return None on failure.
        """
        # For now, we just hardcode the manual version you observed (2.1.1).
        # This can be made smarter later by parsing the manual or `-h` output.
        return "CIRI-full 2.1.1"

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Run CIRI-full on a single cell using the adapter script.

        We call the adapter via bash, passing environment variables:

          R1, R2, REF_FA, GTF, OUT_DIR, OUT_BASENAME, THREADS

        The adapter is responsible for:
          - creating a run directory
          - invoking java -jar CIRI-full.jar Pipeline ...
          - copying / linking the final .ciri/.txt into OUT_TSV

        Here we set OUT_TSV = <outdir>/<cell_id>.tsv
        """
        outdir = inputs.outdir
        outdir.mkdir(parents=True, exist_ok=True)

        r1 = inputs.r1
        r2 = inputs.r2
        ref_fa = inputs.ref_fa
        gtf = inputs.gtf
        threads = inputs.threads
        cell_id = inputs.cell_id

        if ref_fa is None:
            raise ValueError("CiriFullDetector requires ref_fa")
        if gtf is None:
            raise ValueError("CiriFullDetector requires gtf")

        out_tsv = outdir / f"{cell_id}.tsv"
        run_dir = outdir / f"{cell_id}.ciri_full_run"
        log_path = outdir / f"{cell_id}.ciri_full.log"

        env = {
            "R1": str(r1),
            "R2": str(r2) if r2 is not None else "",
            "REF_FA": str(ref_fa),
            "GTF": str(gtf),
            "OUT_DIR": str(outdir),
            "OUT_BASENAME": cell_id,
            "OUT_TSV": str(out_tsv),
            "THREADS": str(threads),
        }

        # Inherit current environment but override with our vars
        real_env = dict(**env)
        real_env.update(**env)

        cmd = [
            "bash",
            "-lc",
            f'"{self.adapter_script}"',
        ]

        # Run and let errors propagate
        subprocess.run(
            " ".join(cmd),
            shell=True,
            check=True,
            env=real_env,
        )

        # We expect the adapter to have written out_tsv and log_path
        return DetectorResult(
            detector=self.name,
            cell_id=cell_id,
            outdir=outdir,
            tsv_path=out_tsv,
            run_dir=run_dir,
            log_path=log_path,
            meta={"threads": threads},
        )
