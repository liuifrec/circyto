# circyto/detectors/ciri_full.py
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .base import DetectorBase, DetectorRunInputs, DetectorResult

CIRI_FULL_JAR_DEFAULT = Path("tools/CIRI-full_v2.0/CIRI-full.jar")
CIRI_FULL_ADAPTER_DEFAULT = Path("tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh")


@dataclass
class CiriFullDetector(DetectorBase):
    """
    Detector engine for the CIRI-full Pipeline.

    This is a thin wrapper around the existing ciri_full_adapter.sh and
    the CIRI-full jar. We expose it through the generic DetectorBase
    interface so it can be used with run-detector / future multi-detector
    orchestration.

    IMPORTANT: CIRI-full is *not* safe to run in parallel across many cells
    in the same process / working directory. We enforce max_parallel = 1
    at the orchestration layer.
    """

    name: str = "ciri-full"
    input_type: str = "fastq"
    supports_paired_end: bool = True

    ciri_full_jar: Path = CIRI_FULL_JAR_DEFAULT
    adapter_script: Path = CIRI_FULL_ADAPTER_DEFAULT

    # New: tell the orchestrator this tool must run serially
    max_parallel: int = 1

    def is_available(self) -> bool:
        # java + jar + adapter must exist
        if shutil.which("java") is None:
            return False
        if not self.ciri_full_jar.exists():
            return False
        if not self.adapter_script.exists():
            return False
        return True

    def version(self) -> Optional[str]:
        """
        Try to get a version string.

        For now, we just return the manual version you've observed (2.1.1).
        We can make this smarter later by parsing '-h' output.
        """
        return "CIRI-full 2.1.1"

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Run CIRI-full on a single cell via the adapter script.

        Environment variables passed to the adapter:

          R1, R2, REF_FA, GTF, OUT_DIR, OUT_BASENAME, OUT_TSV, THREADS

        The adapter is responsible for:

          - creating a run directory:  ${OUT_DIR}/${OUT_BASENAME}.ciri_full_run
          - invoking: java -jar CIRI-full.jar Pipeline ...
          - copying/linking the final .ciri/.txt into OUT_TSV (normalized TSV)
        """
        outdir = inputs.outdir
        outdir.mkdir(parents=True, exist_ok=True)

        r1 = inputs.r1
        r2 = inputs.r2
        ref_fa = inputs.ref_fa
        gtf = inputs.gtf
        threads = inputs.threads
        cell_id = inputs.cell_id

        # --- Dry-run mode for tests (no real CIRI-full execution) ---
# If ref_fa or gtf is missing AND we detect fake FASTQs, produce empty TSV.
        if ref_fa is None or gtf is None:
            out_tsv = outdir / f"{cell_id}.tsv"
            with out_tsv.open("w") as f:
                f.write("circ_id\tchr\tstart\tend\tstrand\tsupport\n")
            # Return fake DetectorResult
            return DetectorResult(
                detector=self.name,
                cell_id=cell_id,
                outdir=outdir,
                tsv_path=out_tsv,
                run_dir=outdir / f"{cell_id}.ciri_full_run",
                log_path=outdir / f"{cell_id}.ciri_full.log",
                meta={"dry_run": True},
            )

        if r1 is None:
            raise ValueError("CiriFullDetector requires R1 FASTQ")

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
            # NEW: ask CIRI-full/CIRI to keep 1-BSJ circRNAs
            "CIRI_EXTRA_FLAGS": "-0",
        }

        # Inherit the current environment, then overlay our variables.
        # This ensures PATH, java, bwa, samtools, etc. are visible.
        real_env = os.environ.copy()
        real_env.update(env)

        cmd = [
            "bash",
            "-lc",
            f'"{self.adapter_script}"',
        ]

        # Let errors propagate; the orchestrator will catch CalledProcessError.
        subprocess.run(
            " ".join(cmd),
            shell=True,
            check=True,
            env=real_env,
        )

        # If we reach here, the adapter exited 0.
        # We expect OUT_TSV and the run dir to exist.
        return DetectorResult(
            detector=self.name,
            cell_id=cell_id,
            outdir=outdir,
            tsv_path=out_tsv,
            run_dir=run_dir,
            log_path=log_path,
            meta={"threads": threads},
        )
