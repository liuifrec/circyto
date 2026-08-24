"""Compatibility exports for the optional, chemistry-gated CIRI-long adapter.

CIRI-long is deliberately not registered in the default short-read detector
registry. Its input and execution model are defined by the dedicated
``circyto ciri-long`` command family.
"""

from circyto.pipeline.ciri_long_adapter import (
    build_ciri_long_call_argv,
    build_ciri_long_collapse_argv,
    check_ciri_long_readiness,
    run_ciri_long_call_stage,
    run_ciri_long_collapse_stage,
)

__all__ = [
    "build_ciri_long_call_argv",
    "build_ciri_long_collapse_argv",
    "check_ciri_long_readiness",
    "run_ciri_long_call_stage",
    "run_ciri_long_collapse_stage",
]
