# tests/test_regression_chr21_manifest2_fuzzy.py

from __future__ import annotations

from pathlib import Path

from circyto.analysis.compare_detectors import load_keys, fuzzy_hits, fuzzy_jaccard


def _write_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def test_chr21_manifest2_cirifull_is_fuzzy_recovered_by_fc3(tmp_path: Path):
    """
    Regression test: CIRI-full calls must be fuzzy-recovered by find-circ3 within ±5bp.

    CI-safe behavior:
      - If real outputs exist under work/multidetector_chr21_manifest2/matrices, validate them.
      - Otherwise, run a synthetic mini-regression to ensure fuzzy matching + strand wildcard works.
    """
    base = Path("work/multidetector_chr21_manifest2/matrices")
    a = base / "ciri-full.circ.txt"
    b = base / "find-circ3.circ.txt"

    if a.exists() and b.exists():
        A = load_keys(a)
        B = load_keys(b)
        assert len(A) == 2, f"Expected 2 CIRI-full calls, got {len(A)}"
        assert len(B) >= 2, "find-circ3 produced empty/too-small callset unexpectedly"

        ha = fuzzy_hits(A, B, window=5)
        hb = fuzzy_hits(B, A, window=5)
        fj = fuzzy_jaccard(A, B, window=5)

        assert ha == len(A), f"Expected all CIRI calls fuzzy-hit fc3 (ha={ha}, |A|={len(A)})"
        assert hb >= len(A), f"Expected at least |A| fc3 calls fuzzy-hit back to CIRI (hb={hb})"
        assert fj > 0.0, f"Expected non-zero fuzzy jaccard, got {fj}"
        return

    # ---------------------------
    # CI-safe synthetic fallback
    # ---------------------------
    # CIRI-full often has unknown strand (?), which should act like a wildcard.
    # We craft fc3 calls that match within ±5 bp.
    a2 = tmp_path / "ciri-full.circ.txt"
    b2 = tmp_path / "find-circ3.circ.txt"

    ciri = [
        "chr21:33523436:33527305:?",  # unknown strand
        "chr21:39224408:39238573:?",  # unknown strand
    ]
    fc3 = [
        "chr21:33523435:33527305:+",  # start -1, strand +
        "chr21:39224407:39238573:-",  # start -1, strand -
        "chr21:10105033:10119240:-",  # extra noise call
    ]

    _write_lines(a2, ciri)
    _write_lines(b2, fc3)

    A = load_keys(a2)
    B = load_keys(b2)

    assert len(A) == 2
    assert len(B) >= 2

    ha = fuzzy_hits(A, B, window=5)
    hb = fuzzy_hits(B, A, window=5)
    fj = fuzzy_jaccard(A, B, window=5)

    assert ha == len(A), f"SYN: expected all CIRI calls fuzzy-hit fc3 (ha={ha}, |A|={len(A)})"
    assert hb >= len(A), f"SYN: expected >=|A| fc3 calls fuzzy-hit back to CIRI (hb={hb})"
    assert fj > 0.0, f"SYN: expected non-zero fuzzy jaccard, got {fj}"
