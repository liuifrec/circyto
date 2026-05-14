from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ProtocolPreset:
    name: str
    full_length_total_rna: bool
    assumes_umi: bool
    assumes_three_prime_counting: bool
    default_aligner: str


_ALIASES = {
    "smartseq3": "smartseq3",
    "smart-seq3": "smartseq3",
    "ramda": "ramda",
    "ramda-seq": "ramda",
    "shin-ramda": "shin-ramda",
    "shin_ramda": "shin-ramda",
    "shinramda": "shin-ramda",
    "shin-ramda-seq": "shin-ramda",
}


PRESETS = {
    "smartseq3": ProtocolPreset(
        name="smartseq3",
        full_length_total_rna=False,
        assumes_umi=True,
        assumes_three_prime_counting=False,
        default_aligner="star",
    ),
    "ramda": ProtocolPreset(
        name="ramda",
        full_length_total_rna=True,
        assumes_umi=False,
        assumes_three_prime_counting=False,
        default_aligner="star",
    ),
    "shin-ramda": ProtocolPreset(
        name="shin-ramda",
        full_length_total_rna=True,
        assumes_umi=False,
        assumes_three_prime_counting=False,
        default_aligner="star",
    ),
}


def normalize_protocol_name(value: str | None) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return ""
    return _ALIASES.get(text, text)


def get_protocol_preset(value: str | None) -> ProtocolPreset:
    normalized = normalize_protocol_name(value)
    if normalized not in PRESETS:
        supported = ", ".join(sorted(PRESETS))
        raise ValueError(f"Unsupported protocol '{value}'. Supported: {supported}")
    return PRESETS[normalized]
