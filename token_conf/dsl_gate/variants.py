"""The perturbation variants, in one place.

Two families, one setting each — the two the paper reports:

    k10     retrieval depth. The unperturbed pipeline retrieves top-5; this
            doubles it, so any change is attributable to depth alone.
    para    an LLM paraphrase of the task description, retrieved at the
            unperturbed top-5. Every physical fact is preserved verbatim;
            only the wording moves.

Every stage — generation, judging, summary, tables — imports this module, so
adding a variant means editing one dict. Output files are tagged
``<prefix>_<variant>``, with the prefix defaulting to the split name
(``pert20_k10``, ``pert50_para``, …).
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Description directory for the paraphrase variant; built by
# build_paraphrases.py. `None` means "use the unperturbed descriptions".
PARA_DESC_DIR = ROOT / "token_conf" / "batch" / "desc_para"

# The retrieval depth of the unperturbed pipeline, for reference in tables.
UNPERTURBED_TOP_K = 5

# variant -> (display label, retrieval top-k, description directory)
VARIANTS: dict[str, tuple[str, int, Path | None]] = {
    "k10":  ("k=10",       10, None),
    "para": ("paraphrase",  5, PARA_DESC_DIR),
}

VARIANT_NAMES = list(VARIANTS)


def label(variant: str) -> str:
    return VARIANTS[variant][0]


def top_k(variant: str) -> int:
    return VARIANTS[variant][1]


def desc_dir(variant: str) -> Path | None:
    return VARIANTS[variant][2]


def tag(variant: str, prefix: str = "pert20") -> str:
    """Output tag for one variant, e.g. ``pert20_k10``."""
    if variant not in VARIANTS:
        raise KeyError(f"unknown variant {variant!r} (have: {', '.join(VARIANTS)})")
    return f"{prefix}_{variant}"


def tags(prefix: str = "pert20") -> list[tuple[str, str, str]]:
    """``(variant, display label, output tag)`` for every variant."""
    return [(v, label(v), tag(v, prefix)) for v in VARIANTS]
