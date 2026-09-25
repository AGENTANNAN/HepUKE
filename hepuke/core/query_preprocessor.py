"""Query preprocessing: Unicode-physics normalization + CN→EN term mapping.

Returns `(dense_text, sparse_text)`:
  dense_text  — normalized full question, fed to BGE-M3 multilingual dense
  sparse_text — extracted English keyword string, fed to BGE-M3 sparse (BM25-like)

Ported from rag-besiii `rag/query_preprocessor.py` (2026-05-13). Two BESIII
specifics are removed / pluggable:

  * `_CN_TO_EN` dictionary is loaded from `hepuke/data/cn_en_physics.json`,
    and callers can pass `user_map=` to `Preprocessor` to extend it.
  * BAM-ID normalization (BESIII paper IDs) is now a `DocIdNormalizer`
    Protocol; the default preprocessor does NOT touch IDs. Register your own
    (see `BamIdNormalizer` at the bottom of this file) if you need it.
"""
from __future__ import annotations

import json
import re
import unicodedata
from importlib.resources import files
from typing import Iterable, Protocol, Tuple

__all__ = [
    "DocIdNormalizer",
    "BamIdNormalizer",
    "ArxivIdNormalizer",
    "Preprocessor",
    "normalize_unicode_physics",
    "preprocess",
    "extract_physics_terms_en",
]


# ---------------------------------------------------------------------------
# Unicode → ASCII normalization tables (identical to rag-besiii)
# ---------------------------------------------------------------------------

_GREEK_TO_ASCII: dict[str, str] = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon",
    "ζ": "zeta", "η": "eta", "θ": "theta", "ι": "iota", "κ": "kappa",
    "λ": "lambda", "μ": "mu", "ν": "nu", "ξ": "xi", "ο": "omicron",
    "π": "pi", "ρ": "rho", "σ": "sigma", "ς": "sigma", "τ": "tau",
    "υ": "upsilon", "φ": "phi", "χ": "chi", "ψ": "psi", "ω": "omega",
    "Α": "Alpha", "Β": "Beta", "Γ": "Gamma", "Δ": "Delta", "Ε": "Epsilon",
    "Ζ": "Zeta", "Η": "Eta", "Θ": "Theta", "Ι": "Iota", "Κ": "Kappa",
    "Λ": "Lambda", "Μ": "Mu", "Ν": "Nu", "Ξ": "Xi", "Ο": "Omicron",
    "Π": "Pi", "Ρ": "Rho", "Σ": "Sigma", "Τ": "Tau", "Υ": "Upsilon",
    "Φ": "Phi", "Χ": "Chi", "Ψ": "Psi", "Ω": "Omega",
}

_UNICODE_OPS: dict[str, str] = {
    "→": "->", "↦": "->", "⇒": "->", "⟶": "->", "⟹": "->",
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
    "⁺": "+", "⁻": "-",
    "₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4",
    "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9",
    "×": "x", "·": ".", "≈": "~=", "≠": "!=",
}

_ANTIPARTICLE_BAR_RE = re.compile(r"([A-Za-z])[̄̅]")
_TRANSLATE_TABLE = str.maketrans({**_GREEK_TO_ASCII, **_UNICODE_OPS})


def normalize_unicode_physics(text: str) -> str:
    """Translate Greek letters + physics operators + antiparticle overline to ASCII.

    Λ_c+ → p K⁻ π⁺   →  Lambda_c+ -> p K- pi+
    ψ(2S) → γ χ_cJ  →  psi(2S) -> gamma chi_cJ
    ν̄_e             →  nubar_e
    """
    text = unicodedata.normalize("NFD", text)
    text = text.translate(_TRANSLATE_TABLE)
    text = _ANTIPARTICLE_BAR_RE.sub(r"\1bar", text)
    return text


# ---------------------------------------------------------------------------
# CN → EN physics term map (loaded from packaged JSON)
# ---------------------------------------------------------------------------

def _load_default_cn_en() -> dict[str, str]:
    try:
        raw = files("hepuke.data").joinpath("cn_en_physics.json").read_text(encoding="utf-8")
        return json.loads(raw)
    except Exception:
        return {}


_DEFAULT_CN_TO_EN: dict[str, str] = _load_default_cn_en()


# ---------------------------------------------------------------------------
# DocIdNormalizer — pluggable protocol for domain-specific ID handling
# ---------------------------------------------------------------------------

class DocIdNormalizer(Protocol):
    def normalize(self, text: str) -> str: ...
    def rewrite_query(self, question: str) -> str: ...


class BamIdNormalizer:
    """BESIII BAM-ID normalizer. Matches rag-besiii `normalize_bam_id`."""

    _BAM_RE = re.compile(r"[Bb][Aa][Mm][-_]?(\d+)")

    def normalize(self, text: str) -> str:
        m = self._BAM_RE.search(text)
        return f"BAM-{m.group(1).zfill(5)}" if m else text

    def rewrite_query(self, question: str) -> str:
        return self._BAM_RE.sub(lambda m: f"BAM-{m.group(1).zfill(5)}", question)


class ArxivIdNormalizer:
    """arXiv 4/5-digit ID: `2401.01234` → `arXiv:2401.01234`."""

    _ARXIV_RE = re.compile(r"\b(\d{4}\.\d{4,5})(v\d+)?\b")

    def normalize(self, text: str) -> str:
        m = self._ARXIV_RE.search(text)
        return f"arXiv:{m.group(1)}" if m else text

    def rewrite_query(self, question: str) -> str:
        return self._ARXIV_RE.sub(
            lambda m: f"arXiv:{m.group(1)}", question,
        )


# ---------------------------------------------------------------------------
# Preprocessor
# ---------------------------------------------------------------------------

class Preprocessor:
    """Configurable query preprocessor.

    Parameters
    ----------
    user_map : optional extra CN→EN entries merged over the default map.
    id_normalizer : optional DocIdNormalizer applied to the question before
        Unicode normalization.
    """

    def __init__(
        self,
        *,
        user_map: dict[str, str] | None = None,
        id_normalizer: DocIdNormalizer | None = None,
    ) -> None:
        merged = dict(_DEFAULT_CN_TO_EN)
        if user_map:
            merged.update(user_map)
        self.cn_to_en = merged
        self.id_normalizer = id_normalizer

    def extract_physics_terms_en(self, question: str) -> str:
        terms: list[str] = []
        remaining = question
        for cn, en in self.cn_to_en.items():
            if cn in remaining:
                terms.extend(en.split())
                remaining = remaining.replace(cn, " ")
        en_tokens = re.findall(r"[A-Za-z][A-Za-z0-9/\(\)\+\-]*", remaining)
        terms.extend(en_tokens)

        seen: set[str] = set()
        unique: list[str] = []
        for t in terms:
            if t.lower() not in seen:
                seen.add(t.lower())
                unique.append(t)
        return " ".join(unique)

    def __call__(self, question: str) -> Tuple[str, str]:
        if self.id_normalizer is not None:
            question = self.id_normalizer.rewrite_query(question)
        question = normalize_unicode_physics(question)
        dense_text = question
        sparse_text = self.extract_physics_terms_en(question)
        if not sparse_text.strip():
            sparse_text = question
        return dense_text, sparse_text


# ---------------------------------------------------------------------------
# Module-level convenience API (backwards compatible with rag-besiii)
# ---------------------------------------------------------------------------

_default_preprocessor = Preprocessor()


def preprocess(question: str) -> Tuple[str, str]:
    """Preprocess without any doc-id normalization (generic default)."""
    return _default_preprocessor(question)


def extract_physics_terms_en(question: str) -> str:
    return _default_preprocessor.extract_physics_terms_en(question)
