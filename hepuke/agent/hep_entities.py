"""High-energy-physics (HEP) named-entity heuristics for the RAG router.

Pure regex, zero deps. Splits mentions into two buckets:

  * "generic"  — collaboration / accelerator / ubiquitous parent decays
                 (BESIII, BEPCII, J/psi, psi(2S)). Present in almost every
                 BAM memo — matching them adds no routing signal.
  * "specific" — narrower states, exotic hadrons, paper ids, specific
                 processes with kinematic tags. Two or more of these in
                 one question is a strong signal that the question bundles
                 independent sub-topics and should be decomposed.

Called from `router.classify()` post-LLM: if `count(specific) >= 2` we
override the LLM's decision to `decompose`, one sub-question per specific
entity. The router LLM's classification is respected only when the
heuristic doesn't fire — this is a hard mechanical rule, not a hint.
"""
from __future__ import annotations

import re
from dataclasses import dataclass


# Parenthesised particle states: Zc(3900), X(3872), Y(4260), Zcs(3985),
# also Z_c(...) / Z_c^+(...) with LaTeX-ish decoration and Chinese/Latin
# spacing. `\d{3,4}` catches typical mass tags (~3-4 digits).
_STATE_RE = re.compile(
    r"""
    (?:
        [ZXY]                       # letters commonly used for exotic states
        (?:_?(?:cs|c|s))?           # optional subscript c / s / cs (order matters)
        (?:[+\-\^]?[+\-0])?         # optional charge tag
    )
    \s*\(?\s*(\d{3,5})\s*\)?        # (3900)
    """,
    re.VERBOSE,
)

# BAM paper ids — always specific.
_BAM_RE = re.compile(r"\bBAM-\d{3,5}\b", re.IGNORECASE)

# Charmonium-like process signatures with a hadron-final-state marker.
# Match on multi-particle-final-state hints (e.g. pi+pi-hc, pi+pi-J/psi).
_PROCESS_RE = re.compile(
    r"""
    (?:
        pi\s*\+?\s*pi\s*\-?\s*      # pi+pi-
        (?:h_?c|J/?psi|eta_?c|X|Y|Z) # followed by a heavy state marker
    )
    """,
    re.IGNORECASE | re.VERBOSE,
)

# Generic collaboration / accelerator / ubiquitous parents. Matched to
# EXCLUDE from the specific-entity count (they still show up in text).
_GENERIC_TERMS = {
    "besiii", "bes iii", "bes-iii", "bes3",
    "bepcii", "bepc-ii", "bepc2",
    "belle", "babar", "lhcb", "cleo",
    # J/psi and psi(2S) alone are ubiquitous parent decays. When they
    # appear WITH another specific state they don't add routing signal.
    "j/psi", "j/ψ", "jpsi",
    "psi(2s)", "psi'", "psi(3686)", "ψ(2s)", "ψ'",
}


@dataclass
class EntityScan:
    specific: list[str]
    generic: list[str]

    def n_specific(self) -> int:
        return len(self.specific)

    def to_dict(self) -> dict:
        return {"specific": list(self.specific), "generic": list(self.generic)}


def _norm_state(match: re.Match) -> str:
    """Canonicalise 'Z_c ( 3900 )' → 'Zc(3900)' for stable dedup."""
    whole = match.group(0)
    mass = match.group(1)
    # Strip subscript/underscore/^/space, keep the leading letter + optional c/s
    head = re.sub(r"\s|_|\^|\{|\}", "", whole.split("(")[0] if "(" in whole else whole)
    # Head now looks like `Zc` or `Zcs` or `X`; drop trailing digits (they're
    # the mass we already captured).
    head = re.sub(r"\d.*$", "", head)
    return f"{head}({mass})"


def _norm_bam(m: re.Match) -> str:
    return m.group(0).upper()


def extract(text: str) -> EntityScan:
    """Extract HEP entity mentions from a natural-language question.

    Returns:
        EntityScan with `specific` (routing-relevant) and `generic`
        (background) lists, each de-duplicated preserving first-seen order.
    """
    if not text:
        return EntityScan(specific=[], generic=[])

    lower = text.lower()
    generic: list[str] = []
    seen_g: set[str] = set()
    for term in _GENERIC_TERMS:
        if term in lower and term not in seen_g:
            generic.append(term)
            seen_g.add(term)

    specific: list[str] = []
    seen_s: set[str] = set()

    for m in _STATE_RE.finditer(text):
        canon = _norm_state(m)
        # Skip states that are actually part of a generic parent marker like
        # psi(2S) / psi(3686). We already captured those as generic; the
        # regex is greedy about parenthesised numerals.
        if canon.lower() in seen_g or canon.lower() in _GENERIC_TERMS:
            continue
        # Also skip if the "head" is empty (rare regex artefact).
        if canon.startswith("("):
            continue
        if canon not in seen_s:
            specific.append(canon)
            seen_s.add(canon)

    for m in _BAM_RE.finditer(text):
        canon = _norm_bam(m)
        if canon not in seen_s:
            specific.append(canon)
            seen_s.add(canon)

    for m in _PROCESS_RE.finditer(text):
        canon = re.sub(r"\s+", "", m.group(0)).lower()
        if canon not in seen_s:
            specific.append(canon)
            seen_s.add(canon)

    return EntityScan(specific=specific, generic=generic)


__all__ = ["extract", "EntityScan"]
