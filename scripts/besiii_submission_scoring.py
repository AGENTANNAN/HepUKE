"""Conservative, deterministic scoring for the frozen BESIII QA benchmark.

Unlike SQuAD normalization this preserves identifiers, signs, particle charges,
and resonance masses. ``answer`` is authoritative: legacy ``answer_num`` was
extracted with a regex that misreads identifiers and digit separators.
Numeric tolerances are taken verbatim from the row and reported. Missing unit
metadata is never guessed from the question or from the prediction.
"""
from __future__ import annotations

import math
import re
import unicodedata
from typing import Any


SCORER_VERSION = "besiii_submission_v1"


def _literal(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "`\"'":
        text = text[1:-1].strip()
    return re.sub(r"\s+", " ", text)


_NUMBER = re.compile(
    r"(?P<op>>=|<=|==|!=|>|<|=)?\s*"
    r"(?P<value>[+-]?(?:(?:\d{1,3}(?:,\d{3})+|\d+(?:_\d+)*)(?:\.\d*)?|\.\d+)"
    r"(?:[eE][+-]?\d+)?)\s*(?P<unit>[A-Za-zμµ°%/][A-Za-z0-9μµ°%/^* -]*)?"
)
_UNITS = {
    "ev": ("energy", 1.0), "kev": ("energy", 1e3),
    "mev": ("energy", 1e6), "gev": ("energy", 1e9),
    "tev": ("energy", 1e12),
    "mm": ("length", 1e-3), "cm": ("length", 1e-2),
    "m": ("length", 1.0), "um": ("length", 1e-6),
    "degree": ("angle", 1.0), "degrees": ("angle", 1.0),
    "deg": ("angle", 1.0), "°": ("angle", 1.0),
    "rad": ("angle", 180.0 / math.pi),
    "%": ("fraction", 0.01),
}


def _number(text: str) -> dict[str, Any] | None:
    text = text.replace("≥", ">=").replace("≤", "<=").replace("−", "-")
    phrases = ((r"at least\s+", ">="), (r"at most\s+", "<="),
               (r"exactly\s+", "=="), (r"greater than\s+", ">"),
               (r"less than\s+", "<"))
    for pattern, op in phrases:
        text = re.sub(r"^" + pattern, op, text, flags=re.I)
    match = _NUMBER.fullmatch(text)
    if not match:
        return None
    value = float(match["value"].replace("_", "").replace(",", ""))
    if not math.isfinite(value):
        return None
    op = match["op"] or "=="
    return {"value": value, "op": "==" if op == "=" else op,
            "unit": (match["unit"] or "").strip()}


def _unit_key(unit: str) -> tuple[str, float]:
    # Preserve unknown unit symbols exactly; only known units have aliases.
    key = unit.replace("μ", "u").replace("µ", "u").lower().replace(" ", "")
    return _UNITS.get(key, ("literal:" + unit, 1.0))


_PHYSICS = str.maketrans({
    "π": "pi", "μ": "mu", "τ": "tau", "ρ": "rho", "γ": "gamma",
    "η": "eta", "ω": "omega", "Ω": "Omega", "φ": "phi", "χ": "chi",
    "ψ": "psi", "Λ": "Lambda", "Σ": "Sigma", "Ξ": "Xi", "Δ": "Delta",
    "ν": "nu", "θ": "theta", "α": "alpha", "β": "beta",
    "⁺": "+", "⁻": "-", "⁰": "0", "−": "-", "→": "->", "⟶": "->",
})
_TECHNIQUES = {
    "double tag": "double_tag",
    "double tag (full reconstruction)": "double_tag",
    "single tag": "single_tag",
    "full reconstruction": "full_reconstruction",
    "partial reconstruction": "partial_reconstruction",
    "missing mass": "missing_mass",
    "missing mass reconstruction": "missing_mass",
    "kalman kinematic fit": "kalman_kinematic_fit",
    "secondary vertex fit": "secondary_vertex_fit",
    "isolated photon selection": "isolated_photon_selection",
}


def _item(text: str) -> str:
    text = _literal(str(text).translate(_PHYSICS))
    technique = re.sub(r"[-_]", " ", text.lower())
    technique = re.sub(r"\s+(?:technique|method)s?$", "", technique)
    technique = re.sub(r"\s+", " ", technique).strip()
    if technique in _TECHNIQUES:
        return _TECHNIQUES[technique]
    # Typographical spaces around charges/parentheses have no meaning;
    # their contents, case, multiplicity, and particle ordering are kept.
    text = re.sub(r"\s*([()+])\s*", r"\1", text)
    return text


def _items(value: Any) -> set[str]:
    parts = value if isinstance(value, (list, tuple, set)) else [value]
    result: set[str] = set()
    for part in parts:
        # A spaced slash in this dataset denotes alternative techniques;
        # J/psi and slash-bearing identifiers remain one item.
        for token in re.split(r"[,;\n、]|\s+/\s+|\s+and\s+", str(part or "")):
            item = _item(token)
            if item:
                result.add(item)
    return result


def score_answer(prediction: str, row: dict[str, Any]) -> dict[str, Any]:
    """Return ``score, em, f1, precision, recall`` and scoring diagnostics.

    ``score`` is tolerance accuracy for numeric, exact match for exact, and
    set F1 for set_f1 rows. ``em`` is tolerance accuracy for numeric rows;
    ``strict_numeric_match`` separately records exact numeric equality.
    """
    pred = _literal(prediction)
    gold = _literal(row.get("answer", ""))
    kind = row.get("scorer", "exact")
    base = {"scorer_version": SCORER_VERSION, "scorer_requested": kind,
            "scorer_used": kind, "score": 0.0, "em": 0, "f1": 0.0,
            "precision": 0.0, "recall": 0.0}
    if kind == "set_f1":
        ps = _items(prediction)
        gs = _items(row.get("answer_set", row.get("answer", "")))
        hits = len(ps & gs)
        precision = hits / len(ps) if ps else 0.0
        recall = hits / len(gs) if gs else 0.0
        f1 = 2 * precision * recall / (precision + recall) if hits else 0.0
        return {**base, "score": f1, "f1": f1, "em": int(bool(gs) and ps == gs),
                "precision": precision, "recall": recall,
                "pred_set": sorted(ps), "gold_set": sorted(gs)}
    if kind == "numeric":
        expected = _number(gold)
        # A digit-bearing identifier is not a numeric literal even if the
        # legacy answer_num field happens to contain an extracted number.
        if expected is None:
            base["scorer_used"] = "exact_fallback"
            base["fallback_reason"] = "gold_is_not_numeric_literal"
        else:
            observed = _number(pred)
            tolerance = row.get("tolerance") or {}
            rtol = float(tolerance.get("rtol", 0.0))
            atol = float(tolerance.get("atol", 0.0))
            if not math.isfinite(rtol + atol) or min(rtol, atol) < 0:
                raise ValueError("Numeric tolerances must be finite and nonnegative")
            base.update(numeric_match=False, strict_numeric_match=False,
                        numeric_gold=expected, numeric_prediction=observed,
                        tolerance={"rtol": rtol, "atol": atol})
            if expected["value"] and atol >= abs(expected["value"]):
                base["tolerance_warning"] = "atol_at_least_gold_magnitude"
            if observed is None:
                return {**base, "reason": "prediction_is_not_numeric_literal"}
            gold_unit = expected["unit"] or str(row.get("unit") or "")
            pred_unit = observed["unit"]
            if bool(gold_unit) != bool(pred_unit):
                return {**base, "reason": "unit_unspecified_on_one_side"}
            value = observed["value"]
            if gold_unit:
                gold_dimension, gold_scale = _unit_key(gold_unit)
                pred_dimension, pred_scale = _unit_key(pred_unit)
                if gold_dimension != pred_dimension:
                    return {**base, "reason": "incompatible_units"}
                value *= pred_scale / gold_scale
            if expected["op"] != observed["op"]:
                return {**base, "reason": "different_comparison_operator"}
            if value * expected["value"] < 0:
                return {**base, "reason": "opposite_numeric_sign"}
            delta = abs(value - expected["value"])
            match = delta <= atol + rtol * abs(expected["value"])
            return {**base, "score": float(match), "em": int(match),
                    "f1": float(match), "precision": float(match), "recall": float(match),
                    "numeric_match": match, "strict_numeric_match": delta == 0,
                    "numeric_delta": delta, "numeric_prediction_in_gold_units": value}
    elif kind != "exact":
        raise ValueError(f"Unknown scorer: {kind!r}")
    match = bool(gold) and pred == gold
    return {**base, "score": float(match), "em": int(match), "f1": float(match),
            "precision": float(match), "recall": float(match)}
