"""Set-based scorer for BESIII multi-hop QA (class B/C/E).

gold_ans / pred_ans are comma-separated collections of physics items
(decay channels, reconstruction techniques, EvtGen model names, …).
This scorer:

  1. Splits pred/gold on `, ; \n and / 、` and the arrow "→".
  2. Normalises each item: unicode → ASCII physics tokens (π→pi, Λ→Lambda,
     ⁺→+, etc.), lower-case, punctuation stripped, whitespace collapsed.
  3. Reports set-level Precision / Recall / F1 (macro over items) plus
     Exact-Match (pred_set == gold_set).

Aliases live in a small canonical map — extend as new mismatches surface.
"""
from __future__ import annotations

import re
import unicodedata

# ----- unicode → physics ASCII ---------------------------------------- #
_UNICODE_MAP = {
    "π": "pi", "Π": "Pi",
    "μ": "mu", "Μ": "mu",
    "τ": "tau", "Τ": "tau",
    "ρ": "rho", "Ρ": "rho",
    "γ": "gamma", "Γ": "Gamma",
    "η": "eta", "Η": "eta",
    "ω": "omega", "Ω": "Omega",
    "φ": "phi", "Φ": "phi",
    "χ": "chi", "Χ": "chi",
    "ψ": "psi", "Ψ": "psi",
    "Λ": "Lambda", "λ": "lambda",
    "Σ": "Sigma", "σ": "sigma",
    "Ξ": "Xi", "ξ": "xi",
    "Δ": "Delta", "δ": "delta",
    "θ": "theta", "Θ": "theta",
    "ν": "nu", "Ν": "nu",
    "α": "alpha", "β": "beta",
    "⁺": "+", "⁻": "-", "⁰": "0",
    "→": " to ", "⟶": " to ",
    "'": "'", "'": "'",
    "‐": "-", "–": "-", "—": "-",
    " ": " ",   # nbsp
}

# ----- canonical alias map (particle / model names) ------------------- #
_ALIAS = {
    # antiparticle notation
    "anti-lambda0": "anti_lambda0",
    "anti-nu_e": "anti_nu_e",
    "anti-nu_mu": "anti_nu_mu",
    "anti_p": "anti_proton",
    "anti-p": "anti_proton",
    "p_bar": "anti_proton",
    "pbar": "anti_proton",
    # baryons / mesons
    "lambda_0": "lambda0",
    "sigma_+": "sigma+",
    "sigma_-": "sigma-",
    "xi_-": "xi-",
    "xi_0": "xi0",
    # kaons
    "k_s0": "k_s0",
    "k0_s": "k_s0",
    "ks0": "k_s0",
    "kshort": "k_s0",
    "k_l0": "k_l0",
    "klong": "k_l0",
    # tag / recon techniques
    "double-tag": "double_tag",
    "double tag": "double_tag",
    "single-tag": "single_tag",
    "single tag": "single_tag",
    "full-reconstruction": "full_reconstruction",
    "full reconstruction": "full_reconstruction",
    "partial-reconstruction": "partial_reconstruction",
    "partial reconstruction": "partial_reconstruction",
    "missing-mass": "missing_mass",
    "missing mass": "missing_mass",
    "missing-mass reconstruction": "missing_mass",
    "kalman-kinematic-fit": "kalman_kinematic_fit",
    "kalman kinematic fit": "kalman_kinematic_fit",
    "secondary-vertex-fit": "secondary_vertex_fit",
    "secondary vertex fit": "secondary_vertex_fit",
    "isolated-photon-selection": "isolated_photon_selection",
    "isolated photon selection": "isolated_photon_selection",
}

_SPLIT_RE = re.compile(r"[,;\n、]|(?:\s+and\s+)|(?:\s+or\s+)|(?:\s/\s)")


def _apply_unicode_map(s: str) -> str:
    for k, v in _UNICODE_MAP.items():
        if k in s:
            s = s.replace(k, v)
    # nfkd fold anything else (e.g. combining accents) then drop non-ascii
    s = unicodedata.normalize("NFKD", s)
    s = s.encode("ascii", "ignore").decode("ascii")
    return s


def _normalise_item(x: str) -> str:
    s = _apply_unicode_map(x)
    s = s.lower().strip()
    # strip filler words that don't carry meaning in a set item
    s = re.sub(r"\b(the|a|an|via|using|technique|method|methods)\b", " ", s)
    # remove parenthetical asides that carry no physics content
    s = re.sub(r"\([^)]*\)", " ", s)
    # collapse whitespace + strip trailing punctuation
    s = re.sub(r"[.\?!]+$", "", s)
    s = re.sub(r"\s+", " ", s).strip()
    # look up alias (dashes/spaces both variants) then normalise separators
    if s in _ALIAS:
        s = _ALIAS[s]
    else:
        s2 = s.replace(" ", "_").replace("-", "_")
        if s2 in _ALIAS:
            s = _ALIAS[s2]
    # final separator normalisation: spaces → underscores for stability
    # Do the hyphen-collapse BEFORE spaces are replaced, so charge markers
    # (which are followed by a space in raw text) never look like a hyphen
    # between letters.
    s = re.sub(r"(?<=[a-z])-(?=[a-z])", "_", s)
    s = s.replace(" ", "_")
    # strip leading/trailing underscores or plus/minus fluff
    s = s.strip("_")
    return s


def parse_set(text: str) -> set[str]:
    """Best-effort: split a comma/newline/and-separated string into items."""
    if not text:
        return set()
    parts = _SPLIT_RE.split(text)
    out: set[str] = set()
    for p in parts:
        n = _normalise_item(p)
        if n:
            out.add(n)
    return out


def _item_tokens(item: str) -> set[str]:
    """Core content tokens of a normalised item (after `_normalise_item`).

    Split only on underscore (the canonical separator). Charge suffixes
    like "+", "-", "0" stay attached to their particle name so that
    "k+" and "k-" are DIFFERENT tokens (charge-sensitivity is required
    for physics correctness — pi+ ≠ pi- in a decay channel).
    """
    return {t for t in item.split("_") if t}


def _covers(gold_item: str, pred_item: str) -> bool:
    """True iff gold_item's core tokens are contained in pred_item's."""
    g = _item_tokens(gold_item)
    p = _item_tokens(pred_item)
    if not g:
        return False
    return g.issubset(p)


def score(pred: str, gold: str) -> dict:
    """Return {em, precision, recall, f1, pred_set, gold_set}.

    Matching rule (subset-covers): a gold item is credited when any pred
    item's tokens are a super-set of the gold item's tokens (so
    "Partial reconstruction with D_s+-tag methods" credits the gold item
    "partial reconstruction"). A pred item is credited symmetrically when
    it super-sets at least one gold item — so pred items that add novel
    physics (not just modifiers) still count as false positives.
    Exact-match is the strict equality of the token-sets.
    """
    ps = parse_set(pred)
    gs = parse_set(gold)
    if not gs:
        return {"em": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0,
                "pred_set": sorted(ps), "gold_set": sorted(gs)}

    # gold hits
    gold_hit = sum(1 for g in gs if any(_covers(g, p) for p in ps))
    # pred credits: a pred item that covers at least one gold item is "aligned"
    pred_hit = sum(1 for p in ps if any(_covers(g, p) for g in gs))

    prec = pred_hit / len(ps) if ps else 0.0
    rec = gold_hit / len(gs)
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    em = int(ps == gs)
    return {"em": em, "precision": round(prec, 4), "recall": round(rec, 4),
            "f1": round(f1, 4), "pred_set": sorted(ps), "gold_set": sorted(gs)}


# short-ans extraction prompt geared for set-valued questions
SHORT_SET_PROMPT = (
    "Extract the FULL set of items answering the question from the long "
    "answer. Return ONLY a comma-separated list, no prose, no bullets, no "
    "explanation. Use standard ASCII physics notation (pi+ not π⁺; Lambda "
    "not Λ; K_S0 not Kshort). If the question asks about decay channels of "
    "a parent particle X, list ONLY the final states (e.g. 'K+ pi0, K0 pi+') "
    "— DO NOT prefix each item with 'X to' or 'X -> '. If the long answer "
    "does not contain a definitive set, return an empty string."
)


# ------ scalar scorers for single-hop numeric / exact -------------------- #

_NUM_RE = re.compile(r"[<>]=?|==|=|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|\S+")


def _canon_scalar(x: str) -> str:
    """Lowercase, ASCII-fold, strip fluff — for exact-string comparison."""
    if x is None:
        return ""
    s = _apply_unicode_map(str(x))
    s = s.lower().strip()
    s = re.sub(r"\s+", "_", s)
    s = s.strip("_")
    return s


def _numeric_tokens(x: str) -> list[str]:
    """Extract numeric tokens (with optional comparison operator prefix)."""
    if not x:
        return []
    s = _apply_unicode_map(x).lower()
    tokens: list[str] = []
    for m in _NUM_RE.finditer(s):
        tok = m.group(0)
        # keep either comparison ops (>=, <=, ==) or bare numerics
        if re.fullmatch(r"[<>]=?|==|=", tok):
            tokens.append(tok)
        elif re.fullmatch(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", tok):
            tokens.append(tok)
    return tokens


def score_exact(pred: str, gold: str) -> dict:
    """Exact-string scorer with light normalisation.

    EM = canonicalised strings are equal. F1 = EM (there's no partial credit
    for exact scorers). Returns the same dict shape as `score()`.
    """
    p = _canon_scalar(pred)
    g = _canon_scalar(gold)
    em = int(p == g and p != "")
    # partial credit: if the canonical gold is a token inside pred, half
    # credit (helps when LLM wraps answer in extra prose)
    partial = 0.0
    if not em and g:
        if g in p or g.replace("_", " ") in p.replace("_", " "):
            partial = 0.5
    return {"em": em, "precision": float(em) if em else partial,
            "recall": float(em) if em else partial,
            "f1": float(em) if em else partial,
            "pred": p, "gold": g}


def score_numeric(pred: str, gold: str) -> dict:
    """Numeric scorer.

    EM = a numeric token (with matching comparison operator, if any) from
    gold appears in pred. F1 = EM. This is deliberately loose because
    predictions embed the number in prose ("the minimum is 50 MeV/c").
    """
    if not gold:
        return {"em": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0,
                "pred": pred, "gold": gold}
    g_toks = _numeric_tokens(gold)
    p_toks = _numeric_tokens(pred)
    # gold might not parse as pure numeric — fall back to exact match
    if not g_toks or all(re.fullmatch(r"[<>]=?|==|=", t) for t in g_toks):
        return score_exact(pred, gold)
    # canonicalise operator+number by concatenating adjacent op+num pairs
    def _pair(tokens: list[str]) -> set[str]:
        out: set[str] = set()
        i = 0
        while i < len(tokens):
            t = tokens[i]
            if re.fullmatch(r"[<>]=?|==|=", t) and i + 1 < len(tokens):
                out.add(t + tokens[i + 1]); i += 2
            else:
                out.add(t); i += 1
        return out
    gset = _pair(g_toks)
    pset = _pair(p_toks)
    hits = gset & pset
    em = int(bool(hits) and len(hits) == len(gset))
    # F1: fraction of gold's numeric tokens present in pred
    f1 = len(hits) / max(1, len(gset))
    return {"em": em, "precision": float(len(hits) / len(pset)) if pset else 0.0,
            "recall": f1, "f1": f1, "pred": pred, "gold": gold,
            "gold_tokens": sorted(gset), "pred_tokens": sorted(pset)}


def dispatch_score(pred: str, gold: str, scorer: str) -> dict:
    """Route to the correct scorer based on the QA row's `scorer` field."""
    if scorer == "set_f1":
        return score(pred, gold)
    if scorer == "numeric":
        return score_numeric(pred, gold)
    if scorer == "exact":
        return score_exact(pred, gold)
    # unknown → set_f1 as most permissive default
    return score(pred, gold)


# Short-answer extraction prompt for single-hop questions (numeric/exact)
SHORT_SCALAR_PROMPT = (
    "Extract the single short answer to the question from the long answer. "
    "Return ONLY the answer, no explanation, no unit words. "
    "For numeric answers include comparison operators when present in the "
    "long answer (e.g. '>=4', '==2', '50'). For dataset/model names return "
    "the exact identifier as written (e.g. 'PHSP', 'off_data', 'data_4009', "
    "'psip_pipi_Jpsi_gamma_etac_gammagamma'). If the answer cannot be "
    "determined from the long answer, return an empty string."
)
