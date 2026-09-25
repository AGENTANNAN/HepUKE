"""Pure-Python retrieval + answer metrics. No numpy / pandas dependency."""
from __future__ import annotations

import math
import re
import string
from collections import Counter
from statistics import mean
from typing import Iterable, Sequence


# ---------------------------------------------------------------------------
# retrieval metrics
# ---------------------------------------------------------------------------

def recall_at_k(retrieved: Sequence[str], relevant: Iterable[str], k: int) -> float:
    rel = {r for r in relevant if r}
    if not rel:
        return 0.0
    top = [d for d in retrieved[:k] if d]
    return len(rel & set(top)) / len(rel)


def mrr(retrieved: Sequence[str], relevant: Iterable[str]) -> float:
    rel = {r for r in relevant if r}
    if not rel:
        return 0.0
    for i, d in enumerate(retrieved, 1):
        if d in rel:
            return 1.0 / i
    return 0.0


def ndcg_at_k(retrieved: Sequence[str], relevant: Iterable[str], k: int) -> float:
    """Binary-relevance nDCG@k."""
    rel = {r for r in relevant if r}
    if not rel:
        return 0.0
    dcg = 0.0
    for i, d in enumerate(retrieved[:k], 1):
        if d in rel:
            dcg += 1.0 / math.log2(i + 1)
    ideal_n = min(len(rel), k)
    idcg = sum(1.0 / math.log2(i + 1) for i in range(1, ideal_n + 1))
    return dcg / idcg if idcg > 0 else 0.0


def aggregate_retrieval(
    per_case: Sequence[dict],
    ks: Sequence[int],
) -> dict:
    """Given per-case metric rows, return mean recall/mrr/ndcg."""
    if not per_case:
        return {
            "recall": {k: 0.0 for k in ks},
            "mrr": 0.0,
            "ndcg": {k: 0.0 for k in ks},
        }
    return {
        "recall": {k: mean(r["recall"][k] for r in per_case) for k in ks},
        "mrr": mean(r["mrr"] for r in per_case),
        "ndcg": {k: mean(r["ndcg"][k] for r in per_case) for k in ks},
    }


# ---------------------------------------------------------------------------
# answer metrics (SQuAD-style)
# ---------------------------------------------------------------------------

_ARTICLES = re.compile(r"\b(a|an|the)\b", re.IGNORECASE)
_WS = re.compile(r"\s+")
_PUNCT_TABLE = str.maketrans("", "", string.punctuation)


def normalize_answer(s: str) -> str:
    """Lower / strip punctuation / drop English articles / squash whitespace."""
    if s is None:
        return ""
    s = s.lower()
    s = s.translate(_PUNCT_TABLE)
    s = _ARTICLES.sub(" ", s)
    s = _WS.sub(" ", s).strip()
    return s


def exact_match(pred: str, golds: Iterable[str]) -> float:
    if not pred:
        return 0.0
    p = normalize_answer(pred)
    return float(any(p == normalize_answer(g) for g in golds if g))


def token_f1(pred: str, golds: Iterable[str]) -> float:
    if not pred:
        return 0.0
    p_toks = normalize_answer(pred).split()
    if not p_toks:
        return 0.0
    best = 0.0
    for g in golds:
        if not g:
            continue
        g_toks = normalize_answer(g).split()
        if not g_toks:
            continue
        common = Counter(p_toks) & Counter(g_toks)
        n_same = sum(common.values())
        if n_same == 0:
            continue
        precision = n_same / len(p_toks)
        recall = n_same / len(g_toks)
        f1 = 2 * precision * recall / (precision + recall)
        best = max(best, f1)
    return best


def aggregate_answer(per_case: Sequence[dict]) -> dict:
    if not per_case:
        return {"em": 0.0, "f1": 0.0}
    return {
        "em": mean(r["em"] for r in per_case),
        "f1": mean(r["f1"] for r in per_case),
    }
