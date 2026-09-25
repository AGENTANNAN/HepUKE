"""Three look-back observation streams (entropy / agreement / novelty).

Signals are stateless functions of the current round's artifacts; the
:class:`~hepuke.agent.confidence.bocpd.BocpdState` keeps their history.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence


def compute_entropy(posterior: Sequence[float]) -> float:
    """Shannon entropy of the top-B answer posterior ``p_t`` (natural log).

    ``posterior`` is the length-B probability vector produced by softmax over
    the parent LLM's top-B answer-token logprobs. Zero-mass entries are
    skipped rather than contributing ``0·log0 = NaN``.
    """

    total = 0.0
    h = 0.0
    for p in posterior:
        if p <= 0.0:
            continue
        total += p
        h -= p * math.log(p)
    if total <= 0.0:
        return 0.0
    return h


def compute_top1_certainty(posterior: Sequence[float], top_b: int = 8) -> float:
    """Certainty scalar ∈ [0,1] used as BOCPD's one-dimensional observation.

    Defined as ``1 - H(p_t) / log(top_b)``: normalises entropy against its
    maximum on the top-B simplex so the observation is scale-free across
    corpora. Approaches 1 as ``p_t`` concentrates on one candidate (LLM is
    certain) and 0 as ``p_t`` spreads uniformly (LLM is confused). Uses the
    whole top-B distribution rather than just ``max p_t`` so that flat but
    top-1-stable posteriors (a known blind spot of a max-only observation)
    still register as low-certainty.
    """

    if top_b <= 1:
        return 1.0
    h = compute_entropy(posterior)
    h_max = math.log(top_b)
    val = 1.0 - h / h_max
    # Clip to [0,1] — floating point can push very slightly outside.
    if val < 0.0:
        return 0.0
    if val > 1.0:
        return 1.0
    return val


def compute_agreement(embeddings: Sequence[Sequence[float]]) -> float:
    """Mean pairwise cosine similarity across sub-agent answer embeddings.

    Cosine proxy for the KL-agreement of eq.~(4) in the paper (see
    App.~B-agreement-proxy). Returns 1.0 when fewer than two embeddings are
    given (no disagreement can be observed).
    """

    vecs = [list(v) for v in embeddings]
    if len(vecs) < 2:
        return 1.0

    def _cos(a: Sequence[float], b: Sequence[float]) -> float:
        dot = sum(x * y for x, y in zip(a, b))
        na = math.sqrt(sum(x * x for x in a))
        nb = math.sqrt(sum(y * y for y in b))
        if na == 0.0 or nb == 0.0:
            return 0.0
        return dot / (na * nb)

    total = 0.0
    n = 0
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            total += _cos(vecs[i], vecs[j])
            n += 1
    return total / n if n else 1.0


def compute_novelty(
    current: Iterable[str],
    history: Iterable[str],
) -> float:
    """Jaccard complement between current-round ``doc_id`` set and history.

    ``1.0`` means the round returned entirely fresh documents; ``0.0`` means
    the same set was already retrieved before (echo-chamber failure).
    """

    cur = {d for d in current if d}
    hist = {d for d in history if d}
    if not cur and not hist:
        return 1.0
    union = cur | hist
    inter = cur & hist
    if not union:
        return 1.0
    return 1.0 - len(inter) / len(union)
