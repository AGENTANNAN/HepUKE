"""Replay collector for the look-ahead KL-gain predictor (T7).

Given a corpus of queries and a fully-configured :class:`RagAgent`, this
module runs each query with the confidence controller **disabled** and a
fixed round budget ``T_max``, snapshotting per-round features into JSONL
so downstream training (T8), BOCPD PIT calibration (T9), and fusion
calibration (T10) can consume them offline.

Design notes (see entropy/section_method.tex §4.3):

* The collector never enables the controller — that would shorten
  trajectories and skew the empirical distribution of ``t`` (paper's
  train/test mismatch caveat).
* Each record contains the *round's* posterior ``p_t``, the union of
  ``doc_id``s newly returned this round, the running features ``φ(H_t)``
  and — after the trajectory ends — the realised KL gain
  ``g_{t+1} = KL(p_{t+1} ‖ p_t)`` as the self-supervised label.
* This module builds on the ``round_observer`` hook in ``loop.py``; it
  does not fork the loop.
"""

from __future__ import annotations

import json
import math
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from ...config import ConfidenceConfig
from .signals import compute_entropy, compute_novelty


# ---------------------------------------------------------------------------
# per-round snapshot
# ---------------------------------------------------------------------------

@dataclass
class RoundSnapshot:
    """One round's raw observations. Contains only the fields the trainer
    actually needs; everything else is available on the full session."""

    t: int
    posterior: list[float]
    doc_ids: list[str]  # new IDs added THIS round
    all_dry: bool
    tokens_spent: int
    search_calls: int  # cumulative `search` tool calls up to and incl. this round


# ---------------------------------------------------------------------------
# features φ(H_t)
# ---------------------------------------------------------------------------

def extract_features(
    round_index: int,
    posterior: Sequence[float],
    current_doc_ids: Sequence[str],
    history_doc_ids: Sequence[str],
    past_gains: Sequence[float],
    expected_docs_per_round: int = 5,
) -> dict[str, float]:
    """Lightweight per-round features. Corresponds to φ(H_t) in eq.~(8).

    Groups:
      * posterior descriptors  — entropy, top1_prob, posterior_var
      * per-round retrieval    — n_unique_docs, novelty
      * historical dynamics    — gain_mean / gain_var / gain_slope over
                                 the last ``L`` realised gains
      * round index            — normalised ``t / T_max`` handled downstream

    Note: query embeddings are the fourth φ group in the paper but they're
    per-query, not per-round; the trainer joins them from a separate cache.
    """
    n = len(posterior)
    if n:
        top1 = max(posterior)
        mean = sum(posterior) / n
        var = sum((p - mean) ** 2 for p in posterior) / n
    else:
        top1 = 0.0
        var = 0.0
    entropy = compute_entropy(posterior)
    novelty = compute_novelty(current_doc_ids, history_doc_ids)
    n_unique = len(set(current_doc_ids))

    # Rolling stats over past L gains (L=3 by default — smallest window
    # that still captures slope; larger windows added later once we see
    # real HotpotQA distributions).
    L = 3
    tail = list(past_gains)[-L:]
    if tail:
        g_mean = sum(tail) / len(tail)
        g_var = sum((g - g_mean) ** 2 for g in tail) / len(tail)
        if len(tail) >= 2:
            # Least-squares slope with x = [0..len-1].
            xs = list(range(len(tail)))
            xm = sum(xs) / len(xs)
            num = sum((xs[i] - xm) * (tail[i] - g_mean) for i in range(len(tail)))
            den = sum((x - xm) ** 2 for x in xs)
            g_slope = num / den if den else 0.0
        else:
            g_slope = 0.0
    else:
        g_mean = 0.0
        g_var = 0.0
        g_slope = 0.0

    denom = max(1, int(expected_docs_per_round))
    exhaustion = 1.0 - min(1.0, n_unique / denom)

    return {
        "t": float(round_index),
        "entropy": entropy,
        "top1_prob": float(top1),
        "posterior_var": float(var),
        "n_unique_docs": float(n_unique),
        "novelty": novelty,
        "gain_mean": g_mean,
        "gain_var": g_var,
        "gain_slope": g_slope,
        "exhaustion": float(exhaustion),
    }


# ---------------------------------------------------------------------------
# KL(p_{t+1} ‖ p_t)
# ---------------------------------------------------------------------------

def kl_divergence(p: Sequence[float], q: Sequence[float]) -> float:
    """KL(p ‖ q) with matched top-B support.

    When ``p`` and ``q`` differ in length, only the shorter prefix is
    used — top-B logprobs at successive rounds are over the same K best
    continuations at each round, so a length mismatch means the
    continuation set shifted; a coarse prefix-KL is the honest fallback
    without token-level alignment (paper §4.2 uses length-normalised
    softmax, so this stays consistent).

    Returns 0.0 when either vector is empty. Zero-mass ``p`` entries are
    skipped; zero-mass ``q`` entries under nonzero ``p`` yield ``+inf``
    (real behaviour of KL) — the trainer clips this to a large finite
    value before Huber loss.

    Kept for API back-compat and diagnostics. New collector runs use
    :func:`js_divergence` (bounded in [0, log 2] and symmetric — see
    memory/plan_confidence_controller.md for why KL exploded on the
    one-hot posteriors HotpotQA produces).
    """
    if not p or not q:
        return 0.0
    n = min(len(p), len(q))
    total = 0.0
    for i in range(n):
        pi = p[i]
        qi = q[i]
        if pi <= 0.0:
            continue
        if qi <= 0.0:
            return math.inf
        total += pi * math.log(pi / qi)
    return total


def js_divergence(p: Sequence[float], q: Sequence[float]) -> float:
    """Jensen–Shannon divergence ``JSD(p, q) = ½·KL(p‖m) + ½·KL(q‖m)``, ``m = (p+q)/2``.

    Bounded in ``[0, log 2] ≈ [0, 0.693]``, symmetric, and finite even when
    one distribution is one-hot — the properties KL lacks that made the
    initial HotpotQA collector emit ``g_next`` values up to 8.18 whenever
    the parent LLM's top-B collapsed to a single token (posterior
    ``[1.0, 0.0, ..., 0.0]``). See :mod:`docs` for the derivation.

    Zero-mass entries under nonzero mass in the *other* distribution are
    handled correctly by construction: ``m`` averages them so no term
    divides by zero.

    Returns 0.0 when either vector is empty. When the vectors differ in
    length the shorter prefix is used (same rationale as
    :func:`kl_divergence`). Mass renormalisation is skipped — top-B
    posteriors from the parent LLM already sum to 1 (within numerical
    slack); slack < 1e-6 is absorbed by the arithmetic mean.
    """
    if not p or not q:
        return 0.0
    n = min(len(p), len(q))

    total = 0.0
    for i in range(n):
        pi = p[i]
        qi = q[i]
        mi = 0.5 * (pi + qi)
        if mi <= 0.0:
            # Both p and q vanish at this index → contributes 0 to JSD.
            continue
        if pi > 0.0:
            total += 0.5 * pi * math.log(pi / mi)
        if qi > 0.0:
            total += 0.5 * qi * math.log(qi / mi)
    return total


# ---------------------------------------------------------------------------
# collector
# ---------------------------------------------------------------------------

@dataclass
class TrajectoryCollector:
    """Callable observer that snapshots one query's trajectory.

    Usage::

        col = TrajectoryCollector()
        agent = RagAgent(..., round_observer=col)
        result = agent.run(query, max_steps=T_max)
        col.finalize()  # computes g_{t+1} labels
        col.records     # list[dict], one per round

    Then serialise ``col.records`` to a JSONL file via
    :func:`dump_records`.
    """

    query_id: str = ""
    query: str = ""
    T_max: int = 8
    _snapshots: list[RoundSnapshot] = field(default_factory=list)
    _past_gains: list[float] = field(default_factory=list)
    records: list[dict] = field(default_factory=list)

    def __call__(self, step_index: int, session: Any) -> None:
        """`round_observer` entry point. `step_index` is 0-based; the
        paper's t is 1-based, so we store t = step_index + 1."""
        posterior = (
            list(session.posteriors[-1]) if getattr(session, "posteriors", None) else []
        )
        # Doc IDs newly seen at THIS step (observed_hits rows are stamped
        # with `step` == step_index at the time of dispatch).
        doc_ids: list[str] = []
        for hit in getattr(session, "observed_hits", []) or []:
            if not isinstance(hit, dict):
                continue
            if hit.get("step") != step_index:
                continue
            did = hit.get("doc_id")
            if did:
                doc_ids.append(str(did))

        # Cumulative `search` tool-call count from the trace (parallel
        # dispatches count separately, matching the eval-side field). This
        # is the cost signal fusion consumes; kept alongside tokens_spent so
        # comparisons stay possible.
        search_calls = 0
        for entry in getattr(session, "trace", []) or []:
            if isinstance(entry, dict) and entry.get("tool") == "search":
                search_calls += 1

        snap = RoundSnapshot(
            t=step_index + 1,
            posterior=posterior,
            doc_ids=doc_ids,
            all_dry=False,  # collector-side heuristic left to finalize()
            tokens_spent=int(getattr(session, "tokens_spent", 0)),
            search_calls=int(search_calls),
        )
        self._snapshots.append(snap)

    def finalize(self) -> None:
        """Compute realised KL gains and pack ``records`` for dumping."""
        history: list[str] = []
        self.records = []
        for i, snap in enumerate(self._snapshots):
            # Realised divergence of THIS round vs. the next round's posterior.
            # Using JSD (bounded, symmetric, finite on one-hot) instead of the
            # old KL — see :func:`js_divergence` for the failure mode KL had.
            # The final round has no successor → label = None; the trainer
            # skips those rows or treats them as "unlabelled".
            if i + 1 < len(self._snapshots):
                g_label: float | None = js_divergence(
                    self._snapshots[i + 1].posterior, snap.posterior
                )
            else:
                g_label = None

            feats = extract_features(
                round_index=snap.t,
                posterior=snap.posterior,
                current_doc_ids=snap.doc_ids,
                history_doc_ids=history,
                past_gains=self._past_gains,
            )
            row = {
                "query_id": self.query_id,
                "query": self.query,
                "t": snap.t,
                "T_max": self.T_max,
                "posterior": snap.posterior,
                "doc_ids": snap.doc_ids,
                "history_doc_ids_count": len(history),
                "features": feats,
                "g_next": g_label,
                "tokens_spent": snap.tokens_spent,
                "search_calls": snap.search_calls,
            }
            self.records.append(row)

            if g_label is not None and math.isfinite(g_label):
                self._past_gains.append(g_label)
            history.extend(snap.doc_ids)


def dump_records(records: Iterable[Mapping[str, Any]], out_path: str | Path) -> int:
    """Append records to a JSONL file. Returns the number of rows written.

    Parent directories are created on demand. Each call is atomic per
    line (small enough that partial writes never cross a newline)."""
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with p.open("a", encoding="utf-8") as f:
        for row in records:
            f.write(json.dumps(row, ensure_ascii=False))
            f.write("\n")
            n += 1
    return n


# ---------------------------------------------------------------------------
# batch driver — thin wrapper so a script can call one function
# ---------------------------------------------------------------------------

def collect_trajectories(
    agent_factory,
    queries: Iterable[Mapping[str, str]],
    *,
    T_max: int = 8,
    out_path: str | Path,
    on_query_end=None,
) -> int:
    """Drive ``queries`` through ``agent_factory()``-produced agents.

    ``agent_factory`` is a zero-arg callable returning a **fresh**
    ``RagAgent`` with the ``round_observer`` slot free — the collector
    installs its own. It's called once per query so pooled state
    (session, tokens_spent, observed_hits) resets cleanly.

    Ensures the confidence controller is OFF on each agent (otherwise
    trajectories get shortened and the t-distribution is skewed —
    §4.3 train/test mismatch caveat).

    Returns total rows written.
    """
    total = 0
    for spec in queries:
        qid = str(spec.get("id") or spec.get("query_id") or "")
        qtext = str(spec.get("query") or spec.get("question") or "")
        if not qtext:
            continue

        agent = agent_factory()
        if getattr(agent, "_confidence_controller", None) is not None:
            raise RuntimeError(
                "collector requires the confidence controller OFF; "
                "pass confidence=None or ConfidenceConfig(enabled=False)."
            )
        col = TrajectoryCollector(query_id=qid, query=qtext, T_max=T_max)
        agent.round_observer = col

        # Mint a stable session_id per query so retrievals / tool_results_full
        # / per_turn_summary rows for this trajectory all share one join key.
        sid = uuid.uuid4().hex
        try:
            agent.run(qtext, max_steps=T_max, session_id=sid)
        finally:
            col.finalize()
            total += dump_records(col.records, out_path)
            if on_query_end is not None:
                try:
                    on_query_end(qid, len(col.records))
                except Exception:
                    pass
    return total


__all__ = [
    "RoundSnapshot",
    "TrajectoryCollector",
    "collect_trajectories",
    "dump_records",
    "extract_features",
    "js_divergence",
    "kl_divergence",
]
