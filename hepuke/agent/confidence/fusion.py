"""Fusion layer: bind look-back Stab_t and look-ahead ĝ into a joint policy.

Implements the meta-cognitive rule of eq.~(11)–(12) in the paper: a
calibrated logistic classifier on ``[ĝ, Stab_t, ψ_t]`` thresholded at
``z*``. When the classifier is not yet fitted (``fusion_ckpt`` empty), the
controller falls back to the stand-alone rules of eq.~(9) and eq.~(10):

* stopping — ``ĝ < ε``
* triggering — ``Stab_t > δ``

Cold start: for ``t < r_min`` we force ``Stab_t = 0`` (paper §4.2.3), so the
triggering rule votes RETRIEVE unconditionally.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from ...config import ConfidenceConfig
from .bocpd import BocpdState
from .predictor import GainPredictor
from .signals import (
    compute_agreement,
    compute_entropy,
    compute_novelty,
    compute_top1_certainty,
)


@dataclass
class ControllerDecision:
    """Joint decision emitted at each RAG round."""

    trigger: str  # "retrieve" | "skip"
    stop: bool
    stab: float
    gain: float
    z: float
    reason: str


@dataclass
class RoundObservation:
    """What the loop feeds into :meth:`ConfidenceController.step` each round."""

    round_index: int  # 1-based ``t``
    posterior: Sequence[float] = field(default_factory=list)
    sub_agent_embeddings: Sequence[Sequence[float]] = field(default_factory=list)
    current_doc_ids: Iterable[str] = field(default_factory=list)
    history_doc_ids: Iterable[str] = field(default_factory=list)
    # Cumulative retrieval cost proxy. The loop populates this with the
    # count of ``search`` tool calls observed on the trace up to and
    # including the just-completed round (parallel fan-out counted as
    # separate calls). The field name is kept for ckpt/back-compat; only
    # the semantics changed from ``session.tokens_spent`` to search-call
    # count. Fusion ckpts trained before this change must be rebuilt.
    mean_cost: float = 0.0
    all_dry: bool = False
    extra_features: Mapping[str, float] = field(default_factory=dict)


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


class ConfidenceController:
    """Round-level controller: consumes RAG-round observations, emits policy.

    ``config.enabled`` gates the whole thing: with ``enabled=False`` the
    controller should not be constructed by the loop at all (the loop keeps
    its fixed ``max_steps`` behaviour). This class is safe to construct even
    when the checkpoints are missing — the stub predictor and BOCPD keep the
    fallback rules well-defined.
    """

    def __init__(self, config: ConfidenceConfig) -> None:
        self.config = config
        self.bocpd = BocpdState(hazard=config.bocpd_hazard, r_min=config.r_min)
        self.predictor = GainPredictor(ckpt_path=config.predictor_ckpt)
        self._fusion_params = self._load_fusion(config.fusion_ckpt)
        self._history_gains: list[float] = []

    @staticmethod
    def _load_fusion(ckpt: str) -> dict[str, Any] | None:
        """Load ``{"alpha": [...], "bias": float, "z_star": float, "meta": {...}}``
        from the pickle produced by ``scripts/grid_search_fusion.py``.

        Returns ``None`` when the path is empty, missing, or the file fails
        to unpickle — the controller falls back to the stand-alone rules
        (ε for stop, δ for trigger) in that case. This matches the paper's
        cold-start policy: no fusion params → stand-alone.
        """
        if not ckpt:
            return None
        import logging
        import pickle
        p = Path(ckpt)
        if not p.exists():
            logging.getLogger(__name__).warning(
                "fusion ckpt %s not found — falling back to stand-alone rules", p,
            )
            return None
        try:
            with p.open("rb") as f:
                bundle = pickle.load(f)
        except Exception as e:
            logging.getLogger(__name__).warning(
                "failed to load fusion ckpt %s: %s — falling back", p, e,
            )
            return None
        # Basic schema validation.
        alpha = bundle.get("alpha")
        bias = bundle.get("bias")
        if not isinstance(alpha, (list, tuple)) or bias is None:
            logging.getLogger(__name__).warning(
                "fusion ckpt %s missing alpha/bias fields — falling back", p,
            )
            return None
        logging.getLogger(__name__).info(
            "loaded fusion ckpt %s (alpha=%s bias=%s meta=%s)",
            p, list(alpha), bias, bundle.get("meta"),
        )
        return bundle

    def reset(self) -> None:
        self.bocpd.reset()
        self._history_gains.clear()

    def step(self, obs: RoundObservation) -> ControllerDecision:
        entropy = compute_entropy(obs.posterior)
        agreement = compute_agreement(obs.sub_agent_embeddings)
        novelty = compute_novelty(obs.current_doc_ids, obs.history_doc_ids)
        # BOCPD sees a single scalar per round: entropy-derived certainty on
        # the top-B posterior. Agreement/novelty stay in the predictor's
        # feature dict below rather than entering the changepoint model.
        certainty = compute_top1_certainty(obs.posterior, top_b=self.config.top_b)

        stab_raw = self.bocpd.update(certainty)
        # Cold start: force Stab_t=0 for t < r_min (paper §4.2.3).
        if obs.round_index < self.config.r_min:
            stab = 0.0
        else:
            stab = stab_raw

        # Feature dict MUST match the 10 keys the predictor was trained on
        # (see collect_predictor_data.extract_features / predictor.DEFAULT_FEATURE_KEYS).
        # Reuse extract_features directly so online and training features
        # stay in lock-step: t, entropy, top1_prob, posterior_var,
        # n_unique_docs, novelty, gain_mean, gain_var, gain_slope, exhaustion.
        from .collect_predictor_data import extract_features
        features = extract_features(
            round_index=int(obs.round_index),
            posterior=list(obs.posterior),
            current_doc_ids=list(obs.current_doc_ids),
            history_doc_ids=list(obs.history_doc_ids),
            past_gains=list(self._history_gains),
            expected_docs_per_round=int(self.config.expected_docs_per_round),
        )
        # Auxiliary signals fusion still wants downstream (BOCPD stab enters
        # via psi, not features; agreement / pit / top1_certainty are exposed
        # for extra_features overrides only — the predictor itself ignores
        # any keys outside its trained schema).
        features["agreement"] = agreement
        features["top1_certainty"] = certainty
        features["pit"] = self.bocpd.pit()
        features.update(obs.extra_features)
        gain = self.predictor.predict(features)

        # Stand-alone rules (used when fusion checkpoint is absent).
        trig = "skip" if stab > self.config.delta else "retrieve"
        stop_rule = gain < self.config.epsilon

        # Fusion classifier (deferred; see T10).
        if self._fusion_params is not None:
            psi = self._auxiliary_features(obs, stab, gain)
            z = _sigmoid(self._linear(psi))
            reason = "fusion"
            # z_star from the trained ckpt takes precedence over the config
            # default — the pickle written by scripts/grid_search_fusion.py
            # stores the operating point chosen on the frontier, and
            # ConfidenceConfig.z_star is only a cold-start fallback (0.5).
            # An explicit `z_star_override` (ablation/sweep) wins over both.
            z_star = self._fusion_params.get("z_star", self.config.z_star)
            override = getattr(self.config, "z_star_override", None)
            if override is not None:
                z_star = override
            stop = z > z_star
            if stop:
                trig = "skip"
        else:
            z = 0.0
            reason = "fallback"
            stop = stop_rule

        # Safety net: when every sub-agent returned dry, don't stop — the
        # salvage layer downstream may still recover evidence next round.
        if obs.all_dry and stop:
            stop = False
            reason = f"{reason}+dry-block"

        # Cold-start: force RETRIEVE for the first r_min rounds (paper §4.2.3
        # + config comment "force RETRIEVE for the first r_min rounds"). Only
        # clamping stab=0 above isn't enough — the predictor MLP is untrained
        # on t=1 features (collector data starts at t=2), so gain(t=1) is an
        # out-of-distribution output that varies wildly across predictor
        # seeds. Blocking stop at t<r_min makes cold-start behaviour
        # deterministic across seeds and matches the config-level contract.
        if obs.round_index < self.config.r_min and stop:
            stop = False
            reason = f"{reason}+cold-start-block"

        # Persist realised gain proxy for future rolling stats (populated by
        # the loop once it computes KL(p_{t+1}||p_t)); the stub tracks the
        # predicted gain instead, which is stable enough for tests.
        self._history_gains.append(gain)

        return ControllerDecision(
            trigger=trig,
            stop=stop,
            stab=stab,
            gain=gain,
            z=z,
            reason=reason,
        )

    def _auxiliary_features(
        self,
        obs: RoundObservation,
        stab: float,
        gain: float,
    ) -> tuple[float, ...]:
        return (
            gain,
            stab,
            float(obs.round_index),
            float(obs.mean_cost),
            1.0 if obs.all_dry else 0.0,
        )

    def _linear(self, psi: tuple[float, ...]) -> float:
        params = self._fusion_params or {}
        alpha: Sequence[float] = params.get("alpha", (0.0,) * len(psi))
        bias: float = params.get("bias", 0.0)
        return sum(a * x for a, x in zip(alpha, psi)) + bias
