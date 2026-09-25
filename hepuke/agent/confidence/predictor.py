"""Look-ahead KL-gain predictor ``ĝ_θ(φ(H_t))``.

Loads the pickle bundle produced by ``scripts/train_predictor.py`` (see
`plan_confidence_controller.md T8`) and exposes two heads:

* **KL head** — ``MLPRegressor`` predicting ``KL(p_{t+1} || p_t)`` from φ,
  followed by an ``IsotonicRegression`` calibrator so the raw output is
  monotonically mapped to a comparable scale. This is the primary ``ĝ``.
* **Early-stop head** — ``LogisticRegression`` predicting the probability
  that the current φ belongs to a one-shot terminating trajectory. Exposed
  as :attr:`early_stop_prob` for the fusion layer to combine with ``ĝ``.

When no checkpoint is present (empty path or missing file) the predictor
falls back to a constant ``default_gain`` so the stand-alone stopping rule
never fires prematurely — this is the safe cold-start behaviour used before
any data has been collected.
"""

from __future__ import annotations

import logging
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

log = logging.getLogger(__name__)


# 9-dim φ feature order MUST match ``scripts/train_predictor.py:FEATURE_KEYS``
# and ``collect_predictor_data.extract_features``. Any drift here silently
# feeds the MLP mis-ordered inputs.
DEFAULT_FEATURE_KEYS: tuple[str, ...] = (
    "t",
    "entropy",
    "top1_prob",
    "posterior_var",
    "n_unique_docs",
    "novelty",
    "gain_mean",
    "gain_var",
    "gain_slope",
    "exhaustion",
)


@dataclass
class GainPredictor:
    """Predict the expected KL gain of the next retrieval round.

    ``ckpt_path`` points at the pickle bundle from ``scripts/train_predictor.py``.
    When empty or missing, :meth:`predict` returns ``default_gain`` (chosen
    high so the stand-alone stop rule never fires without a trained model).
    """

    ckpt_path: str = ""
    default_gain: float = 1.0
    _feature_keys: tuple[str, ...] = field(default_factory=lambda: DEFAULT_FEATURE_KEYS)
    _kl_mlp: Any = None
    _kl_isotonic: Any = None
    _early_stop_clf: Any = None
    _meta: dict = field(default_factory=dict)
    _loaded: bool = False
    _last_early_stop_prob: float = 0.0

    def __post_init__(self) -> None:
        if not self.ckpt_path:
            return
        path = Path(self.ckpt_path)
        if not path.exists():
            log.warning("predictor ckpt %s does not exist — falling back to default_gain", path)
            return
        try:
            with path.open("rb") as f:
                bundle = pickle.load(f)
        except Exception as e:
            log.warning("failed to load predictor ckpt %s: %s — falling back", path, e)
            return
        keys = bundle.get("feature_keys")
        if keys:
            self._feature_keys = tuple(keys)
        self._kl_mlp = bundle.get("kl_mlp")
        self._kl_isotonic = bundle.get("kl_isotonic")
        self._early_stop_clf = bundle.get("early_stop_clf")
        self._meta = dict(bundle.get("meta") or {})
        self._loaded = self._kl_mlp is not None
        log.info(
            "loaded predictor ckpt from %s (kl_mlp=%s isotonic=%s early_stop=%s meta=%s)",
            path,
            self._kl_mlp is not None,
            self._kl_isotonic is not None,
            self._early_stop_clf is not None,
            self._meta,
        )

    @property
    def is_loaded(self) -> bool:
        return self._loaded

    @property
    def early_stop_prob(self) -> float:
        """Probability of one-shot termination from the auxiliary classifier.

        Set on every :meth:`predict` call. ``0.0`` when the auxiliary head is
        unavailable (checkpoint missing or classifier not fit).
        """
        return self._last_early_stop_prob

    def _to_vector(self, features: Mapping[str, float]) -> list[float]:
        """Assemble φ in the exact order the models expect. Missing → 0.0."""
        vec: list[float] = []
        for k in self._feature_keys:
            try:
                v = float(features.get(k, 0.0) or 0.0)
            except (TypeError, ValueError):
                v = 0.0
            vec.append(v)
        return vec

    def predict(self, features: Mapping[str, float]) -> float:
        """Return calibrated ``ĝ`` (predicted next-round KL gain).

        Features accepted (all optional; missing → 0.0):
          * ``t``, ``entropy``, ``top1_prob``, ``posterior_var``,
            ``n_unique_docs``, ``novelty``, ``gain_mean``, ``gain_var``,
            ``gain_slope``.

        Also computes :attr:`early_stop_prob` as a side-effect so callers
        can pull it without a second predict call.
        """
        vec = self._to_vector(features)

        # Auxiliary head: always runs when present.
        self._last_early_stop_prob = 0.0
        if self._early_stop_clf is not None:
            try:
                # LogisticRegression: predict_proba → [n_samples, n_classes]
                proba = self._early_stop_clf.predict_proba([vec])[0]
                # Class order comes from `.classes_`; we want P(y=1).
                classes = list(getattr(self._early_stop_clf, "classes_", [0, 1]))
                if 1 in classes:
                    self._last_early_stop_prob = float(proba[classes.index(1)])
                else:
                    self._last_early_stop_prob = float(proba[-1])
            except Exception as e:
                log.debug("early-stop head predict failed: %s", e)

        # KL head: return default when not loaded.
        if not self._loaded or self._kl_mlp is None:
            return self.default_gain

        try:
            raw = float(self._kl_mlp.predict([vec])[0])
        except Exception as e:
            log.warning("KL head predict failed: %s — falling back to default_gain", e)
            return self.default_gain

        if self._kl_isotonic is not None:
            try:
                calibrated = float(self._kl_isotonic.transform([raw])[0])
                return calibrated
            except Exception as e:
                log.debug("isotonic transform failed: %s — returning raw MLP output", e)

        return raw
