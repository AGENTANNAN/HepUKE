"""Bayesian Online Changepoint Detection (Adams & MacKay, 2007).

Beta-conjugate observation model on a one-dimensional certainty signal
``x_t ∈ [0,1]`` (see :func:`hepuke.agent.confidence.signals.compute_top1_certainty`).
Emits a stability scalar per round plus a probability-integral-transform (PIT)
value computed against the mixture predictive.

Recursion (per round t)::

    joint growth :  P(r_t = r+1, x_{1:t}) = P(r_{t-1}=r, x_{1:t-1}) * f(x_t | α_r, β_r) * (1 - H)
    joint reset  :  P(r_t = 0, x_{1:t})   = Σ_r P(r_{t-1}=r, x_{1:t-1}) * f(x_t | α_r, β_r) * H
    normalize joint → P(r_t | x_{1:t})     (the run-length posterior)
    posterior    :  α_r ← α_r + x_t
                     β_r ← β_r + (1 - x_t)      (Beta conjugate update)

Stability is the tail probability that the current run has already lasted
``r_min`` rounds::

    Stab_t = P(r_t ≥ r_min | x_{1:t})   →   fires the look-back stop signal

The mixture predictive PIT is::

    u_t = Σ_r P(r_{t-1}=r) · Beta_CDF(x_t | α_r, β_r)

integrated over the run-length posterior *before* the current observation is
absorbed. ``u_t`` is Uniform[0,1] under the null (model is well-specified,
observations are stationary). Large deviations from Uniform signal
mis-specification or an ongoing changepoint. ``u_t`` is retained so downstream
consumers (fusion, diagnostics) can use it without re-computing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence


def _log_beta(a: float, b: float) -> float:
    """log B(a, b) = lgamma(a) + lgamma(b) - lgamma(a+b)."""
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _beta_logpdf(x: float, a: float, b: float) -> float:
    """Beta(a, b) log density at ``x`` ∈ (0,1). Endpoints clipped for stability."""
    if x <= 0.0:
        x = 1e-12
    elif x >= 1.0:
        x = 1.0 - 1e-12
    return (a - 1.0) * math.log(x) + (b - 1.0) * math.log(1.0 - x) - _log_beta(a, b)


def _beta_pdf(x: float, a: float, b: float) -> float:
    return math.exp(_beta_logpdf(x, a, b))


def _beta_cdf(x: float, a: float, b: float) -> float:
    """Beta(a,b) CDF via the regularised incomplete beta function."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    # scipy has a stable implementation; keeping the import local so importing
    # this module doesn't require scipy at collection time.
    from scipy.special import betainc  # type: ignore
    return float(betainc(a, b, x))


@dataclass
class BocpdState:
    """Rolling BOCPD state fed one scalar observation ``x_t ∈ [0,1]`` per round.

    Parameters
    ----------
    hazard:
        Geometric hazard ``H = 1/λ_h``. Higher values make the model more
        willing to declare a changepoint.
    r_min:
        Minimum run-length considered "stable"; :meth:`stability` returns
        the tail ``P(r_t ≥ r_min | x_{1:t})``.
    prior_alpha, prior_beta:
        Beta prior for a freshly opened run (``r = 0``). ``(1, 1)`` = uniform.
    """

    hazard: float = 0.1
    r_min: int = 2
    prior_alpha: float = 1.0
    prior_beta: float = 1.0

    # Internal state
    _r_probs: list[float] = field(default_factory=list)
    _alphas: list[float] = field(default_factory=list)
    _betas: list[float] = field(default_factory=list)
    _last_stability: float = 0.0
    _last_pit: float = 0.5  # Uniform mean → neutral

    def __post_init__(self) -> None:
        if not self._r_probs:
            self._initialise()

    def _initialise(self) -> None:
        # Before any observation, r=0 with probability 1; parameters at prior.
        self._r_probs = [1.0]
        self._alphas = [self.prior_alpha]
        self._betas = [self.prior_beta]
        self._last_stability = 0.0
        self._last_pit = 0.5

    # ------------------------------------------------------------------ API
    def update(self, s_t) -> float:
        """Ingest one observation ``x_t`` and return the current stability.

        ``s_t`` may be:
        * a scalar in [0,1] — used directly (the T9+ default), or
        * a sequence — the first entry is used as the scalar (legacy: the
          controller used to pass a ``(entropy, agreement, novelty)`` tuple).
        """
        # Backward-compat: accept either a scalar or a sequence.
        if isinstance(s_t, (int, float)):
            x = float(s_t)
        else:
            try:
                x = float(next(iter(s_t)))
            except (StopIteration, TypeError):
                x = 0.5
        # Clip observation into (0,1) to keep Beta density finite.
        if x < 0.0:
            x = 0.0
        elif x > 1.0:
            x = 1.0

        # --- PIT computed BEFORE absorbing this observation
        # u_t = Σ_r P(r_{t-1}=r) * Beta_CDF(x_t | α_r, β_r)
        pit = 0.0
        for p, a, b in zip(self._r_probs, self._alphas, self._betas):
            if p <= 0.0:
                continue
            pit += p * _beta_cdf(x, a, b)
        self._last_pit = max(0.0, min(1.0, pit))

        # --- Predictive density per run-length (for run-length update)
        pdfs = [_beta_pdf(x, a, b) for a, b in zip(self._alphas, self._betas)]

        # --- Growth / change-point mass
        H = self.hazard
        # growth[r+1] = r_probs[r] * pdf[r] * (1 - H)
        growth = [p * f * (1.0 - H) for p, f in zip(self._r_probs, pdfs)]
        # change[0] = Σ_r r_probs[r] * pdf[r] * H
        change = sum(p * f * H for p, f in zip(self._r_probs, pdfs))

        # --- Compose the new (unnormalised) joint
        new_joint: list[float] = [change] + growth  # index 0 = r=0
        total = sum(new_joint)
        if total <= 0.0:
            # numerical underflow — restart from prior. Rare with the log-pdf
            # path, but a hard reset keeps things well-defined.
            self._initialise()
            self._last_stability = 0.0
            return 0.0
        new_r_probs = [v / total for v in new_joint]

        # --- Update Beta parameters
        # new run r=0 starts fresh at prior; old runs r → r+1 absorb x_t.
        new_alphas: list[float] = [self.prior_alpha]
        new_betas: list[float] = [self.prior_beta]
        for a, b in zip(self._alphas, self._betas):
            new_alphas.append(a + x)
            new_betas.append(b + (1.0 - x))

        # Optional truncation: keep only run-lengths whose posterior mass is
        # non-negligible. Guards against the arrays growing linearly with t
        # (they only do in principle; in practice the tail decays fast).
        THRESH = 1e-6
        kept: list[int] = [i for i, p in enumerate(new_r_probs) if p >= THRESH]
        if not kept:
            kept = [int(max(range(len(new_r_probs)), key=lambda i: new_r_probs[i]))]
        r_probs = [new_r_probs[i] for i in kept]
        alphas = [new_alphas[i] for i in kept]
        betas = [new_betas[i] for i in kept]
        # Re-normalise after truncation.
        s = sum(r_probs)
        if s > 0.0:
            r_probs = [v / s for v in r_probs]
        self._r_probs, self._alphas, self._betas = r_probs, alphas, betas

        # --- Stability: mass on runs of length ≥ r_min
        # Reconstruct actual run-lengths from `kept` indices (which reference
        # positions in the pre-truncation array; index i ↔ r = i).
        stab = 0.0
        for i, p in zip(kept, r_probs):
            if i >= self.r_min:
                stab += p
        self._last_stability = max(0.0, min(1.0, stab))
        return self._last_stability

    def stability(self) -> float:
        """Posterior tail probability ``P(r_t ≥ r_min | x_{1:t})``."""
        return self._last_stability

    def pit(self) -> float:
        """Last computed PIT ``u_t = Σ_r P(r_{t-1}=r) · Beta_CDF(x_t | α_r, β_r)``.

        Neutral value ``0.5`` before any observation. Under a well-specified
        model on stationary data the sequence ``{u_t}`` is Uniform[0,1].
        """
        return self._last_pit

    def reset(self) -> None:
        self._initialise()
