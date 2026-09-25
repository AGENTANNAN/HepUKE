"""Confidence-based retrieval controller.

Implements the meta-cognitive controller from entropy/section_method.tex §4:
three-stream look-back signal (entropy / agreement / novelty) fed through
BOCPD, a look-ahead KL-gain predictor, and a fusion classifier that issues
the joint (trigger, stop) decision.

Public entry point is :class:`ConfidenceController` in ``fusion``.
"""

from .fusion import ConfidenceController, ControllerDecision
from .signals import compute_agreement, compute_entropy, compute_novelty

__all__ = [
    "ConfidenceController",
    "ControllerDecision",
    "compute_agreement",
    "compute_entropy",
    "compute_novelty",
]
