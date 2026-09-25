"""Actionable errors for pipeline stages that have not been run yet.

The analysis scripts consume files produced by ``generate_candidates`` /
``generate_baseline`` / ``judge_phys_acc``. On a fresh clone those files do not exist, and a bare
``FileNotFoundError`` does not tell you which command to run. ``require()``
turns that into a message naming the producing command.
"""
from __future__ import annotations

from pathlib import Path

# file name -> command that produces it
PRODUCERS = {
    "all.jsonl":
        "python token_conf/dsl_gate/generate_candidates.py --split dev "
        "--n-samples 5 --temperature 0.1 --out-tag all",
    "all_baseline.jsonl":
        "python token_conf/dsl_gate/generate_baseline.py --split dev "
        "--n-samples 5 --temperature 0.1 --desc-tier bare "
        "--skill-mode no-example --out-tag all",
    "all_summary.jsonl":
        "python token_conf/dsl_gate/strategy_matrix.py --split dev "
        "--per-qid-out token_conf/dsl_gate/results/",
    "judgments.jsonl":
        "python token_conf/dsl_gate/judge_phys_acc.py --workers 8",
    "splits.json":
        "python token_conf/dsl_gate/splits.py",
    "gold_blocks.json":
        "python scripts/index_besiii_dsl_blocks_v1.py --dry-run "
        "--allow-unresolved-gold",
}


def require(path: Path, produced_by: str | None = None) -> Path:
    """Return ``path``, or raise with the command that produces it."""
    if path.exists():
        return path
    cmd = produced_by or PRODUCERS.get(path.name)
    hint = f"\n  produce it with:\n    {cmd}" if cmd else ""
    raise FileNotFoundError(
        f"required input not found: {path}"
        f"{hint}\n  (run `python scripts/check_setup.py` for a full status list)")
