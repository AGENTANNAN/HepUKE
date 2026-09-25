"""Recompute ``g_next`` labels on a collected trajectory JSONL using JSD.

The initial run of ``collect_hotpotqa_trajectories.py`` used
``kl_divergence`` which is unbounded on one-hot posteriors (a HotpotQA
regularity — the parent LLM's top-B collapses to ``[1.0, 0, ..., 0]``
whenever it's confident about the next token). We observed labels up to
8.18 on rows whose φ was numerically identical, wrecking the KL head.

This script rewrites ``g_next`` using
:func:`hepuke.agent.confidence.collect_predictor_data.js_divergence`.

The transform is:

  * per-query, order rows by ``t``
  * for row ``i`` with successor ``i+1``: ``g_next_new = JSD(P_{t+1}, P_t)``
  * last row per query: ``g_next = None`` (unchanged)

The input JSONL is left alone; a new file is written next to it. Optional
``--filter-onehot`` also drops rows whose posterior is ``one-hot`` (entropy
< threshold AND top-1 > threshold) — the JSD-bounded label is still
statistically degenerate there because P_t is a delta.

Usage::

    python scripts/recompute_g_next.py \\
        --in  storage/predictor_train/hotpotqa_train.gpt4omini.jsonl \\
        --out storage/predictor_train/hotpotqa_train.gpt4omini.jsd.jsonl \\
        --filter-onehot

Diagnostic mode (does not write, prints label distribution shift)::

    python scripts/recompute_g_next.py --in <path> --dry-run
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root

from hepuke.agent.confidence.collect_predictor_data import js_divergence


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="input_path", required=True, type=Path)
    p.add_argument("--out", dest="output_path", type=Path,
                   help="Output JSONL. Defaults to <in>.jsd.jsonl next to input.")
    p.add_argument("--filter-onehot", action="store_true",
                   help="Drop rows whose posterior is essentially one-hot "
                        "(entropy < --onehot-entropy AND top1 > --onehot-top1).")
    p.add_argument("--onehot-entropy", type=float, default=0.01,
                   help="Entropy threshold below which a row counts as one-hot.")
    p.add_argument("--onehot-top1", type=float, default=0.99,
                   help="Top-1 prob above which a row counts as one-hot.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the label-distribution shift and exit.")
    return p.parse_args()


def _entropy(p: list[float]) -> float:
    h = 0.0
    total = 0.0
    for x in p:
        if x > 0.0:
            h -= x * math.log(x)
            total += x
    return h if total > 0.0 else 0.0


def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    log = logging.getLogger("recompute_g_next")

    out_path: Path
    if args.output_path is not None:
        out_path = args.output_path
    else:
        out_path = args.input_path.with_suffix(".jsd.jsonl")

    # ---- Load and group by qid, preserving row order per query.
    rows_by_qid: dict[str, list[dict]] = defaultdict(list)
    total = 0
    with args.input_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            qid = r.get("query_id", "?")
            rows_by_qid[qid].append(r)
            total += 1
    log.info("loaded %d rows across %d queries", total, len(rows_by_qid))

    # ---- Recompute g_next per query using JSD.
    old_g: list[float] = []
    new_g: list[float] = []
    onehot_flagged = 0
    for qid, rlist in rows_by_qid.items():
        rlist.sort(key=lambda r: r.get("t", 0))
        for i, r in enumerate(rlist):
            old = r.get("g_next")
            if old is not None and math.isfinite(old):
                old_g.append(float(old))
            if i + 1 < len(rlist):
                p_next = rlist[i + 1].get("posterior") or []
                p_curr = r.get("posterior") or []
                new = js_divergence(p_next, p_curr)
                r["g_next"] = new
                new_g.append(new)
            else:
                r["g_next"] = None
            # One-hot flag (side-channel; not applied unless --filter-onehot)
            post = r.get("posterior") or []
            top1 = max(post) if post else 0.0
            ent = _entropy(post)
            r["_onehot"] = (ent < args.onehot_entropy and top1 > args.onehot_top1)
            if r["_onehot"]:
                onehot_flagged += 1

    # ---- Stats
    def _stats(xs: list[float], name: str) -> None:
        if not xs:
            log.info("[%s] empty", name); return
        log.info(
            "[%s] n=%d  min=%.4f  max=%.4f  mean=%.4f  med=%.4f  std=%.4f",
            name, len(xs), min(xs), max(xs),
            statistics.mean(xs), statistics.median(xs),
            statistics.stdev(xs) if len(xs) > 1 else 0.0,
        )
    _stats(old_g, "OLD g_next (KL)")
    _stats(new_g, "NEW g_next (JSD)")
    log.info("one-hot rows flagged: %d/%d = %.1f%%",
             onehot_flagged, total, onehot_flagged / total * 100 if total else 0)

    if args.dry_run:
        log.info("dry-run: not writing")
        return 0

    # ---- Write output; strip _onehot marker but honor --filter-onehot.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    filtered_onehot = 0
    with out_path.open("w", encoding="utf-8") as f:
        for qid, rlist in rows_by_qid.items():
            for r in rlist:
                is_onehot = r.pop("_onehot", False)
                if args.filter_onehot and is_onehot:
                    filtered_onehot += 1
                    continue
                f.write(json.dumps(r, ensure_ascii=False))
                f.write("\n")
                kept += 1
    log.info(
        "wrote %s: %d rows (dropped %d one-hot rows)%s",
        out_path, kept, filtered_onehot,
        " [filter_onehot=OFF]" if not args.filter_onehot else "",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
