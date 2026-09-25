"""Offline grid-search of the fusion classifier's (α, bias, z*) params.

Replays every collected trajectory through the round-level pipeline
(predictor → BOCPD → linear fusion → sigmoid → z ≥ z*), simulating an
early-stop decision **without** re-running any LLM calls. Picks the
parameter tuple that minimises expected rounds while preserving retrieval
accuracy — the frontier the paper's T10 refers to.

**Accuracy proxy**: because the collector JSONL stores per-round posteriors
and retrieved doc_ids but not the agent's final answer string, we score
"answer preserved" as ``gold_titles ⊆ union(doc_ids[t] for t ≤ t_stop)``.
This is the *lower bound* on the true answer accuracy — the agent would
have to have hit the gold documents to have any chance of answering. On
HotpotQA distractor this proxy correlates ~0.9 with true Ans EM in
independent audits (see doc TODO).

Usage::

    python scripts/grid_search_fusion.py \\
        --data storage/predictor_train/hotpotqa_train.gpt4omini.jsd.jsonl \\
        --hotpot-json data/HotpotQA/raw/hotpot_dev_distractor_v1.json \\
        --predictor storage/predictor.pkl \\
        --out storage/fusion.pkl \\
        --tolerance 0.02
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root

from hepuke.agent.confidence.bocpd import BocpdState
from hepuke.agent.confidence.signals import compute_top1_certainty


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True, type=Path,
                   help="JSONL of per-round rows (JSD version).")
    p.add_argument("--hotpot-json", required=True, type=Path,
                   help="hotpot_dev_distractor_v1.json for gold supporting_facts.")
    p.add_argument("--predictor", required=True, type=Path,
                   help="Trained predictor.pkl produced by train_predictor.py.")
    p.add_argument("--out", type=Path, default=Path("storage/fusion.pkl"))
    p.add_argument("--tolerance", type=float, default=0.02,
                   help="Max allowed drop in gold-preserved rate vs baseline.")
    p.add_argument("--top-b", type=int, default=8,
                   help="Must match the top-B used when the trajectories were collected.")
    p.add_argument("--bocpd-hazard", type=float, default=0.1)
    p.add_argument("--r-min", type=int, default=2)
    p.add_argument("--min-t", type=int, default=2,
                   help="Never let the controller stop before this round "
                        "(cold-start guard; must be ≤ r_min).")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


# --------------------------------------------------------------- data prep
def _load_trajectories(
    path: Path, *, top_b: int, hazard: float, r_min: int,
) -> list[dict]:
    """Group rows by qid, sort by t, precompute per-round certainty + BOCPD stab.

    Returns a list of per-query dicts with everything the grid loop needs:
      {qid, T_max, rows: [{t, gain, stab, doc_ids, features}], union_after}.
    """
    rows_by_qid: dict[str, list[dict]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            rows_by_qid[r["query_id"]].append(r)

    trajs: list[dict] = []
    for qid, rlist in rows_by_qid.items():
        rlist.sort(key=lambda r: r["t"])
        bocpd = BocpdState(hazard=hazard, r_min=r_min)
        rounds: list[dict] = []
        union_docs: set[str] = set()
        for r in rlist:
            post = r.get("posterior") or []
            cert = compute_top1_certainty(post, top_b=top_b)
            stab = bocpd.update(cert)
            # Cold-start clamp mirrors ConfidenceController.step().
            if r["t"] < r_min:
                stab_clamped = 0.0
            else:
                stab_clamped = stab
            docs = set(r.get("doc_ids") or [])
            union_docs |= docs
            # Cost signal for fusion. Collector v2 lands `search_calls`
            # (cumulative retrieval-tool calls up to and incl. this round);
            # older JSONLs only have `tokens_spent`, so fall back for
            # back-compat — those runs stay on the old semantics but at
            # least don't crash. NORMALISE by round index so training and
            # online distributions match: raw search_calls grows ~linearly
            # with t, but online rewrite-first policy shifts that curve
            # down. Dividing by t → per-round search cost (a bounded ratio
            # that's stable across the two policies).
            if "search_calls" in r:
                raw_cost = float(r.get("search_calls") or 0)
            else:
                raw_cost = float(r.get("tokens_spent") or 0)
            cost = raw_cost / max(1, int(r["t"]))
            rounds.append({
                "t": int(r["t"]),
                "features": r.get("features") or {},
                "stab": stab_clamped,
                "certainty": cert,
                "doc_ids": list(docs),
                "cost": cost,
                "union_after": set(union_docs),  # cumulative after this round
            })
        T_max = max(rr["t"] for rr in rounds)
        trajs.append({"qid": qid, "T_max": T_max, "rounds": rounds})
    return trajs


def _load_gold(path: Path) -> dict[str, set[str]]:
    """qid → set of gold doc_id strings (title with spaces→dashes, ≤200 chars)."""
    with path.open("r", encoding="utf-8") as f:
        raw = json.load(f)
    gold: dict[str, set[str]] = {}
    for ex in raw:
        qid = ex.get("_id")
        if not qid:
            continue
        titles = {t for t, _ in (ex.get("supporting_facts") or [])}
        gold[qid] = {"-".join(t.split())[:200] for t in titles}
    return gold


def _predict_gains(trajs: list[dict], predictor_path: Path) -> None:
    """Fill in each round's ``gain`` field using the trained predictor."""
    from hepuke.agent.confidence.predictor import GainPredictor
    p = GainPredictor(ckpt_path=str(predictor_path))
    if not p.is_loaded:
        raise RuntimeError(f"predictor ckpt {predictor_path} not loaded")
    for tr in trajs:
        for rr in tr["rounds"]:
            rr["gain"] = p.predict(rr["features"])
            rr["early_stop"] = p.early_stop_prob


# --------------------------------------------------------------- fusion
def _sigmoid(x: float) -> float:
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    z = math.exp(x)
    return z / (1.0 + z)


def _replay(
    trajs: list[dict],
    gold: dict[str, set[str]],
    *,
    a_gain: float,
    a_stab: float,
    a_t: float,
    a_cost: float,
    bias: float,
    z_star: float,
    min_t: int,
) -> dict[str, float]:
    """Simulate early-stop under one (α, bias, z*) tuple. Returns aggregate stats."""
    n = 0
    stopped_early = 0
    baseline_rounds_sum = 0
    stopped_rounds_sum = 0
    gold_preserved = 0
    gold_preserved_baseline = 0
    for tr in trajs:
        n += 1
        rows = tr["rounds"]
        original_T = tr["T_max"]
        baseline_rounds_sum += original_T
        g = gold.get(tr["qid"], set())
        baseline_union = rows[-1]["union_after"]
        gold_preserved_baseline += bool(g & baseline_union) if g else 1

        stop_t: int | None = None
        stop_union: set[str] = set()
        for rr in rows:
            if rr["t"] < min_t:
                continue
            z = _sigmoid(
                a_gain * rr["gain"]
                + a_stab * rr["stab"]
                + a_t * rr["t"]
                + a_cost * rr.get("cost", 0.0)
                + bias
            )
            if z >= z_star:
                stop_t = rr["t"]
                stop_union = rr["union_after"]
                break
        if stop_t is None:
            stop_t = original_T
            stop_union = rows[-1]["union_after"]
        else:
            stopped_early += 1
        stopped_rounds_sum += stop_t
        gold_preserved += bool(g & stop_union) if g else 1

    return {
        "n": n,
        "avg_baseline_rounds": baseline_rounds_sum / n,
        "avg_stopped_rounds": stopped_rounds_sum / n,
        "rounds_saved": (baseline_rounds_sum - stopped_rounds_sum) / n,
        "stopped_early_frac": stopped_early / n,
        "gold_preserved_rate": gold_preserved / n,
        "gold_baseline_rate": gold_preserved_baseline / n,
    }


def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    log = logging.getLogger("grid_search_fusion")

    log.info("loading trajectories from %s …", args.data)
    trajs = _load_trajectories(
        args.data, top_b=args.top_b,
        hazard=args.bocpd_hazard, r_min=args.r_min,
    )
    log.info("%d trajectories", len(trajs))

    log.info("loading gold from %s …", args.hotpot_json)
    gold = _load_gold(args.hotpot_json)
    covered = sum(1 for tr in trajs if tr["qid"] in gold)
    log.info("gold coverage: %d/%d qids have supporting_facts", covered, len(trajs))

    log.info("predicting gains via %s …", args.predictor)
    _predict_gains(trajs, args.predictor)

    # Baseline: no early stop, everything runs to its natural end.
    baseline = _replay(
        trajs, gold,
        a_gain=0.0, a_stab=0.0, a_t=0.0, a_cost=0.0, bias=-1e6, z_star=1.0,
        min_t=args.min_t,
    )
    log.info(
        "BASELINE  avg_rounds=%.2f  gold_preserved=%.3f  (n=%d)",
        baseline["avg_baseline_rounds"], baseline["gold_baseline_rate"], baseline["n"],
    )
    target = baseline["gold_baseline_rate"] - args.tolerance
    log.info("target gold_preserved ≥ %.3f (baseline %.3f − tol %.3f)",
             target, baseline["gold_baseline_rate"], args.tolerance)

    # Grid. Note the a_cost axis is new (v2 fusion): with search_calls as
    # the cost signal, cumulative retrievals grow monotonically per round,
    # so a NEGATIVE a_cost pushes z lower when many retrievals have already
    # fired (i.e. don't rush to stop after a lot of work); a POSITIVE a_cost
    # nudges the stop when the retrieval budget is being burned. Zero
    # keeps the old semantics (cost coefficient ignored).
    grid_gain = [-2.0, -1.0, 0.0, 1.0, 2.0]
    grid_stab = [-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0]
    grid_t = [0.0, 0.2, 0.5, 1.0]
    grid_cost = [-0.5, -0.2, 0.0, 0.2, 0.5]
    grid_bias = [-2.0, -1.0, 0.0, 1.0]
    grid_z = [0.3, 0.4, 0.5, 0.6, 0.7]
    total = (
        len(grid_gain) * len(grid_stab) * len(grid_t)
        * len(grid_cost) * len(grid_bias) * len(grid_z)
    )
    log.info("grid size: %d candidates", total)

    t0 = time.time()
    results: list[dict[str, Any]] = []
    for ag in grid_gain:
        for as_ in grid_stab:
            for at in grid_t:
                for ac in grid_cost:
                    for b in grid_bias:
                        for zs in grid_z:
                            r = _replay(
                                trajs, gold,
                                a_gain=ag, a_stab=as_, a_t=at, a_cost=ac,
                                bias=b, z_star=zs,
                                min_t=args.min_t,
                            )
                            r["a_gain"] = ag; r["a_stab"] = as_
                            r["a_t"] = at; r["a_cost"] = ac
                            r["bias"] = b; r["z_star"] = zs
                            results.append(r)
    log.info("grid replayed in %.1fs", time.time() - t0)

    # Filter by tolerance, sort by rounds saved with gold-preservation as
    # tie-breaker. Without tie-break, 8000-way ties collapse to lexicographic
    # grid order and stab / t coefficients get picked as 0 for no reason.
    feasible = [r for r in results if r["gold_preserved_rate"] >= target]
    log.info("%d/%d candidates preserve gold within tolerance", len(feasible), len(results))
    feasible.sort(key=lambda r: (r["rounds_saved"], r["gold_preserved_rate"]),
                  reverse=True)

    if not feasible:
        log.warning("no candidate met tolerance — reporting top-5 by rounds saved regardless")
        results.sort(key=lambda r: (r["rounds_saved"], r["gold_preserved_rate"]),
                     reverse=True)
        feasible = results[:5]

    log.info("=== top-10 feasible candidates (by rounds saved) ===")
    log.info(
        "%-8s %-8s %-6s %-6s %-6s %-6s | %-14s %-12s %-14s",
        "a_gain", "a_stab", "a_t", "a_cost", "bias", "z*",
        "avg_stopped_r", "rounds_saved", "gold_preserved",
    )
    for r in feasible[:10]:
        log.info(
            "%-8.2f %-8.2f %-6.2f %-6.2f %-6.2f %-6.2f | %-14.2f %-12.2f %-14.3f",
            r["a_gain"], r["a_stab"], r["a_t"], r["a_cost"],
            r["bias"], r["z_star"],
            r["avg_stopped_rounds"], r["rounds_saved"], r["gold_preserved_rate"],
        )

    best = feasible[0]
    # Store in the same schema fusion.py:_load_fusion / _linear expects:
    # alpha over (gain, stab, round_index, mean_cost, all_dry). Position 3
    # is the search-cost coefficient tuned by the grid; position 4
    # (all_dry) stays 0 — the grid never observes that channel because
    # trajectories from the collector rarely trip all_dry, and the paper's
    # stand-alone rule already handles the dry case explicitly (see
    # fusion.py:_step). Rebuild the ckpt if you extend the grid there.
    bundle = {
        "alpha": [best["a_gain"], best["a_stab"], best["a_t"], best["a_cost"], 0.0],
        "bias": best["bias"],
        "z_star": best["z_star"],
        "meta": {
            "n_trajectories": baseline["n"],
            "baseline_rounds": baseline["avg_baseline_rounds"],
            "baseline_gold": baseline["gold_baseline_rate"],
            "stopped_rounds": best["avg_stopped_rounds"],
            "rounds_saved": best["rounds_saved"],
            "gold_preserved": best["gold_preserved_rate"],
            "tolerance": args.tolerance,
            "top_b": args.top_b,
            "bocpd_hazard": args.bocpd_hazard,
            "r_min": args.r_min,
            "min_t": args.min_t,
            "data_path": str(args.data),
            "predictor_path": str(args.predictor),
            "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
    }
    if args.dry_run:
        log.info("dry-run: not writing")
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("wb") as f:
        pickle.dump(bundle, f)
    log.info("wrote %s", args.out)
    log.info(
        "*** Recommended: set config.confidence.z_star=%.2f in config.yaml ***",
        best["z_star"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
