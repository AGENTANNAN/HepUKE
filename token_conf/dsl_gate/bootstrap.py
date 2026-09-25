"""Two-level (target outer, sample inner) paired bootstrap.

Estimates the sampling variance of every strategy on every metric while
accounting for BOTH task-level and generator-sample-level variance.

Inputs
------
--cand-file     runs_gate/<tag>.jsonl           (generate_candidates summary rows per cand)
--baseline-file runs_gate/<tag>_baseline.jsonl  (generate_candidates summary rows per baseline)
--raw-file      runs_gate/<tag>_raw.jsonl.gz    (per-sample raw answers)
--judgments     results/phys_acc100_<tag>/judgments.jsonl  (judge_phys_acc per-sample verdicts)

Protocol
--------
Every qid has, at each stage, N samples (default N=5). One bootstrap iteration:
  1. resample qids with replacement (size = original qid count)
  2. for each resampled qid, resample its sample_idx set with replacement
  3. for every strategy, pick the candidate that strategy would pick on that
     qid (candidate identity is fixed per qid × strategy — the pick uses generate_candidates
     summary rows, not per-sample rows, because the entropy signal is a
     per-cand mean already); then average the picked candidate's per-sample
     metric across the resampled sample_idx set
  4. average across resampled qids → one bootstrapped mean per strategy
Repeat B iterations, report:
  - mean, 2.5 / 97.5 percentiles of each strategy's mean
  - two-sided percentile CI of Δ = strategy − baseline, per baseline
    (baseline = no_rag AND baseline = reranker_top1, both directions kept)

Metrics
-------
struct_mean / strict_struct — recomputed per (qid, sample_idx) with evaluator.
phys100 — pulled from --judgments cache; needs judgments per sample_idx (judge_phys_acc
--raw-file mode).

Usage
-----
    python token_conf/dsl_gate/bootstrap.py \\
        --cand-file runs_gate/dev20_bare_gold_T01.jsonl \\
        --baseline-file runs_gate/dev20_bare_gold_T01_baseline.jsonl \\
        --raw-file  runs_gate/dev20_bare_gold_T01_raw.jsonl.gz \\
        --judgments results/phys_acc100_dev20_bare_gold_T01/judgments.jsonl \\
        --seed 20260922 --iters 10000 \\
        --out-md results/bootstrap_dev20_bare_gold_T01.md
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

import pipeline_common as pc                     # noqa: E402
from pipeline_common import strip_code_fence      # noqa: E402
from evaluator import score_pair              # noqa: E402

GOLD = ROOT / "data" / "corpus" / "generated_dsl"

STRAT_ALPHAS = [0.0, 0.3, 0.5, 0.7, 1.0]


# --------------------------- loading -------------------------------------

def _load_cand(path: Path) -> dict[str, list[dict]]:
    out = defaultdict(list)
    for l in path.open():
        r = json.loads(l)
        out[r["qid"]].append(r)
    return out


def _load_raw(path: Path) -> dict[tuple, dict[int, str]]:
    """key = (qid, kind, cand_stem) -> {sample_idx: ruby_text}"""
    out: dict[tuple, dict[int, str]] = defaultdict(dict)
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        for l in opener(path, "rt", encoding="utf-8"):
            r = json.loads(l)
            key = (r["qid"], r["kind"], r.get("cand_stem"))
            out[key][int(r["sample_idx"])] = strip_code_fence(r.get("answer", "") or "")
    except (EOFError, OSError) as e:
        print(f"[boot] raw file truncated: {type(e).__name__}: {e}", flush=True)
    return out


def _load_phys(path: Path | None) -> dict[tuple, dict[int, float]]:
    """key = (qid, kind, cand_stem) -> {sample_idx: total_score}"""
    out: dict[tuple, dict[int, float]] = defaultdict(dict)
    if path is None or not path.exists():
        return out
    for l in path.open():
        try:
            r = json.loads(l)
        except json.JSONDecodeError:
            continue
        if r.get("alt_judge") or r.get("total") is None:
            continue
        key = (r["qid"], r["kind"], r["cand_stem"])
        out[key][int(r.get("sample_idx", 0))] = float(r["total"])
    return out


# --------------------------- strategy pickers -----------------------------

def _real(cands: list[dict]) -> list[dict]:
    return [c for c in cands if not c.get("is_gold") and c.get("cand_rank") != 0]


def _gold_cand(cands: list[dict]) -> dict | None:
    for c in cands:
        if c.get("is_gold") or c.get("cand_rank") == 0:
            return c
    return None


def pick_stem(strat: str, qid: str, cands: list[dict],
              phys: dict[tuple, dict[int, float]] | None = None) -> tuple[str, str]:
    """Return (kind, cand_stem_or_None) for this strategy.

    kind ∈ {"baseline", "cand"}. For "baseline" the cand_stem is None.
    Picker uses generate_candidates summary fields (dH_full / dH_sel), so identity is stable
    across sample-idx resampling (only the average metric on that identity
    changes)."""
    real = _real(cands)
    if strat == "no_rag":
        return "baseline", None
    if strat == "reranker_top1":
        return "cand", min(real or cands, key=lambda c: c.get("cand_rank", 999))["cand_stem"]
    if strat.startswith("entropy(a="):
        a = float(strat.split("=")[1].rstrip(")"))
        return "cand", max(real or cands,
                           key=lambda c: a * (c.get("dH_full") or 0.0)
                                          + (1 - a) * (c.get("dH_sel") or 0.0))["cand_stem"]
    if strat == "random":
        return "cand", random.Random(42).choice(real or cands)["cand_stem"]
    if strat == "oracle_gold":
        g = _gold_cand(cands)
        return ("cand", g["cand_stem"]) if g else (None, None)
    if strat == "oracle_phys100":
        if not phys:
            return (None, None)
        best, bs = None, -1.0
        for c in (real or cands):
            samples = phys.get((qid, "cand", c["cand_stem"]), {})
            if not samples:
                continue
            m = sum(samples.values()) / len(samples)
            if m > bs:
                bs, best = m, c["cand_stem"]
        return ("cand", best) if best else (None, None)
    raise ValueError(strat)


# --------------------------- per-sample metrics ---------------------------

def _score_metric(ruby: str, gold: str, metric: str) -> float:
    m = score_pair(ruby, gold)
    return m[metric]


def build_metric_table(
    cands_by_qid: dict[str, list[dict]],
    raw: dict[tuple, dict[int, str]],
    phys: dict[tuple, dict[int, float]],
    strategies: list[str],
    metrics: list[str],
) -> dict[tuple, dict[int, float]]:
    """Precompute metric[strategy, qid, sample_idx] once so the inner
    bootstrap loop is pure numeric aggregation.

    Returns: (strategy, qid, metric) -> {sample_idx -> value}
    Missing values fall back to None; strategies that can't pick on a qid
    are absent from the map so bootstrap can skip that qid for that strategy.
    """
    tbl: dict[tuple, dict[int, float]] = {}
    for qid, cands in cands_by_qid.items():
        gold = (GOLD / f"{qid}.rb").read_text(encoding="utf-8")
        for strat in strategies:
            kind, cs = pick_stem(strat, qid, cands, phys=phys)
            if kind is None:
                continue
            sample_map = raw.get((qid, kind, cs), {})
            for metric in metrics:
                key = (strat, qid, metric)
                if metric == "phys100":
                    verdicts = phys.get((qid, kind, cs), {})
                    tbl[key] = dict(verdicts)  # copy
                else:
                    tbl[key] = {si: _score_metric(ruby, gold, metric)
                                for si, ruby in sample_map.items() if ruby}
    return tbl


# --------------------------- bootstrap -----------------------------------

def bootstrap(
    tbl: dict[tuple, dict[int, float]],
    strategies: list[str],
    metrics: list[str],
    qids: list[str],
    B: int,
    seed: int,
) -> dict[tuple, list[float]]:
    """Return (strategy, metric) -> list of B bootstrapped means.

    Pairing preserved: within iteration b, ALL (strategy, metric) share the
    same resampled qid list AND the same per-qid sample_idx resamples. This is
    what lets a downstream Δ(strat, base) = boot[strat,m][b] − boot[base,m][b]
    be a valid paired difference."""
    rng = random.Random(seed)
    out: dict[tuple, list[float]] = {(s, m): [] for s in strategies for m in metrics}
    n = len(qids)
    for _ in range(B):
        # 1. one shared qid resample
        boot_qids = [qids[rng.randrange(n)] for _ in range(n)]
        # 2. one shared sample_idx resample per qid (indices into 0..N-1,
        #    per-qid N inferred from whichever (strat, qid, metric) row is
        #    present; assumes N is same across strategies for a qid, which
        #    holds because raw stores all N samples per (qid, kind))
        sample_boot: dict[str, list[int]] = {}
        for q in boot_qids:
            if q in sample_boot:
                continue
            # find any strategy that has samples for this qid to learn N
            keys = None
            for strat in strategies:
                for metric in metrics:
                    s = tbl.get((strat, q, metric))
                    if s:
                        keys = list(s.keys())
                        break
                if keys:
                    break
            if not keys:
                sample_boot[q] = []
                continue
            k = len(keys)
            sample_boot[q] = [keys[rng.randrange(k)] for _ in range(k)]
        # 3. compute per-(strat, metric) mean under this shared resample
        for strat in strategies:
            for metric in metrics:
                per_qid = []
                for q in boot_qids:
                    samples = tbl.get((strat, q, metric), {})
                    if not samples:
                        continue
                    idxs = sample_boot.get(q, [])
                    if not idxs:
                        continue
                    boot_vals = [samples[i] for i in idxs if i in samples]
                    if not boot_vals:
                        continue
                    per_qid.append(sum(boot_vals) / len(boot_vals))
                if per_qid:
                    out[(strat, metric)].append(sum(per_qid) / len(per_qid))
    return out


def _pct(xs: list[float], p: float) -> float:
    if not xs:
        return float("nan")
    xs2 = sorted(xs)
    i = max(0, min(len(xs2) - 1, int(round(p * (len(xs2) - 1)))))
    return xs2[i]


# --------------------------- main -----------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand-file", type=Path, required=True)
    ap.add_argument("--baseline-file", type=Path, required=True)
    ap.add_argument("--raw-file", type=Path, required=True,
                    help="gzipped per-sample raw jsonl from generate_candidates")
    ap.add_argument("--judgments", type=Path, default=None,
                    help="judge_phys_acc judgments.jsonl (with per-sample scores if available)")
    ap.add_argument("--seed", type=int, default=20260922)
    ap.add_argument("--iters", type=int, default=10000)
    ap.add_argument("--metrics", nargs="+",
                    default=["struct_mean", "strict_struct_mean", "phys100"])
    ap.add_argument("--baselines-vs", nargs="+",
                    default=["no_rag", "reranker_top1"],
                    help="strategies to test Δ AGAINST (two-sided CIs kept)")
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args()

    strategies = (["no_rag", "reranker_top1"]
                  + [f"entropy(a={a:g})" for a in STRAT_ALPHAS]
                  + ["random", "oracle_phys100", "oracle_gold"])

    cands_by_qid = _load_cand(args.cand_file)
    # baseline rows just supply the qid list — raw already carries baseline samples
    base_qids = [json.loads(l)["qid"] for l in args.baseline_file.open()]
    qids = sorted(q for q in cands_by_qid if q in set(base_qids))
    print(f"[boot] qids={len(qids)} strategies={len(strategies)} metrics={args.metrics}")

    raw = _load_raw(args.raw_file)
    phys = _load_phys(args.judgments)
    print(f"[boot] raw groups={len(raw)}  phys groups={len(phys)}")

    tbl = build_metric_table(cands_by_qid, raw, phys, strategies, args.metrics)

    # point estimates: strategy × metric, averaged across sample_idx then qids
    point: dict[tuple, float] = {}
    for strat in strategies:
        for metric in args.metrics:
            per_qid = []
            for q in qids:
                s = tbl.get((strat, q, metric), {})
                if not s:
                    continue
                per_qid.append(sum(s.values()) / len(s))
            if per_qid:
                point[(strat, metric)] = sum(per_qid) / len(per_qid)

    print(f"[boot] running {args.iters} iters ...")
    boot = bootstrap(tbl, strategies, args.metrics, qids, args.iters, args.seed)

    # --- format ---
    md: list[str] = [
        f"# Two-level paired bootstrap (dev n={len(qids)}, seed={args.seed}, iters={args.iters})",
        "",
        "Outer resample: qids with replacement. Inner resample: per-qid, "
        "resample the N generator samples with replacement, then average. "
        "Both variance sources are propagated into every CI below.",
        "",
    ]

    # per-strategy per-metric mean + 95% CI
    md += ["## Point estimate + 95% percentile CI",
           "",
           "| strategy | " + " | ".join(f"{m}" for m in args.metrics) + " |",
           "|" + ("---|" * (1 + len(args.metrics)))]
    for strat in strategies:
        cells = [strat]
        for metric in args.metrics:
            pt = point.get((strat, metric))
            xs = boot.get((strat, metric), [])
            if pt is None or not xs:
                cells.append("—")
                continue
            lo, hi = _pct(xs, 0.025), _pct(xs, 0.975)
            cells.append(f"{pt:.4f} [{lo:.4f}, {hi:.4f}]")
        md.append("| " + " | ".join(cells) + " |")

    # paired Δ CIs — both directions kept, plus one-sided win rate P(strategy > baseline)
    for baseline in args.baselines_vs:
        md += ["", f"## Δ vs {baseline} — two-sided 95% CI + one-sided P(strategy > {baseline})", ""]
        md += ["| strategy | " + " | ".join(f"{m}" for m in args.metrics) + " |",
               "|" + ("---|" * (1 + len(args.metrics)))]
        for strat in strategies:
            if strat == baseline:
                continue
            cells = [strat]
            for metric in args.metrics:
                xs_s = boot.get((strat, metric), [])
                xs_b = boot.get((baseline, metric), [])
                if not xs_s or not xs_b or len(xs_s) != len(xs_b):
                    cells.append("—")
                    continue
                diffs = [a - b for a, b in zip(xs_s, xs_b)]
                lo, hi = _pct(diffs, 0.025), _pct(diffs, 0.975)
                med = _pct(diffs, 0.5)
                p_gt = sum(1 for a, b in zip(xs_s, xs_b) if a > b) / len(diffs)
                p_eq = sum(1 for a, b in zip(xs_s, xs_b) if a == b) / len(diffs)
                sig = "*" if lo > 0 or hi < 0 else " "
                cells.append(f"{med:+.4f} [{lo:+.4f}, {hi:+.4f}]{sig}  P(>)={p_gt:.3f}"
                             + (f" (eq={p_eq:.3f})" if p_eq > 0.02 else ""))
            md.append("| " + " | ".join(cells) + " |")

    md += ["",
           "\\* two-sided 95% percentile CI excludes 0.",
           "P(>) = fraction of the 10 000 paired bootstrap iterations in which the "
           "strategy strictly outperforms the baseline; equivalent to a one-sided "
           "posterior probability. `eq` shown only when ties > 2%.",
           ]
    txt = "\n".join(md)
    print(txt)
    if args.out_md:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text(txt)
        print(f"[boot] wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
