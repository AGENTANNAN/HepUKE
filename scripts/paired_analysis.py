"""Paired analysis of two (or more) hotpot eval runs on the same seed.

Aligns rows by query_id, drops any qid that errored on either side, and
reports absolute + Δ metrics with paired bootstrap 95% CIs and
Wilcoxon signed-rank p-values. Splits by:
  * gold sp count (2 vs >=3)
  * OFF-side rounds (1 vs 2-3 vs >=4)  — a proxy for "easy vs hard"
    that doesn't leak the ON-side controller decision.
  * max_t bucket (t=1-3 easy vs t=4-8 hard) — where controller gains
    are expected to concentrate.

Multi-seed support:
  * --on accepts K files (predictor seeds 0..K-1). Per-qid metrics
    are averaged across seeds BEFORE bootstrapping. This is the
    correct order — bootstrap resamples over qids, and each qid's
    "true" score under the treatment is the seed-mean.
  * Bootstrap uses the SAME resampled index set for ON and OFF so
    the paired structure survives.
  * Wilcoxon runs on per-qid (seed-averaged) deltas.

Usage — single seed (backward compatible):
    python scripts/paired_analysis.py \\
        --on  storage/eval_v4_on.jsonl \\
        --off storage/eval_v4_off.jsonl

Usage — 3-seed (the warm-raven B4 setup):
    python scripts/paired_analysis.py \\
        --on  storage/eval_v4_on_s0.jsonl \\
              storage/eval_v4_on_s1.jsonl \\
              storage/eval_v4_on_s2.jsonl \\
        --off storage/eval_v4_off.jsonl \\
        --bootstrap-B 1000 --seed 0

Output (per stratum): both absolute mean±SE and paired Δ w/ CI + p.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Iterable


# ---------------------------------------------------------------------- #
# Metric primitives (unchanged from v1 for backward-compat).             #
# ---------------------------------------------------------------------- #
def _norm(s: str) -> str:
    import re, string
    if not s: return ""
    s = s.lower()
    s = re.sub(r"\b(a|an|the)\b", " ", s)
    s = "".join(ch for ch in s if ch not in set(string.punctuation))
    return " ".join(s.split())


def _f1(pred: str, gold: str) -> tuple[float, float, float]:
    p = _norm(pred).split()
    g = _norm(gold).split()
    if not p or not g:
        return 0.0, 0.0, 0.0
    common: dict = {}
    # Increment (not reset to 0) so token counts survive to the g-loop
    # below. Original code shipped with `common[t] = common.get(t, 0)`
    # which zeroed everything and made ans_f1 / sp_f1 / joint_f1 always 0.
    for t in p:
        common[t] = common.get(t, 0) + 1
    common_count = 0
    for t in g:
        if common.get(t, 0) > 0:
            common_count += 1
            common[t] -= 1
    if common_count == 0: return 0.0, 0.0, 0.0
    prec = common_count / len(p)
    rec = common_count / len(g)
    return 2*prec*rec/(prec+rec), prec, rec


def _sp_metrics(pred: list, gold: list) -> tuple[float, float, float, float]:
    cur = set(map(tuple, pred or []))
    au = set(map(tuple, gold or []))
    tp = len(cur & au); fp = len(cur - au); fn = len(au - cur)
    prec = tp/(tp+fp) if tp+fp else 0.0
    rec  = tp/(tp+fn) if tp+fn else 0.0
    f1 = 2*prec*rec/(prec+rec) if prec+rec else 0.0
    em = 1.0 if fp+fn==0 else 0.0
    return em, f1, prec, rec


def _row_scores(row: dict) -> dict:
    pred_ans = row.get("short_ans") or row.get("pred_ans") or ""
    gold_ans = row.get("gold_ans") or ""
    ans_em = 1.0 if _norm(pred_ans) == _norm(gold_ans) and pred_ans else 0.0
    ans_f1, ans_p, ans_r = _f1(pred_ans, gold_ans)
    sp_em, sp_f1, sp_p, sp_r = _sp_metrics(row.get("pred_sp") or [], row.get("gold_sp") or [])
    joint_p = ans_p * sp_p; joint_r = ans_r * sp_r
    joint_f1 = 2*joint_p*joint_r/(joint_p+joint_r) if (joint_p+joint_r) else 0.0
    joint_em = ans_em * sp_em
    return {
        "ans_em": ans_em*100, "ans_f1": ans_f1*100,
        "sp_em": sp_em*100, "sp_f1": sp_f1*100, "sp_prec": sp_p*100, "sp_rec": sp_r*100,
        "joint_em": joint_em*100, "joint_f1": joint_f1*100,
        "rounds": float(row.get("rounds") or 0),
        "n_gold_sp": len(row.get("gold_sp") or []),
    }


_METRIC_KEYS = ("ans_em","ans_f1","sp_em","sp_f1","sp_prec","sp_rec","joint_em","joint_f1","rounds")


# ---------------------------------------------------------------------- #
# Loading                                                                #
# ---------------------------------------------------------------------- #
def _load(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for line in path.open("r", encoding="utf-8"):
        line = line.strip()
        if not line: continue
        r = json.loads(line)
        if r.get("error"): continue
        qid = r.get("query_id")
        if qid: out[qid] = r
    return out


def _load_multi(paths: list[Path]) -> list[dict[str, dict]]:
    """Load K files → list of {qid: row}."""
    return [_load(p) for p in paths]


def _seed_averaged_scores(rows_per_seed: list[dict]) -> dict:
    """Given the SAME qid's row from K seeds, return metric-wise mean
    across seeds. rounds is averaged as a float. All entries must have
    parseable pred_ans / pred_sp; missing seeds are silently dropped."""
    per_metric_vals = defaultdict(list)
    for row in rows_per_seed:
        if row is None: continue
        s = _row_scores(row)
        for m in _METRIC_KEYS:
            per_metric_vals[m].append(s[m])
    if not per_metric_vals:
        return {m: 0.0 for m in _METRIC_KEYS}
    return {m: statistics.mean(per_metric_vals[m]) if per_metric_vals[m] else 0.0
            for m in _METRIC_KEYS}


# ---------------------------------------------------------------------- #
# Stats: bootstrap (paired) + Wilcoxon.                                  #
# ---------------------------------------------------------------------- #
def _bootstrap_paired(
    on_vals: list[float],
    off_vals: list[float],
    B: int,
    seed: int,
) -> dict:
    """Paired bootstrap: resample the SAME qid index set for both on and
    off, so paired structure is preserved. Returns SE + 95% CI for on,
    off, and delta."""
    assert len(on_vals) == len(off_vals)
    n = len(on_vals)
    if n == 0:
        return {"on_mean": 0.0, "on_se": 0.0,
                "off_mean": 0.0, "off_se": 0.0,
                "delta_mean": 0.0, "delta_se": 0.0,
                "delta_ci_lo": 0.0, "delta_ci_hi": 0.0}
    rng = random.Random(seed)
    on_boots: list[float] = []
    off_boots: list[float] = []
    d_boots: list[float] = []
    for _ in range(B):
        idxs = [rng.randrange(n) for _ in range(n)]
        on_sum = off_sum = d_sum = 0.0
        for i in idxs:
            on_sum += on_vals[i]
            off_sum += off_vals[i]
            d_sum += on_vals[i] - off_vals[i]
        on_boots.append(on_sum / n)
        off_boots.append(off_sum / n)
        d_boots.append(d_sum / n)
    d_boots.sort()
    lo_ix = int(0.025 * B); hi_ix = int(0.975 * B)
    return {
        "on_mean": statistics.mean(on_vals),
        "on_se":   statistics.stdev(on_boots) if len(on_boots) > 1 else 0.0,
        "off_mean": statistics.mean(off_vals),
        "off_se":   statistics.stdev(off_boots) if len(off_boots) > 1 else 0.0,
        "delta_mean": statistics.mean([o - f for o, f in zip(on_vals, off_vals)]),
        "delta_se":   statistics.stdev(d_boots) if len(d_boots) > 1 else 0.0,
        "delta_ci_lo": d_boots[lo_ix],
        "delta_ci_hi": d_boots[hi_ix],
    }


def _wilcoxon_p(deltas: list[float]) -> float:
    """Wilcoxon signed-rank two-sided p. Falls back to scipy; if scipy
    is unavailable return NaN. Zero deltas dropped (Wilcoxon
    convention). Not defined for n<6."""
    try:
        from scipy.stats import wilcoxon
    except ImportError:
        return float("nan")
    nz = [d for d in deltas if d != 0.0]
    if len(nz) < 6:
        return float("nan")
    try:
        _, p = wilcoxon(nz, alternative="two-sided", zero_method="wilcox")
        return float(p)
    except Exception:
        return float("nan")


def _p_str(p: float) -> str:
    if math.isnan(p): return "n/a"
    if p < 0.001: return "<.001"
    if p < 0.01:  return f"{p:.3f}"
    if p < 0.05:  return f"{p:.3f}"
    return f"{p:.2f} n.s."


# ---------------------------------------------------------------------- #
# Reporting                                                              #
# ---------------------------------------------------------------------- #
def _report(
    paired: list[tuple[dict, dict]],
    label: str,
    B: int = 1000,
    boot_seed: int = 0,
) -> None:
    """paired: list of (on_scores_dict, off_scores_dict) where each
    dict is already the seed-averaged _row_scores output."""
    n = len(paired)
    print(f"\n=== {label} (n={n}) ===")
    if n == 0:
        print("  (empty)")
        return
    header = f"{'metric':<10} {'ctrl_on':>15} {'ctrl_off':>15} {'Δ (ci)':>25} {'p_wilcox':>10}"
    print(header)
    print("-" * len(header))
    for m in _METRIC_KEYS:
        on_vals = [a[m] for a, _ in paired]
        off_vals = [b[m] for _, b in paired]
        deltas = [o - f for o, f in zip(on_vals, off_vals)]
        boot = _bootstrap_paired(on_vals, off_vals, B=B, seed=boot_seed)
        p = _wilcoxon_p(deltas)
        on_str = f"{boot['on_mean']:>7.2f} ± {boot['on_se']:.2f}"
        off_str = f"{boot['off_mean']:>7.2f} ± {boot['off_se']:.2f}"
        d_str = (f"{boot['delta_mean']:>+6.2f} "
                 f"[{boot['delta_ci_lo']:>+.2f},{boot['delta_ci_hi']:>+.2f}]")
        print(f"{m:<10} {on_str:>15} {off_str:>15} {d_str:>25} {_p_str(p):>10}")


# ---------------------------------------------------------------------- #
# Cohort building                                                        #
# ---------------------------------------------------------------------- #
def _build_paired(
    on_maps: list[dict[str, dict]],
    off_map: dict[str, dict],
    *,
    drop_sp_error: bool = False,
) -> list[tuple[dict, dict, dict]]:
    """Return list of (on_seed_avg_scores, off_scores, off_row) tuples
    for all qids present in ALL sides (off + every on seed).

    Third element `off_row` is kept so downstream stratifiers can read
    off-side fields (rounds, gold_sp) without re-parsing.
    """
    common = set(off_map.keys())
    for on_map in on_maps:
        common &= set(on_map.keys())

    paired: list[tuple[dict, dict, dict]] = []
    n_sp_dropped = 0
    for qid in common:
        off_row = off_map[qid]
        on_rows = [m[qid] for m in on_maps]
        if drop_sp_error:
            if off_row.get("sp_error"):
                n_sp_dropped += 1; continue
            if any(r.get("sp_error") for r in on_rows):
                n_sp_dropped += 1; continue
        on_scores = _seed_averaged_scores(on_rows)
        off_scores = _row_scores(off_row)
        paired.append((on_scores, off_scores, off_row))
    if drop_sp_error and n_sp_dropped:
        print(f"dropped for sp_err: {n_sp_dropped}")
    return paired


def _strip(paired3):
    """Drop the third element (off_row) for reporting."""
    return [(a, b) for a, b, _ in paired3]


# ---------------------------------------------------------------------- #
# CLI                                                                    #
# ---------------------------------------------------------------------- #
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--on", required=True, type=Path, nargs="+",
                   help="One or more ctrl_on eval jsonls (one per predictor seed). "
                        "Per-qid metrics are averaged across seeds before "
                        "bootstrap/Wilcoxon so paired structure is preserved.")
    p.add_argument("--off", required=True, type=Path,
                   help="Single ctrl_off eval jsonl (same seed=1 sample). "
                        "Shared across all --on files.")
    p.add_argument("--drop-sp-error", action="store_true",
                   help="Drop qids where sp_err was set on either side "
                        "(picker failure produced Sp=0 fallback).")
    p.add_argument("--bootstrap-B", type=int, default=1000,
                   help="Number of bootstrap resamples (default 1000).")
    p.add_argument("--seed", type=int, default=0,
                   help="Bootstrap RNG seed.")
    args = p.parse_args()

    on_maps = _load_multi(args.on)
    off_map = _load(args.off)
    print(f"loaded: {len(args.on)} ctrl_on file(s), sizes={[len(m) for m in on_maps]}")
    print(f"        1 ctrl_off file, size={len(off_map)}")

    paired3 = _build_paired(on_maps, off_map, drop_sp_error=args.drop_sp_error)
    paired = _strip(paired3)
    print(f"paired qids: {len(paired)}")

    # --- Overall ---
    _report(paired, "ALL PAIRED", B=args.bootstrap_B, boot_seed=args.seed)

    # --- Stratum: gold sp count ---
    by_nsp: dict[str, list] = defaultdict(list)
    for a, b, off_row in paired3:
        n_sp = len(off_row.get("gold_sp") or [])
        bucket = "n_sp==2" if n_sp == 2 else ("n_sp>=3" if n_sp >= 3 else "n_sp<2")
        by_nsp[bucket].append((a, b))
    for k in sorted(by_nsp):
        _report(by_nsp[k], f"stratum {k}", B=args.bootstrap_B, boot_seed=args.seed)

    # --- Stratum: OFF rounds (proxy for hardness) ---
    # 3-way split to expose where controller helps most:
    #   r=1     : one-shot easy — controller should NOT hurt these
    #   r=2-3   : mid            — small gains expected
    #   r=4+    : long tail       — biggest gains (t=8 posterior-collapsed batch lives here)
    by_r: dict[str, list] = defaultdict(list)
    for a, b, off_row in paired3:
        rb = int(off_row.get("rounds") or 0)
        if rb <= 1:      bucket = "OFF r==1"
        elif rb <= 3:    bucket = "OFF r=2-3"
        else:            bucket = "OFF r>=4"
        by_r[bucket].append((a, b))
    for k in sorted(by_r):
        _report(by_r[k], f"stratum {k}", B=args.bootstrap_B, boot_seed=args.seed)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
