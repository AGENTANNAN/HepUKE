"""Sweet-spot gating ablation (offline, no LLM calls).

Motivated by the ∩-shape in fig_calibration.png: argmax(dH_sel) can select
a candidate whose ΔH_sel is *too large* (~0.04+), where the LLM imitates the
reference too closely and downstream match collapses.

We evaluate several "sweet-spot" gating rules on the same generate_candidates outputs used by
strategy_matrix, comparing struct_mean and oracle-recovery to plain entropy_top1(α=0.5)
and the reranker top-1 baseline.

Rules:
  clip@0.02        : discard cand if dH_sel > 0.02 (fallback: rerank rank-1)
  clip@0.02_max    : same, but among remaining candidates pick max(dH_sel)
  soft_penalty     : score = dH_sel if dH_sel <= 0.02 else 0.02 - k*(dH_sel-0.02)
                     (k=1, i.e. reflection back)
  interval_prefer  : score = 1 if 0.005 <= dH_sel < 0.02 else 0.5 if dH_sel > 0
                     else 0
  reranker_then_ok : take reranker top-1, but if its dH_sel > 0.02 fall back to
                     rank-2 (or best rank-k with dH_sel <= 0.02)

Baselines echoed: reranker_top1, entropy_top1(a=0.5), oracle(struct_mean).

Reads:  runs_gate/all.jsonl, runs_gate/all_baseline.jsonl, batch/runs/<qid>.rb
Writes: results/sweet_spot_summary.md
"""
from __future__ import annotations

import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _prereq import require  # noqa: E402

from evaluator import score_pair               # noqa: E402
from pipeline_common import strip_code_fence      # noqa: E402


GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"

CAND_ROWS = [json.loads(l) for l in require(GATE / "all.jsonl").open()]
BASE_ROWS = [json.loads(l) for l in require(GATE / "all_baseline.jsonl").open()]

# ---- index -----
by_qid: dict[str, list[dict]] = defaultdict(list)
for r in CAND_ROWS:
    by_qid[r["qid"]].append(r)
baselines: dict[str, dict] = {b["qid"]: b for b in BASE_ROWS}


def _score(pred_ruby: str, gold_text: str) -> float:
    if not pred_ruby or not gold_text:
        return 0.0
    return float(score_pair(pred_ruby, gold_text)["struct_mean"])


def _gold(qid: str) -> str:
    p = GOLD / f"{qid}.rb"
    return p.read_text() if p.exists() else ""


# ---- pickers -----------


def pick_no_rag(qid, cands, base):
    return "baseline", base.get("gen_ruby_first", "") or ""


def pick_reranker(qid, cands, base):
    r = min(cands, key=lambda c: c.get("cand_rank", 999))
    return r["cand_stem"], r.get("gen_ruby_first", "") or ""


def pick_entropy_top1(alpha: float):
    def _pick(qid, cands, base):
        def s(c):
            f, sv = c.get("dH_full") or 0.0, c.get("dH_sel") or 0.0
            return alpha * f + (1 - alpha) * sv
        r = max(cands, key=s)
        return r["cand_stem"], r.get("gen_ruby_first", "") or ""
    return _pick


def pick_clip_max(thr: float = 0.02):
    """From cands with dH_sel <= thr pick the one with largest dH_sel; else
    reranker top-1."""
    def _pick(qid, cands, base):
        ok = [c for c in cands if (c.get("dH_sel") or -1e9) <= thr]
        if not ok:
            return pick_reranker(qid, cands, base)
        r = max(ok, key=lambda c: c.get("dH_sel") or -1e9)
        return r["cand_stem"], r.get("gen_ruby_first", "") or ""
    return _pick


def pick_soft_penalty(sweet: float = 0.02, k: float = 1.0):
    """score = dH_sel if dH_sel <= sweet else sweet - k*(dH_sel-sweet)."""
    def _pick(qid, cands, base):
        def s(c):
            d = c.get("dH_sel") or 0.0
            return d if d <= sweet else sweet - k * (d - sweet)
        r = max(cands, key=s)
        return r["cand_stem"], r.get("gen_ruby_first", "") or ""
    return _pick


def pick_interval_prefer(lo: float = 0.005, hi: float = 0.02):
    """Prefer candidates whose dH_sel is inside [lo, hi]; break ties by dH_sel.
    If none inside, prefer positive dH_sel; else reranker top-1."""
    def _pick(qid, cands, base):
        inside = [c for c in cands if lo <= (c.get("dH_sel") or -1e9) < hi]
        if inside:
            r = max(inside, key=lambda c: c.get("dH_sel") or -1e9)
            return r["cand_stem"], r.get("gen_ruby_first", "") or ""
        pos = [c for c in cands if (c.get("dH_sel") or -1e9) > 0]
        if pos:
            r = max(pos, key=lambda c: c.get("dH_sel") or -1e9)
            return r["cand_stem"], r.get("gen_ruby_first", "") or ""
        return pick_reranker(qid, cands, base)
    return _pick


def pick_reranker_then_ok(thr: float = 0.02):
    """Reranker top-1 if its dH_sel <= thr, else next-best cand with dH_sel<=thr."""
    def _pick(qid, cands, base):
        cs = sorted(cands, key=lambda c: c.get("cand_rank", 999))
        for c in cs:
            if (c.get("dH_sel") or -1e9) <= thr:
                return c["cand_stem"], c.get("gen_ruby_first", "") or ""
        return cs[0]["cand_stem"], cs[0].get("gen_ruby_first", "") or ""
    return _pick


def pick_oracle():
    """cheat: choose cand with highest downstream struct_mean."""
    def _pick(qid, cands, base):
        gold = _gold(qid)
        best_c, best_s = None, -1.0
        for c in cands:
            s = _score(strip_code_fence(c.get("gen_ruby_first", "")), gold)
            if s > best_s:
                best_s, best_c = s, c
        return best_c["cand_stem"], best_c.get("gen_ruby_first", "") or ""
    return _pick


def pick_random(seed: int = 42):
    import random
    rng = random.Random(seed)
    def _pick(qid, cands, base):
        c = rng.choice(cands)
        return c["cand_stem"], c.get("gen_ruby_first", "") or ""
    return _pick


# ---- eval loop -----------


def evaluate(name, pick):
    scores = []
    for qid in sorted(by_qid):
        cands = by_qid[qid]
        base = baselines.get(qid, {})
        _, ruby = pick(qid, cands, base)
        gold = _gold(qid)
        scores.append(_score(strip_code_fence(ruby), gold))
    return {
        "name": name,
        "mean_struct": st.mean(scores),
        "n": len(scores),
    }


def main():
    strategies = [
        ("no_rag", pick_no_rag),
        ("reranker_top1", pick_reranker),
        ("random", pick_random()),
        ("entropy_top1(a=0.5)", pick_entropy_top1(0.5)),
        ("entropy_top1(a=1)", pick_entropy_top1(1.0)),
        ("clip@0.02_max", pick_clip_max(0.02)),
        ("clip@0.015_max", pick_clip_max(0.015)),
        ("clip@0.010_max", pick_clip_max(0.010)),
        ("clip@0.008_max", pick_clip_max(0.008)),
        ("soft_penalty(0.02,k=1)", pick_soft_penalty(0.02, 1.0)),
        ("soft_penalty(0.01,k=2)", pick_soft_penalty(0.01, 2.0)),
        ("interval_prefer[0.005,0.02)", pick_interval_prefer(0.005, 0.02)),
        ("interval_prefer[0.000,0.015)", pick_interval_prefer(0.0, 0.015)),
        ("reranker_then_clip@0.02", pick_reranker_then_ok(0.02)),
        ("oracle(struct_mean)", pick_oracle()),
    ]

    rows = []
    for name, pick in strategies:
        rows.append(evaluate(name, pick))

    no_rag = next(r["mean_struct"] for r in rows if r["name"] == "no_rag")
    oracle = next(r["mean_struct"] for r in rows if r["name"] == "oracle(struct_mean)")
    rr = next(r["mean_struct"] for r in rows if r["name"] == "reranker_top1")
    gap = oracle - no_rag

    out = RESULTS / "sweet_spot_summary.md"
    lines = ["# Sweet-spot gating ablation (all n=267)", ""]
    lines.append(f"no_rag={no_rag:.4f}  reranker_top1={rr:.4f}  oracle={oracle:.4f}  gap={gap:.4f}")
    lines.append("")
    lines.append("| strategy | struct_mean | Δ vs reranker (pp) | oracle-recovery |")
    lines.append("|---|---|---|---|")
    for r in rows:
        v = r["mean_struct"]
        rec = (v - no_rag) / gap * 100 if gap > 0 else 0
        dvr = (v - rr) * 100
        lines.append(f"| {r['name']} | {v:.4f} | {dvr:+.2f} | {rec:+.1f}% |")
    out.write_text("\n".join(lines))
    print(out.read_text())
    print(f"\n[+] saved {out}")


if __name__ == "__main__":
    main()
