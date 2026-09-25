"""Interval-gated confidence for knowledge-base self-evolution.

kb_evolution used a ONE-SIDED confidence gate `c >= tau_write` (reject only low
confidence). §6.4's ∩-shape calibration found the downstream metric drops at the
HIGH end of dH_sel too (pathological right tail: a reference that drastically
rewrites the selection is often a bad precedent), and kb_evolution found dH_full is
INVERTED (defective generations get HIGHER blended confidence). Both point the
same way: a one-sided gate admits high-confidence garbage.

interval_gate replaces the gate with an INTERVAL  `tau_low <= c <= tau_high`, rejecting
BOTH the low-confidence (uncertain) and high-confidence (pathological) tails,
and sweeps (alpha, tau_high) to measure whether the interval gate cuts the
false-positive rate (admitted entries with struct_mean < 0.5) without gutting
admission.

Reuses kb_evolution's validation / dedup / confidence machinery; no LLM calls (offline).

Usage:
    python token_conf/dsl_gate/interval_gate.py
"""
from __future__ import annotations

import json
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _prereq import require  # noqa: E402

from evaluator import score_pair                     # noqa: E402
from pipeline_common import strip_code_fence            # noqa: E402
from splits import load_splits                          # noqa: E402
from kb_evolution import (                           # noqa: E402
    validation_gate, ast_hash, desc_similarity, confidence_from_dH,
)

GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"
DESC = ROOT / "data" / "corpus" / "descriptions"
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results" / "kb_evolution"


def _entropy_top1(cands: list[dict], alpha: float) -> dict:
    def score(c):
        f = c.get("dH_full") or 0.0
        s = c.get("dH_sel") or 0.0
        return alpha * f + (1 - alpha) * s
    return max(cands, key=score)


def simulate(tau_low, tau_high, alpha, seed=42):
    splits = load_splits()
    train_stems = set(splits["train"])
    query_stems = list(splits["dev"]) + list(splits["test"])
    import random
    rng = random.Random(seed)
    rng.shuffle(query_stems)

    kb_entries: list[tuple[str, str]] = []   # (ast_hash, desc)
    for stem in sorted(train_stems):
        rb = (GOLD / f"{stem}.rb"); dc = (DESC / f"{stem}.txt")
        if not rb.exists() or not dc.exists():
            continue
        kb_entries.append((ast_hash(rb.read_text()), dc.read_text().strip()))

    cand = [json.loads(l) for l in require(GATE / "all.jsonl").open()]
    by_qid: dict[str, list[dict]] = defaultdict(list)
    for r in cand:
        by_qid[r["qid"]].append(r)

    # downstream struct per (qid, cand_stem), on demand
    def struct_of(qid, cs, ruby):
        g = (GOLD / f"{qid}.rb")
        if not g.exists():
            return None
        return float(score_pair(strip_code_fence(ruby), g.read_text())["struct_mean"])

    stats = Counter()
    admitted = []          # (qid, confidence, struct)
    rejected_high = []     # (qid, confidence, struct)
    rejected_low = []      # (qid, confidence, struct)

    for qid in query_stems:
        cs = by_qid.get(qid)
        if not cs:
            stats["missing"] += 1
            continue
        top1 = _entropy_top1(cs, alpha)
        ruby = strip_code_fence(top1.get("gen_ruby_first", "") or "")
        c = confidence_from_dH(top1.get("dH_full"), top1.get("dH_sel"), alpha)
        desc = (DESC / f"{qid}.txt").read_text().strip()
        sm = struct_of(qid, top1.get("cand_stem"), ruby)

        # interval confidence gate
        if c < tau_low:
            stats["reject_low"] += 1
            rejected_low.append((qid, c, sm))
            continue
        if c > tau_high:
            stats["reject_high"] += 1
            rejected_high.append((qid, c, sm))
            continue
        # validation gate
        v_ok, v_reason = validation_gate(ruby)
        if not v_ok:
            stats[f"reject_val:{v_reason.split(':',1)[0]}"] += 1
            continue
        # dedup gate: ast_hash match AND desc similarity > theta_dup
        h = ast_hash(ruby)
        dup = False
        for e_hash, e_desc in kb_entries:
            if e_hash == h and desc_similarity(desc, e_desc) > 0.35:
                dup = True
                break
        if dup:
            stats["reject_dedup"] += 1
            continue
        kb_entries.append((h, desc))
        stats["admit"] += 1
        admitted.append((qid, c, sm))

    def _m(vals, field):
        xs = [v[field] for v in vals if v[field] is not None]
        return st.mean(xs) if xs else float("nan")

    fp = [(q, c, s) for q, c, s in admitted if s is not None and s < 0.5]
    return {
        "admit": stats["admit"],
        "reject_low": stats["reject_low"],
        "reject_high": stats["reject_high"],
        "reject_val": sum(v for k, v in stats.items() if k.startswith("reject_val:")),
        "reject_dedup": stats["reject_dedup"],
        "admitted_struct": _m(admitted, 2),
        "rejected_high_struct": _m(rejected_high, 2),
        "rejected_low_struct": _m(rejected_low, 2),
        "n_fp": len(fp),
        "fp": fp,
    }


def main():
    rows = []
    for alpha in [0.5, 0.0]:
        for tau_high in [1.0, 0.95, 0.9, 0.85, 0.8]:
            r = simulate(tau_low=0.5, tau_high=tau_high, alpha=alpha)
            r["alpha"] = alpha; r["tau_high"] = tau_high
            rows.append(r)

    print(f"\n{'α':>4} {'τ_high':>7} {'admit':>5} {'rejL':>4} {'rejH':>4} {'rejVal':>6} "
          f"{'rejDup':>6} {'admStruct':>9} {'rejHStruct':>10} {'FP<0.5':>6}")
    print("-" * 72)
    for r in rows:
        print(f"{r['alpha']:>4.1f} {r['tau_high']:>7.2f} {r['admit']:>5} "
              f"{r['reject_low']:>4} {r['reject_high']:>4} {r['reject_val']:>6} "
              f"{r['reject_dedup']:>6} {r['admitted_struct']:>9.4f} "
              f"{r['rejected_high_struct']:>10.4f} {r['n_fp']:>6}")

    # which false positives does the interval gate remove?
    print("\n--- false positives (struct<0.5) by config ---")
    for r in rows:
        if r["n_fp"]:
            fps = ", ".join(f"{q.split('v')[0]}v…(c={c:.3f})" for q, c, s in r["fp"])
            print(f"  α={r['alpha']:.1f} τ_high={r['tau_high']:.2f}: {r['n_fp']}  [{fps}]")

    md = ["# interval_gate — interval-gated confidence for KB write-back", "",
          "Gate = `tau_low <= c <= tau_high` (reject both tails), tau_low=0.5 fixed.",
          "α blends dH_full/dH_sel in the confidence sigmoid. FP = admitted entries with",
          "struct_mean < 0.5 (false positives the gate let through).",
          "",
          "| α | τ_high | admit | reject_low | reject_high | reject_val | reject_dedup | admitted struct | rejected-high struct | FP<0.5 |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['alpha']:.1f} | {r['tau_high']:.2f} | {r['admit']} | "
                  f"{r['reject_low']} | {r['reject_high']} | {r['reject_val']} | "
                  f"{r['reject_dedup']} | {r['admitted_struct']:.4f} | "
                  f"{r['rejected_high_struct']:.4f} | {r['n_fp']} |")
    out = RESULTS / "interval_gate_sweep.md"
    out.write_text("\n".join(md) + "\n")
    print(f"\n[+] {out}")


if __name__ == "__main__":
    main()
