"""Strategy matrix over the candidate cache from generate_candidates.py.

Evaluates five selection strategies on the generate_candidates output:

  (a) no_rag           — the baseline (no-reference) generation from
                          runs_gate/<split>_baseline.jsonl
  (b) reranker_top1    — the candidate whose cand_rank == 1 (from BGE-M3
                          hybrid + reranker)
  (c) entropy_top1     — pick the candidate with the largest gating score.
                          Score = alpha * dH_full + (1 - alpha) * dH_sel
                          (default alpha=0.5). --alpha sweeps it.
  (d) oracle           — pick the candidate whose evaluator score is highest
                          (upper bound). Metric configurable via --oracle-key.
  (e) random           — deterministic random pick (seed configurable).

Each strategy's chosen ruby text is scored with evaluator.score_pair against
the gold .rb. We aggregate per-metric means across queries and emit a compact
summary table.

Usage:
    # standard eval on test split
    python token_conf/dsl_gate/strategy_matrix.py --split test
    # sweep alpha for the entropy strategy (ablation)
    python token_conf/dsl_gate/strategy_matrix.py --split dev --alpha 0 0.3 0.5 0.7 1.0
    # write per-qid rows for downstream analysis / calibration plots
    python token_conf/dsl_gate/strategy_matrix.py --split test --per-qid-out runs_gate/test_perqid.jsonl
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

from evaluator import score_pair   # noqa: E402


RUNS_GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD_RB_DIR = ROOT / "data" / "corpus" / "generated_dsl"

METRIC_KEYS = ("dataset_em", "decay_card_em", "sel_method_f1",
               "sel_slot_f1", "struct_mean",
               "strict_method_f1", "strict_slot_f1", "strict_struct_mean")


# ------------- loaders --------------------------------------------------


def _load_jsonl(p: Path) -> list[dict]:
    out: list[dict] = []
    if not p.exists():
        return out
    with p.open() as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            out.append(json.loads(ln))
    return out


def _index_by_qid(rows: list[dict], keep_last: bool = True) -> dict[str, list[dict]]:
    """Group candidate rows by qid. keep_last=True dedupes to the most recent
    entry per (qid, cand_stem) — protects against interrupted --force reruns
    that appended twice."""
    by_qid: dict[str, dict[str, dict]] = defaultdict(dict)
    for r in rows:
        by_qid[r["qid"]][r["cand_stem"]] = r
    return {q: list(cands.values()) for q, cands in by_qid.items()}


def _load_baseline(rows: list[dict]) -> dict[str, dict]:
    """qid -> baseline row (latest wins if duplicated)."""
    out: dict[str, dict] = {}
    for r in rows:
        out[r["qid"]] = r
    return out


def _gold_text(qid: str) -> str | None:
    p = GOLD_RB_DIR / f"{qid}.rb"
    if not p.exists():
        return None
    return p.read_text(encoding="utf-8")


# ------------- strategies -----------------------------------------------


def pick_no_rag(qid: str, cands: list[dict], base: dict) -> tuple[str, str]:
    return "baseline", base.get("gen_ruby_first", "") or ""


def pick_reranker_top1(qid: str, cands: list[dict], base: dict) -> tuple[str, str]:
    real = [c for c in cands if not c.get("is_gold")]
    r = min(real or cands, key=lambda c: c.get("cand_rank", 999))
    return r["cand_stem"], r.get("gen_ruby_first", "") or ""


def pick_oracle_gold(qid: str, cands: list[dict], base: dict) -> tuple[str, str] | None:
    """Pick the synthetic gold-DSL candidate (cand_rank=0 / is_gold=True) if
    generate_candidates was run with --include-gold; otherwise the strategy is skipped."""
    for c in cands:
        if c.get("is_gold") or c.get("cand_rank") == 0:
            return c["cand_stem"], c.get("gen_ruby_first", "") or ""
    return None


def make_pick_entropy(alpha: float, field_a: str = "dH_full",
                      field_b: str = "dH_sel") -> Callable:
    """Pick the candidate maximizing alpha*field_a + (1-alpha)*field_b.

    Defaults reproduce the original token-ΔH gate. --gate-pair swaps in the
    two-component signals (e.g. dH_stat_code,dH14) without touching this math.
    """
    def _pick(qid: str, cands: list[dict], base: dict) -> tuple[str, str]:
        real = [c for c in cands if not c.get("is_gold")]
        def _score(c: dict) -> float:
            f = c.get(field_a)
            s = c.get(field_b)
            f = 0.0 if f is None else f
            s = 0.0 if s is None else s
            return alpha * f + (1 - alpha) * s
        r = max(real or cands, key=_score)
        return r["cand_stem"], r.get("gen_ruby_first", "") or ""
    return _pick


def make_pick_oracle(oracle_key: str) -> Callable:
    def _pick(qid: str, cands: list[dict], base: dict) -> tuple[str, str]:
        real = [c for c in cands if not c.get("is_gold")]
        gold = _gold_text(qid) or ""
        best_r, best_s = None, -1.0
        for c in (real or cands):
            ruby = c.get("gen_ruby_first", "") or ""
            m = score_pair(ruby, gold)
            s = m[oracle_key]
            if s > best_s:
                best_s, best_r = s, c
        return best_r["cand_stem"], best_r.get("gen_ruby_first", "") or ""
    return _pick


def make_pick_random(seed: int) -> Callable:
    rng = random.Random(seed)
    def _pick(qid: str, cands: list[dict], base: dict) -> tuple[str, str]:
        real = [c for c in cands if not c.get("is_gold")]
        r = rng.choice(real or cands)
        return r["cand_stem"], r.get("gen_ruby_first", "") or ""
    return _pick


# ------------- evaluation -----------------------------------------------


def _score_one(pred_ruby: str, gold_text: str) -> dict:
    if not pred_ruby or not gold_text:
        return {k: 0.0 for k in METRIC_KEYS}
    m = score_pair(pred_ruby, gold_text)
    return {k: float(m[k]) for k in METRIC_KEYS}


def evaluate_strategy(
    name: str,
    pick: Callable,
    qids: list[str],
    by_qid: dict[str, list[dict]],
    baselines: dict[str, dict],
) -> tuple[dict[str, float], list[dict]]:
    per_qid: list[dict] = []
    agg: dict[str, list[float]] = {k: [] for k in METRIC_KEYS}
    n_missing = 0
    for qid in qids:
        cands = by_qid.get(qid) or []
        base = baselines.get(qid) or {}
        if not base and not cands:
            n_missing += 1
            continue
        picked = pick(qid, cands, base)
        if picked is None:
            n_missing += 1
            continue
        cand_stem, ruby = picked
        gold = _gold_text(qid) or ""
        m = _score_one(ruby, gold)
        row = {"qid": qid, "strategy": name, "cand_stem": cand_stem, **m}
        per_qid.append(row)
        for k in METRIC_KEYS:
            agg[k].append(m[k])
    means = {k: (statistics.fmean(v) if v else 0.0) for k, v in agg.items()}
    means["n_eval"] = float(len(per_qid))
    means["n_missing"] = float(n_missing)
    return means, per_qid


# ------------- pretty-print ---------------------------------------------


def _print_table(rows: list[tuple[str, dict]]) -> None:
    cols = list(METRIC_KEYS) + ["n_eval"]
    widths = {c: max(len(c), 8) for c in cols}
    header = f"{'strategy':<22}  " + "  ".join(f"{c:>{widths[c]}}" for c in cols)
    print(header)
    print("-" * len(header))
    for name, m in rows:
        cells = []
        for c in cols:
            v = m.get(c, 0.0)
            if c == "n_eval":
                cells.append(f"{int(v):>{widths[c]}d}")
            else:
                cells.append(f"{v:>{widths[c]}.4f}")
        print(f"{name:<22}  " + "  ".join(cells))


# ------------- main -----------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", choices=["dev", "test", "train", "all"], default="test")
    ap.add_argument("--alpha", type=float, nargs="+", default=[0.5],
                    help="alpha values for entropy_top1 (dH_full weight); "
                         "multiple values sweep and print each as its own row")
    ap.add_argument("--oracle-key", default="struct_mean",
                    choices=list(METRIC_KEYS),
                    help="metric maximized by the oracle picker")
    ap.add_argument("--gate-pair", default="dH_full,dH_sel",
                    help="two candidate-row fields combined by --alpha for the "
                         "entropy gate, e.g. 'dH_stat_code,dH14'. Default "
                         "reproduces the original token-ΔH gate exactly.")
    ap.add_argument("--cand-file", type=Path, default=None,
                    help="override the candidate jsonl (default runs_gate/<split>.jsonl); "
                         "use the two-component *_h2.jsonl files to gate on the two-component H")
    ap.add_argument("--baseline-file", type=Path, default=None,
                    help="override the baseline jsonl (default runs_gate/<split>_baseline.jsonl); "
                         "used together with --cand-file when a T-sweep run writes to a separate tag")
    ap.add_argument("--random-seed", type=int, default=42)
    ap.add_argument("--per-qid-out", type=Path, default=None,
                    help="write per-(qid, strategy) rows here (jsonl)")
    ap.add_argument("--summary-out", type=Path, default=None,
                    help="write summary rows here (jsonl)")
    args = ap.parse_args()

    cand_path = args.cand_file or (RUNS_GATE / f"{args.split}.jsonl")
    cand_rows = _load_jsonl(cand_path)
    base_path = args.baseline_file or (RUNS_GATE / f"{args.split}_baseline.jsonl")
    base_rows = _load_jsonl(base_path)
    if not cand_rows:
        print(f"[!] no candidate rows in {cand_path} — run generate_candidates first")
        return 2
    gate_a, _, gate_b = args.gate_pair.partition(",")
    if not gate_a or not gate_b:
        print(f"[!] --gate-pair needs two comma-separated fields, got {args.gate_pair!r}")
        return 2
    by_qid = _index_by_qid(cand_rows)
    baselines = _load_baseline(base_rows)
    qids = sorted(set(by_qid) | set(baselines))
    print(f"[strategy_matrix] split={args.split} qids={len(qids)} "
          f"({len(baselines)} with baseline, {len(by_qid)} with candidates)")

    strategies: list[tuple[str, Callable]] = []
    strategies.append(("no_rag", pick_no_rag))
    strategies.append(("reranker_top1", pick_reranker_top1))
    gate_label = ("entropy_top1" if (gate_a, gate_b) == ("dH_full", "dH_sel")
                  else f"gate[{gate_a}|{gate_b}]")
    for a in args.alpha:
        strategies.append((f"{gate_label}(a={a:g})",
                           make_pick_entropy(a, gate_a, gate_b)))
    strategies.append((f"oracle({args.oracle_key})", make_pick_oracle(args.oracle_key)))
    strategies.append(("oracle_gold", pick_oracle_gold))
    strategies.append(("random", make_pick_random(args.random_seed)))

    summary: list[tuple[str, dict]] = []
    per_qid_rows: list[dict] = []
    for name, pick in strategies:
        means, rows = evaluate_strategy(name, pick, qids, by_qid, baselines)
        summary.append((name, means))
        per_qid_rows.extend(rows)

    print()
    _print_table(summary)

    if args.per_qid_out:
        args.per_qid_out.parent.mkdir(parents=True, exist_ok=True)
        with args.per_qid_out.open("w") as f:
            for r in per_qid_rows:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        print(f"\n[+] per-qid rows -> {args.per_qid_out}")
    if args.summary_out:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_out.open("w") as f:
            for name, m in summary:
                f.write(json.dumps({"strategy": name, **m}, ensure_ascii=False) + "\n")
        print(f"[+] summary -> {args.summary_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
