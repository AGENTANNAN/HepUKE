"""STEP 3 — Compute per-token entropy GOLD for every run, code-tokens-only.

Reuses the verified logic in compute_entropy.py (char-level Ruby lexer that drops
comment + whitespace tokens, tail-bucket entropy). Writes batch/gold/<stem>.json per
paper and a batch/gold_index.csv summary (mean/sum/max entropy, token counts) so the
whole corpus is comparable at a glance. Pure local — no API calls.

    python token_conf/gold_batch.py                 # all runs
    python token_conf/gold_batch.py --only 1310.4101v3
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pipeline_common as pc
from compute_entropy import build_gold


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--keep-all", action="store_true")
    args = ap.parse_args()

    pc.GOLD_DIR.mkdir(parents=True, exist_ok=True)
    runs = ([pc.RUN_DIR / f"{s}.json" for s in args.only] if args.only
            else sorted(pc.RUN_DIR.glob("*.json")))

    rows = []
    ok = fail = 0
    for rj in runs:
        stem = rj.stem
        try:
            d = json.loads(rj.read_text())
            toks = d.get("content_logprobs") or []
            if not toks:
                raise ValueError("no content_logprobs")
            gold = build_gold(toks, {"model": d.get("model")}, rj, keep_all=args.keep_all)
            (pc.GOLD_DIR / f"{stem}.json").write_text(json.dumps(gold, ensure_ascii=False))
            s = gold["summary"]
            rows.append({
                "stem": stem,
                "n_total": s["n_tokens_total"], "n_code": s["n_tokens_code"],
                "n_comment": s["n_tokens_comment"], "n_ws": s["n_tokens_ws"],
                "n_gold": s["n_tokens_gold"],
                "H_mean_nats": round(s["H_tail_nats_mean"], 5),
                "H_sum_nats": round(s["H_tail_nats_sum"], 4),
                "H_max_nats": round(s["H_tail_nats_max"], 4),
                "H_mean_bits": round(s["H_tail_bits_mean"], 5),
                "mean_tail_mass": s["mean_tail_mass"],
            })
            ok += 1
        except Exception as e:
            fail += 1
            print(f"ERR {stem}: {e}")

    rows.sort(key=lambda r: r["H_mean_nats"], reverse=True)
    idx = pc.BATCH / "gold_index.csv"
    if rows:
        with idx.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    print(f"[step3] gold written for {ok} papers (fail={fail}) -> {pc.GOLD_DIR}")
    print(f"[step3] index -> {idx}")
    if rows:
        import statistics as st
        hs = [r["H_mean_nats"] for r in rows]
        print(f"[step3] corpus mean-entropy: min={min(hs):.4f} "
              f"median={st.median(hs):.4f} max={max(hs):.4f} nats")
        print("[step3] highest-entropy papers (candidate hard/ambiguous specs):")
        for r in rows[:5]:
            print(f"    {r['stem']:<16} H_mean={r['H_mean_nats']:.4f} nats  "
                  f"({r['n_gold']} code tok)")


if __name__ == "__main__":
    main()
