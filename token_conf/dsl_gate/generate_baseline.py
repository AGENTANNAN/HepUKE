"""Baseline (No-RAG) generation under a chosen skill mode.

Complements generate_candidates.py: reads the same descriptions, runs ONLY the
baseline call, and writes to <tag>_baseline.jsonl and <tag>_raw.jsonl.gz with
kind="baseline". Use when you want the baseline to run under a leaner skill
context (typically `no-example`) than the candidate calls did — so no_rag
cannot secretly borrow the shape of the three worked-examples embedded in the
default skill package.

Existing candidate rows (with kind="cand") are left untouched; the strategy_matrix loader
merges the two sources by kind.

Usage:
    python token_conf/dsl_gate/generate_baseline.py \\
        --split dev --limit 20 \\
        --n-samples 5 --temperature 0.1 \\
        --desc-tier bare --skill-mode no-example \\
        --out-tag dev20_bare_gold_T01_noexample
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import json
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "token_conf"))

import pipeline_common as pc                                    # noqa: E402
from token_conf.dsl_gate.generate_candidates import (            # noqa: E402
    DSL_MODEL, MAX_TOKENS, REQUEST_TIMEOUT, DEFAULT_N,
    EXTRA_PAYLOAD, SYSTEM_TMPL, USER_BASELINE_TMPL,
    _summarize, _mean,
)
from token_conf.dsl_gate.splits import load_splits              # noqa: E402


OUT_DIR = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
_print_lock = threading.Lock()


def _call_once(system, user, base, key, temperature):
    for _ in range(2):
        res = pc.chat(
            system=system, user=user,
            model=DSL_MODEL, base_url=base, api_key=key,
            max_tokens=MAX_TOKENS, temperature=temperature,
            logprobs=True, top_logprobs=20,
            timeout=REQUEST_TIMEOUT, extra_payload=EXTRA_PAYLOAD,
        )
        if res.get("content_logprobs"):
            return res
    return res


def run_qid(qid, *, system, base, key, n_samples, temperature, desc_tier):
    dd = ROOT / "data" / "corpus" / "descriptions"
    full_desc = (dd / f"{qid}.txt").read_text(encoding="utf-8").strip()
    prompt_desc = pc.bare_desc(full_desc) if desc_tier == "bare" else full_desc

    t0 = time.time()
    user = USER_BASELINE_TMPL.format(desc=prompt_desc)
    calls = [_call_once(system, user, base, key, temperature) for _ in range(n_samples)]
    summaries = [_summarize(c["content_logprobs"]) for c in calls]

    raw = []
    for si, c in enumerate(calls):
        raw.append({
            "qid": qid, "kind": "baseline", "cand_stem": None,
            "cand_rank": None, "sample_idx": si,
            "answer": c.get("answer", ""),
            "content_logprobs": c.get("content_logprobs") or [],
            "usage": c.get("usage"),
        })

    row = {
        "qid": qid,
        "H_full_base_mean": _mean([s["H_full"] for s in summaries]),
        "H_sel_base_mean": _mean([s["H_sel"] for s in summaries]),
        "n_tok_base_mean": _mean([float(s["n_tok"]) for s in summaries]),
        "n_samples": n_samples,
        "gen_ruby_first": pc.strip_code_fence(calls[0]["answer"]),
        "wall_ms_base": (time.time() - t0) * 1000,
    }
    return row, raw


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="dev")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--n-samples", type=int, default=DEFAULT_N)
    ap.add_argument("--temperature", type=float, default=0.7)
    ap.add_argument("--desc-tier", choices=["full", "bare"], default="full")
    ap.add_argument("--skill-mode", choices=list(pc.SKILL_MODES), default="no-example",
                    help="skill package mode for THIS baseline run (default no-example)")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--out-tag", required=True,
                    help="writes <tag>_baseline.jsonl and <tag>_raw.jsonl.gz")
    args = ap.parse_args()

    if args.only:
        qids = list(args.only)
    else:
        splits = load_splits()
        qids = splits[args.split]
        if args.limit:
            qids = qids[: args.limit]

    base_url, key = pc.load_cfg()
    system = SYSTEM_TMPL.format(skill=pc.assemble_skill_context(args.skill_mode))
    print(f"[generate_baseline] skill_mode={args.skill_mode} desc_tier={args.desc_tier} "
          f"N={args.n_samples} T={args.temperature} qids={len(qids)} "
          f"system={len(system):,} chars", flush=True)

    base_out = OUT_DIR / f"{args.out_tag}_baseline.jsonl"
    raw_out = OUT_DIR / f"{args.out_tag}_raw.jsonl.gz"
    fp_base = base_out.open("w", encoding="utf-8")
    fp_raw = gzip.open(raw_out, "wt", encoding="utf-8")

    def _do(q):
        t = time.time()
        try:
            row, raw = run_qid(q, system=system, base=base_url, key=key,
                               n_samples=args.n_samples,
                               temperature=args.temperature,
                               desc_tier=args.desc_tier)
            with _print_lock:
                fp_base.write(json.dumps(row, ensure_ascii=False) + "\n"); fp_base.flush()
                for r in raw:
                    fp_raw.write(json.dumps(r, ensure_ascii=False) + "\n")
                fp_raw.flush()
            return q, time.time() - t, None
        except Exception:
            import traceback
            return q, time.time() - t, traceback.format_exc()

    t0 = time.time()
    ok = fail = 0
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_do, q): q for q in qids}
        for i, fut in enumerate(cf.as_completed(futs), 1):
            q, dt, err = fut.result()
            if err is None:
                ok += 1
                with _print_lock:
                    print(f"[{i}/{len(qids)}] OK {q} {dt:.1f}s | avg {(time.time()-t0)/i:.1f}s/q",
                          flush=True)
            else:
                fail += 1
                with _print_lock:
                    print(f"[{i}/{len(qids)}] ERR {q}  {err}", flush=True)

    fp_base.close(); fp_raw.close()
    print(f"[generate_baseline] done ok={ok} fail={fail} elapsed={time.time()-t0:.0f}s "
          f"-> {base_out} + {raw_out}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
