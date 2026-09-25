"""Candidate generation: sample N DSL programs per target and score their entropy.

For each query stem in a split:
  1. search_dsl(desc, top_k=5)          → 5 candidate ruby files
  2. run baseline prompt (no ruby example) N=3 times     → H_full_base, H_sel_base
  3. for each candidate:
       run +ruby prompt N=3 times                        → H_full_cand, H_sel_cand
       ΔH_full = H_full_base - H_full_cand   (larger = candidate helped more)
       ΔH_sel  = H_sel_base  - H_sel_cand

Output: token_conf/dsl_gate/runs_gate/<split>.jsonl
Row schema (one row per (qid, candidate)):
  {
    "qid": stem,
    "cand_stem": stem,
    "cand_rank": int,              # 1..K reranker rank
    "cand_score": float,
    "H_full_base_mean": float,     # mean over N samples of per-sample mean H
    "H_sel_base_mean":  float,
    "H_full_cand_mean": float,
    "H_sel_cand_mean":  float,
    "dH_full": float,
    "dH_sel":  float,
    "n_samples": N,
    "gen_ruby_first": str,         # first candidate sample (for evaluator)
    "usage": {"prompt_tokens": .., "completion_tokens": ..},
    "wall_ms": float,
    "n_tok_base_mean": float,
    "n_tok_cand_mean": float,
  }

Also writes a "meta" JSON alongside {qid} holding the baseline samples so the
evaluator (evaluator) can score No-RAG later:
  runs_gate/<split>_baseline.jsonl
  {"qid": stem, "H_full_base_mean": .., "H_sel_base_mean": .., "gen_ruby_first": ..}
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import statistics
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "token_conf"))

import pipeline_common as pc          # noqa: E402
from compute_entropy import tag_tokens, token_entropy   # noqa: E402
from segment_modules import tag_modules                 # noqa: E402
from token_conf.dsl_gate.search_dsl import search_dsl   # noqa: E402
from token_conf.dsl_gate.splits import load_splits      # noqa: E402


DSL_MODEL = pc.model_id("deepseek-ai/deepseek-v4.1-flash")
# reasoning_effort="none" (probed) makes flash skip its reasoning stream, so
# max_tokens governs answer tokens only; gold answer p99 ≈ 7.6k tokens.
MAX_TOKENS = 12000
DEFAULT_TEMPERATURE = 0.7  # >0 required for meaningful logprobs on this gateway
REQUEST_TIMEOUT = 600.0
DEFAULT_N = 3
DEFAULT_K = 5

# Extra payload for pc.chat(): disable the reasoning stream on flash.
EXTRA_PAYLOAD = {"reasoning_effort": "none"}

OUT_DIR = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SYSTEM_TMPL = """You are executing the description2DSL-BOSS skill. Follow its rules verbatim.
Output MUST be a single Ruby code block, no prose before or after.

The skill package (SKILL.md + references) follows. Treat it as authoritative:

{skill}
"""

USER_BASELINE_TMPL = """Description:
{desc}

Generate the BOSS DSL for the analysis described above.
"""

USER_WITH_EXAMPLE_TMPL = """A retrieved reference DSL for a related BESIII analysis is shown below.
Use it as a stylistic and structural hint; do NOT copy values that do not match
the target description. Adapt cuts, particle modes, and dataset choices to the
target description.

--- reference DSL (from arXiv:{ref_stem}) ---
{ref_rb}
--- end reference ---

Description:
{desc}

Generate the BOSS DSL for the target description.
"""

_print_lock = threading.Lock()


# -------------------- entropy summarization ------------------------------


def _summarize(content_logprobs: list[dict]) -> dict:
    """One LLM call -> per-sample mean H_full + mean H_selection (nats).

    tag_tokens/tag_modules mutate their input in-place, so we shallow-copy each
    token dict first — safe against concurrent workers sharing the same list
    across threads (which would otherwise race on 'kind'/'kept'/'module' keys
    and produce 'list has no attribute keys' when a partially-tagged token is
    reprocessed)."""
    if not content_logprobs:
        return {"H_full": None, "H_sel": None, "n_code": 0, "n_sel": 0, "n_tok": 0}
    toks = [dict(t) for t in content_logprobs]
    tag_tokens(toks)     # in-place: kind/kept
    tag_modules(toks)    # in-place: module
    Hs_all: list[float] = []
    Hs_sel: list[float] = []
    for t in toks:
        if not t.get("kept"):
            continue
        H = token_entropy(t)["H_tail_nats"]
        Hs_all.append(H)
        if t.get("module") == "selection":
            Hs_sel.append(H)
    return {
        "H_full": (sum(Hs_all) / len(Hs_all)) if Hs_all else None,
        "H_sel": (sum(Hs_sel) / len(Hs_sel)) if Hs_sel else None,
        "n_code": len(Hs_all),
        "n_sel": len(Hs_sel),
        "n_tok": len(toks),
    }


def _mean(xs: list[float | None]) -> float | None:
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


# -------------------- LLM call --------------------------------------------


def _call_once(system: str, user: str, base: str, key: str, temperature: float) -> dict:
    """One call with a single retry on empty content."""
    for _ in range(2):
        res = pc.chat(
            system=system, user=user,
            model=DSL_MODEL, base_url=base, api_key=key,
            max_tokens=MAX_TOKENS, temperature=temperature,
            logprobs=True, top_logprobs=20,
            timeout=REQUEST_TIMEOUT,
            extra_payload=EXTRA_PAYLOAD,
        )
        if res.get("content_logprobs"):
            return res
    return res


def _run_n(system: str, user: str, base: str, key: str, n: int, temperature: float) -> list[dict]:
    """N sequential samples. The gateway serializes streams per API key, so
    intra-query parallelism just queues on the server side — outer per-query
    parallelism (--workers) is what actually scales."""
    return [_call_once(system, user, base, key, temperature) for _ in range(n)]


# -------------------- per-query pipeline ----------------------------------


def process_query(
    qid: str,
    *,
    system: str,
    base: str,
    key: str,
    top_k: int,
    n_samples: int,
    temperature: float,
    desc_dir: Path | None = None,
    desc_tier: str = "full",
    include_gold: bool = False,
) -> tuple[list[dict], dict, list[dict]]:
    """Return (candidate_rows, baseline_row, raw_rows). raw_rows carries the
    full per-sample content_logprobs + answer, one dict per (qid, kind, cand,
    sample_idx). Written to a separate gzipped jsonl so the summary file stays
    small enough for interactive tools.

    desc_tier="full"  — use the reverse-engineered ~200-word description as-is.
    desc_tier="bare"  — collapse to "energy point + reaction arrow" (leak-free).
                        Retrieval still uses the FULL desc; only the prompt to
                        the DSL generator changes.
    include_gold      — additionally emit a synthetic candidate with the query's
                        own gold .rb as reference (cand_rank=0, cand_stem=qid).
                        This is the oracle-gold pass; use for the oracle row.
    """
    dd = desc_dir if desc_dir is not None else (ROOT / "data" / "corpus" / "descriptions")
    desc_path = dd / f"{qid}.txt"
    full_desc = desc_path.read_text(encoding="utf-8").strip()
    prompt_desc = pc.bare_desc(full_desc) if desc_tier == "bare" else full_desc

    # ---- retrieval always uses the full description (dense retriever needs signal) ----
    # The KB (besiii_dsl_desc_v1) contains every gold DSL including this qid's own,
    # so exclude_stems={qid} is what actually enforces the RAG contract.
    hits = search_dsl(full_desc, top_k=top_k, exclude_stems={qid})

    # ---- oracle-gold: prepend the query's own gold DSL as rank-0 candidate ----
    if include_gold:
        gold_path = ROOT / "data" / "corpus" / "generated_dsl" / f"{qid}.rb"
        if gold_path.exists():
            hits = [{"stem": qid, "rb_text": gold_path.read_text(encoding="utf-8"),
                     "score": float("inf"), "is_gold": True}] + hits

    raw_rows: list[dict] = []

    # ---- baseline (no reference): shared across candidates ----
    t0 = time.time()
    base_user = USER_BASELINE_TMPL.format(desc=prompt_desc)
    base_calls = _run_n(system, base_user, base, key, n_samples, temperature)
    base_summaries = [_summarize(c["content_logprobs"]) for c in base_calls]
    base_H_full = _mean([s["H_full"] for s in base_summaries])
    base_H_sel = _mean([s["H_sel"] for s in base_summaries])
    base_ntok = _mean([float(s["n_tok"]) for s in base_summaries])
    for si, c in enumerate(base_calls):
        raw_rows.append({
            "qid": qid, "kind": "baseline", "cand_stem": None,
            "cand_rank": None, "sample_idx": si,
            "answer": c.get("answer", ""),
            "content_logprobs": c.get("content_logprobs") or [],
            "usage": c.get("usage"),
        })
    baseline_row = {
        "qid": qid,
        "H_full_base_mean": base_H_full,
        "H_sel_base_mean": base_H_sel,
        "n_tok_base_mean": base_ntok,
        "n_samples": n_samples,
        "gen_ruby_first": pc.strip_code_fence(base_calls[0]["answer"]),
        "wall_ms_base": (time.time() - t0) * 1000,
    }

    # ---- candidates ----
    cand_rows: list[dict] = []
    # gold candidate (if present) uses rank=0 so ordinary rank-1..K numbering is preserved
    _enum_start = 0 if (include_gold and hits and hits[0].get("is_gold")) else 1
    for rank, h in enumerate(hits, _enum_start):
        t1 = time.time()
        cand_user = USER_WITH_EXAMPLE_TMPL.format(
            ref_stem=h["stem"], ref_rb=h["rb_text"], desc=prompt_desc,
        )
        cand_calls = _run_n(system, cand_user, base, key, n_samples, temperature)
        cand_summaries = [_summarize(c["content_logprobs"]) for c in cand_calls]
        cand_H_full = _mean([s["H_full"] for s in cand_summaries])
        cand_H_sel = _mean([s["H_sel"] for s in cand_summaries])
        cand_ntok = _mean([float(s["n_tok"]) for s in cand_summaries])
        for si, c in enumerate(cand_calls):
            raw_rows.append({
                "qid": qid, "kind": "cand", "cand_stem": h["stem"],
                "cand_rank": rank, "sample_idx": si,
                "answer": c.get("answer", ""),
                "content_logprobs": c.get("content_logprobs") or [],
                "usage": c.get("usage"),
            })
        row = {
            "qid": qid,
            "cand_stem": h["stem"],
            "cand_rank": rank,
            "cand_score": h["score"] if h["score"] != float("inf") else None,
            "H_full_base_mean": base_H_full,
            "H_sel_base_mean": base_H_sel,
            "H_full_cand_mean": cand_H_full,
            "H_sel_cand_mean": cand_H_sel,
            "dH_full": (base_H_full - cand_H_full)
                       if (base_H_full is not None and cand_H_full is not None) else None,
            "dH_sel":  (base_H_sel - cand_H_sel)
                       if (base_H_sel is not None and cand_H_sel is not None) else None,
            "n_samples": n_samples,
            "gen_ruby_first": pc.strip_code_fence(cand_calls[0]["answer"]),
            "n_tok_base_mean": base_ntok,
            "n_tok_cand_mean": cand_ntok,
            "wall_ms_cand": (time.time() - t1) * 1000,
        }
        if h.get("is_gold"):
            row["is_gold"] = True
        cand_rows.append(row)
    return cand_rows, baseline_row, raw_rows


# -------------------- driver ----------------------------------------------


def _already_done(out_path: Path) -> set[str]:
    done: set[str] = set()
    if out_path.exists():
        with out_path.open() as f:
            for line in f:
                try:
                    done.add(json.loads(line)["qid"])
                except Exception:
                    continue
    return done


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="dev",
                    help="split name in splits.json (dev/test/train/all/custom)")
    ap.add_argument("--limit", type=int, default=None,
                    help="only run the first N queries in the split (smoke)")
    ap.add_argument("--only", nargs="*", default=None,
                    help="override split entirely with these stems")
    ap.add_argument("--n-samples", type=int, default=DEFAULT_N)
    ap.add_argument("--top-k", type=int, default=DEFAULT_K)
    ap.add_argument("--temperature", type=float, default=DEFAULT_TEMPERATURE,
                    help="sampling temperature for the generator; default 0.7 "
                         "reproduces the main-table run")
    ap.add_argument("--workers", type=int, default=4,
                    help="parallel queries (each internally serial)")
    ap.add_argument("--force", action="store_true",
                    help="ignore existing rows in output jsonl")
    ap.add_argument("--out-tag", default=None,
                    help="override output filename stem (default = split name); "
                         "writes to runs_gate/<tag>{,_baseline,_raw.jsonl.gz}")
    ap.add_argument("--desc-dir", type=Path, default=None,
                    help="override description directory (default = token_conf/batch/descriptions); "
                         "used by perturbation experiments (paraphrase / etc.)")
    ap.add_argument("--desc-tier", choices=["full", "bare"], default="full",
                    help="'full' = reverse-engineered ~200-word description (main-table default, "
                         "leaks cut values into No-RAG); 'bare' = energy + reaction only, RAG-headroom clean")
    ap.add_argument("--include-gold", action="store_true",
                    help="also emit an oracle-gold candidate (rank=0, cand_stem=qid) whose "
                         "reference DSL is the query's own gold .rb — used for the oracle row")
    ap.add_argument("--skill-mode", choices=list(pc.SKILL_MODES), default="full",
                    help="skill package mode for BOTH baseline and cand system prompts. "
                         "Default 'full' reproduces the main-table run. 'no-example-clean-manual' "
                         "removes example1/2/3 AND redacts the CppTemplate defaults block in "
                         "dsl_manual.md (cos_theta/Vz/Vr/EMC thresholds) into typed placeholders.")
    ap.add_argument("--api-key", default=None,
                    help="override the API key from config.yaml (base_url stays the same). "
                         "Use to shard a large run across multiple keys — each shard is one "
                         "process with its own --only qid list + --api-key + --out-tag.")
    args = ap.parse_args()

    if args.only:
        qids = list(args.only)
        tag = "only"
    else:
        splits = load_splits()
        qids = splits[args.split]
        tag = args.split
        if args.limit:
            qids = qids[: args.limit]

    file_tag = args.out_tag or tag
    out_path = OUT_DIR / f"{file_tag}.jsonl"
    base_out = OUT_DIR / f"{file_tag}_baseline.jsonl"
    raw_out = OUT_DIR / f"{file_tag}_raw.jsonl.gz"
    done = set() if args.force else _already_done(out_path)
    todo = [q for q in qids if q not in done]
    print(f"[generate_candidates] split={tag} qids={len(qids)} done={len(done)} todo={len(todo)} "
          f"K={args.top_k} N={args.n_samples} workers={args.workers}", flush=True)

    base_url, key = pc.load_cfg()
    if args.api_key:
        key = args.api_key
    system = SYSTEM_TMPL.format(skill=pc.assemble_skill_context(args.skill_mode))
    print(f"[generate_candidates] system prompt {len(system):,} chars; model={DSL_MODEL} temp={args.temperature} "
          f"skill_mode={args.skill_mode} key=...{key[-4:]}", flush=True)

    import gzip
    t0 = time.time()
    fp_cand = out_path.open("a", encoding="utf-8")
    fp_base = base_out.open("a", encoding="utf-8")
    fp_raw = gzip.open(raw_out, "at", encoding="utf-8")

    def _do_one(q: str) -> tuple[str, int, float, str | None]:
        t = time.time()
        try:
            cand_rows, base_row, raw_rows = process_query(
                q, system=system, base=base_url, key=key,
                top_k=args.top_k, n_samples=args.n_samples,
                temperature=args.temperature,
                desc_dir=args.desc_dir,
                desc_tier=args.desc_tier,
                include_gold=args.include_gold,
            )
            with _print_lock:
                for r in cand_rows:
                    fp_cand.write(json.dumps(r, ensure_ascii=False) + "\n")
                fp_cand.flush()
                fp_base.write(json.dumps(base_row, ensure_ascii=False) + "\n")
                fp_base.flush()
                for r in raw_rows:
                    fp_raw.write(json.dumps(r, ensure_ascii=False) + "\n")
                fp_raw.flush()
            return q, len(cand_rows), time.time() - t, None
        except Exception:
            import traceback
            return q, 0, time.time() - t, traceback.format_exc()

    ok = fail = 0
    try:
        with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_do_one, q): q for q in todo}
            for i, fut in enumerate(cf.as_completed(futs), 1):
                q, ncand, dt, err = fut.result()
                if err is None:
                    ok += 1
                    with _print_lock:
                        print(f"[{i}/{len(todo)}] OK  {q}  cands={ncand}  {dt:.1f}s "
                              f"| avg {(time.time()-t0)/i:.1f}s/q", flush=True)
                else:
                    fail += 1
                    with _print_lock:
                        print(f"[{i}/{len(todo)}] ERR {q}  {dt:.1f}s  {err}", flush=True)
    finally:
        fp_cand.close()
        fp_base.close()
        fp_raw.close()

    print(f"[generate_candidates] done split={tag}: ok={ok} fail={fail} "
          f"elapsed={time.time()-t0:.0f}s -> {out_path}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
