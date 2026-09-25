"""STEP 2 — Run the description2DSL-BOSS skill on each description, collect logprobs.

Reads batch/descriptions/<stem>.txt, feeds it (with the FULL skill context) to the
DSL model with logprobs on, and writes batch/runs/<stem>.json (answer + reasoning +
content_logprobs + usage) and batch/runs/<stem>.rb (generated DSL).

Resumable: skips a stem whose runs/<stem>.json already exists (and has tokens).
Parallel: --workers N.

    python token_conf/run_skill_batch.py --limit 3          # smoke test
    python token_conf/run_skill_batch.py --workers 6        # full
    python token_conf/run_skill_batch.py --only 1310.4101v3
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import threading
import time
from pathlib import Path

import pipeline_common as pc

DSL_MODEL = "deepseek-ai/deepseek-v4.1-flash"   # HepAI; returns logprobs only when streaming
MAX_TOKENS = 131072

SYSTEM_TMPL = """You are executing the description2DSL-BOSS skill. Follow its rules verbatim.
Output MUST be a single Ruby code block, no prose before or after.

The skill package (SKILL.md + references) follows. Treat it as authoritative:

{skill}
"""

_print_lock = threading.Lock()


def _has_tokens(p: Path) -> bool:
    try:
        return bool(json.loads(p.read_text()).get("content_logprobs"))
    except Exception:
        return False


def one(stem: str, base: str, key: str, system: str) -> tuple[str, int]:
    desc = (pc.DESC_DIR / f"{stem}.txt").read_text().strip()
    res = pc.chat(
        system=system, user=desc,
        model=DSL_MODEL, base_url=base, api_key=key,
        max_tokens=MAX_TOKENS, temperature=0.0,
        logprobs=True, top_logprobs=20,
    )
    ntok = len(res["content_logprobs"])
    if ntok == 0:
        raise RuntimeError(f"no content_logprobs (finish={res['finish_reason']}; "
                           "reasoning may have eaten the budget)")
    out = {
        "stem": stem, "model": DSL_MODEL, "user_description": desc,
        "answer": res["answer"], "reasoning": res["reasoning"],
        "finish_reason": res["finish_reason"],
        "content_logprobs": res["content_logprobs"],
        "usage": res["usage"],
    }
    (pc.RUN_DIR / f"{stem}.json").write_text(json.dumps(out, ensure_ascii=False))
    (pc.RUN_DIR / f"{stem}.rb").write_text(pc.strip_code_fence(res["answer"]))
    return stem, ntok


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--abort-after", type=int, default=30,
                    help="stop the run after this many CONSECUTIVE failures (balance/outage guard)")
    args = ap.parse_args()

    pc.RUN_DIR.mkdir(parents=True, exist_ok=True)
    base, key = pc.load_cfg()       # HepAI (streaming path yields real logprobs)
    system = SYSTEM_TMPL.format(skill=pc.assemble_skill_context())

    cand = [s for s in pc.stems(args.limit, args.only)
            if (pc.DESC_DIR / f"{s}.txt").exists()]
    todo = [s for s in cand
            if args.force or not _has_tokens(pc.RUN_DIR / f"{s}.json")]
    missing_desc = len(pc.stems(args.limit, args.only)) - len(cand)
    print(f"[step2] {len(cand)} have descriptions ({missing_desc} missing desc), "
          f"{len(todo)} to run, {len(cand)-len(todo)} cached, "
          f"workers={args.workers}, model={DSL_MODEL}, sys={len(system):,} chars, "
          f"abort_after={args.abort_after} consec fails")

    done = ok = fail = 0
    consec = 0
    aborted = False
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(one, s, base, key, system): s for s in todo}
        for fut in cf.as_completed(futs):
            s = futs[fut]
            done += 1
            try:
                _, ntok = fut.result()
                ok += 1
                consec = 0
                with _print_lock:
                    print(f"[{done}/{len(todo)}] OK  {s}  ({ntok} tok) "
                          f"| {(time.time()-t0)/done:.1f}s/it")
            except Exception as e:
                fail += 1
                consec += 1
                with _print_lock:
                    print(f"[{done}/{len(todo)}] ERR {s}: {e}")
                if consec >= args.abort_after and not aborted:
                    aborted = True
                    with _print_lock:
                        print(f"[step2] ABORT: {consec} consecutive failures "
                              f"(API balance/outage?). Cancelling remaining work; "
                              f"resume later — completed runs are cached.")
                    for f2 in futs:
                        f2.cancel()
    print(f"[step2] {'ABORTED' if aborted else 'done'}: ok={ok} fail={fail} "
          f"cached={len(cand)-len(todo)} elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
