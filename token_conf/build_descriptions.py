"""STEP 1 — Reverse-engineer a concise Chinese description from each gold .rb.

For every arxiv_papers_dsl/<stem>.rb we ask the LLM to produce ONE compact
paragraph (like the hand-written "样例描述 4") that a physicist would write to
*commission* that DSL: physics channel, datasets, exclusive-MC, and the event-
selection chain with concrete cut values — no Ruby, no method names.

Resumable: skips a stem whose descriptions/<stem>.txt already exists.
Parallel: --workers N threads.

    python token_conf/build_descriptions.py --limit 3            # smoke test
    python token_conf/build_descriptions.py --workers 8          # full 862
    python token_conf/build_descriptions.py --only 1310.4101v3   # one paper
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import threading
import time
from pathlib import Path

import pipeline_common as pc

DESC_MODEL = pc.model_id("deepseek-ai/deepseek-v4.1-flash")

DESC_SYSTEM = """You reverse-engineer a BESIII analysis DESCRIPTION from its BOSS Ruby DSL.

Given the Ruby DSL of one analysis, write the concise natural-language description
a physicist would hand to a colleague to COMMISSION exactly this BOSS event-selection
code. Write in ENGLISH (physics/particle jargon and symbols as usual).

RULES:
- ONE compact paragraph, ~4-8 sentences. No headings, no lists, no Ruby, no DSL
  method names (do NOT write select_track / kinematic_fit / DatasetManager etc.).
- DO include, in physicist's words: the physics channel & decay chain; the datasets
  (energy points, real data / inclusive MC) and exclusive-MC (events per mode, which
  decay modes); and the event selection with CONCRETE cut values — charged-track cuts
  (|cosθ|, |Vz|, Vr, track/charge counts), photon cuts (energy thresholds, count),
  PID strategy (method, lepton/π-K separation, particle counts), and the kinematic
  fit (nC constraint, χ² cut, mass windows).
- Mirror the style of this reference description:
  "On ψ(4260), study e+e− → γ X(3872), X(3872) → π+π−J/ψ, J/ψ → e+e− / μ+μ−. Use
   4.260 GeV real data and inclusive MC; generate 100k-event exclusive MC for each of
   the two J/ψ decay modes. Event selection: charged tracks with |cosθ|<0.93, |Vz|<10
   cm, Vr<1 cm, two positive and two negative tracks, net charge 0; photons with 25
   MeV (barrel) / 50 MeV (endcap) energy threshold, at least one; PID by the
   probability method, high-momentum tracks (p>1.0) treated as leptons (electron if
   EMC energy >0.6, else muon), with π/K separation, requiring 1 π+, 1 π−, 1 l+, 1 l−.
   Perform a 4C kinematic fit (γ π+π− l+l−), χ² cut 60, flagged nominal. Both decay
   modes share the same selection chain."
- Describe only what the DSL actually encodes. Do not invent physics not in the code.
- Output ONLY the paragraph, nothing else."""

_print_lock = threading.Lock()


def one(stem: str, base: str, key: str) -> tuple[str, str]:
    rb = (pc.RB_DIR / f"{stem}.rb").read_text()
    res = pc.chat(
        system=DESC_SYSTEM,
        user=f"Ruby DSL of analysis `{stem}`:\n\n```ruby\n{rb}\n```",
        model=DESC_MODEL, base_url=base, api_key=key,
        max_tokens=131072, temperature=0.0, logprobs=False,  # reasoning model: max headroom for CoT
    )
    desc = res["answer"].strip()
    if not desc:
        raise RuntimeError(f"empty description (finish={res['finish_reason']})")
    (pc.DESC_DIR / f"{stem}.txt").write_text(desc)
    return stem, desc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--force", action="store_true", help="regenerate even if .txt exists")
    ap.add_argument("--abort-after", type=int, default=12,
                    help="stop the run after this many CONSECUTIVE failures (balance/outage guard)")
    args = ap.parse_args()

    pc.DESC_DIR.mkdir(parents=True, exist_ok=True)
    base, key = pc.load_cfg()
    todo = [s for s in pc.stems(args.limit, args.only)
            if args.force or not (pc.DESC_DIR / f"{s}.txt").exists()]
    total = len(pc.stems(args.limit, args.only))
    print(f"[step1] {total} targets, {len(todo)} to build, {total-len(todo)} cached, "
          f"workers={args.workers}, model={DESC_MODEL}, abort_after={args.abort_after} consec fails")

    done = ok = fail = 0
    consec = 0
    aborted = False
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(one, s, base, key): s for s in todo}
        for fut in cf.as_completed(futs):
            s = futs[fut]
            done += 1
            try:
                _, desc = fut.result()
                ok += 1
                consec = 0
                with _print_lock:
                    print(f"[{done}/{len(todo)}] OK  {s}  ({len(desc)} chars) "
                          f"| {(time.time()-t0)/done:.1f}s/it")
            except Exception as e:
                fail += 1
                consec += 1
                with _print_lock:
                    print(f"[{done}/{len(todo)}] ERR {s}: {e}")
                if consec >= args.abort_after and not aborted:
                    aborted = True
                    with _print_lock:
                        print(f"[step1] ABORT: {consec} consecutive failures "
                              f"(API balance/outage?). Cancelling remaining work; "
                              f"resume later — completed descriptions are cached.")
                    for f2 in futs:
                        f2.cancel()
    print(f"[step1] {'ABORTED' if aborted else 'done'}: ok={ok} fail={fail} "
          f"cached={total-len(todo)} elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
