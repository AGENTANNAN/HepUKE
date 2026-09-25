"""Phys-Acc-100: strict multi-aspect LLM-judged physics equivalence.

Replaces the coarse 1-5 rubric of judge_worklist.py (which saturated at ~91% and
could not resolve policy differences) with a seven-dimension review scored 0-100
per dimension. The judge scores dimensions only; the weighted total is computed
HERE, in code, because LLM arithmetic is unreliable and because the per-dimension
breakdown is itself a result (it tells us which part of the specification a
retrieved reference actually helps with).

Dimensions and weights (sum 100):
    decay_chain_and_final_state  25
    dataset_and_energy           15
    topology_and_tagging         15
    track_photon_selection       15
    pid                          10
    fits_and_constraints         10
    mc_and_decay_card            10

A dimension the REFERENCE does not address at all may be scored `null`; weights
are then renormalised over the applicable dimensions, so a missing dimension is
not silently charged as zero. (A dimension the reference DOES address but the
candidate omits is a 0, not a null — this is stated in the prompt.)

Presentation is deliberately ASYMMETRIC (REFERENCE vs CANDIDATE) rather than
randomised A/B: the task is asymmetric, and telling the judge which program is
the accepted specification makes the strict rubric far more stable. The cost is
that the position-bias control of judge_worklist does not carry over; we instead keep the
cross-judge agreement check.

Usage:
    python token_conf/dsl_gate/judge_phys_acc.py --limit 20 --workers 8
    python token_conf/dsl_gate/judge_phys_acc.py --summarize-only
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import random
import re
import statistics as st
import sys
import threading
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

import pipeline_common as pc                      # noqa: E402
from pipeline_common import strip_code_fence      # noqa: E402
from judge_worklist import build_worklist           # noqa: E402  (same task list)


GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"     # GOLD_A
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results" / "phys_acc100"
RESULTS.mkdir(parents=True, exist_ok=True)
CAND_FILE_OVERRIDE: Path | None = None
JUDGMENTS = RESULTS / "judgments.jsonl"

JUDGE_MODEL = pc.model_id("deepseek-ai/deepseek-v4.1-flash")
ALT_JUDGE_MODEL = "zhipu/glm-5.2"
JUDGE_TEMP = 0.0
JUDGE_MAX_TOKENS = 2500          # richer output than judge_worklist: evidence lists + 7 scores
JUDGE_TIMEOUT = 300.0
NO_REASONING = {"reasoning_effort": "none"}

WEIGHTS = {
    "decay_chain_and_final_state": 25,
    "dataset_and_energy": 15,
    "topology_and_tagging": 15,
    "track_photon_selection": 15,
    "pid": 10,
    "fits_and_constraints": 10,
    "mc_and_decay_card": 10,
}
DIMS = list(WEIGHTS)

_lock = threading.Lock()


JUDGE_SYSTEM = """You are a senior BESIII collaboration analyst performing a strict technical review. \
You are given two BOSS analysis specifications written in a Ruby domain-specific language: the \
REFERENCE (the accepted specification) and the CANDIDATE (a specification to be reviewed). You must \
score how faithfully the CANDIDATE reproduces the physics of the REFERENCE, dimension by dimension.

BE STRICT. Your default assumption is NOT that a difference is acceptable. You are looking for \
reasons the two specifications would NOT yield the same published result. A specification that \
"looks like" the reference but selects different events is a failure, not a near-miss.

What counts as a substantive difference:
  * ANY difference in a discrete choice is substantive: a different particle, a different decay \
mode, a different tagged side, a different PID hypothesis set, a different fit type, a different \
generator model (PHSP vs SVS vs HELAMP ...), a different energy point or dataset id.
  * For a continuous threshold (|cos(theta)|, Vz, Vr, energy thresholds, isolation angles, \
probability cuts, chi-square cuts, mass windows), a relative difference above 20% is substantive. \
Below 20% it is a minor difference. An omitted threshold that the reference imposes is \
substantive, not minor.
  * Missing a requirement the reference imposes is always worse than stating it with a slightly \
different value.

What does NOT count and must be ignored entirely:
  * Variable names, comment text and language, indentation, line order among independent \
statements, and the content of free-text `note(...)` blocks.

Score each dimension on 0-100 using these anchors:
  100 - substantively identical; any difference is purely notational.
   80 - equivalent in physics; only minor threshold differences (all under 20% relative).
   60 - same intent, but ONE substantive difference, or several minor ones that together would \
visibly change the selected sample.
   40 - partially overlapping: the dimension is addressed but with multiple substantive errors.
   20 - the dimension is addressed in name only; what is specified is largely wrong.
    0 - the dimension is absent from the CANDIDATE, or directly contradicts the REFERENCE.

The dimensions:
  decay_chain_and_final_state - the parent particle, the full decay chain, and the set of \
final-state particles actually reconstructed.
  dataset_and_energy - which real-data and inclusive-MC samples are used, and at which energy \
point / centre-of-mass energy; dataset identifiers must agree.
  topology_and_tagging - single-tag vs double-tag, which side is tagged versus signal-side, \
recoil-mass tagging versus direct reconstruction.
  track_photon_selection - charged-track quality requirements (|cos(theta)|, Vz, Vr) and \
multiplicity/net-charge requirements; photon requirements (energy thresholds by barrel and \
endcap, EMC timing window, isolation angles, photon multiplicity).
  pid - the identification method, the probability or likelihood cut, which hypotheses are tested \
against which, and the resulting particle multiplicity requirements.
  fits_and_constraints - kinematic fit and/or vertex fit: which fit, over which particle list, \
what is constrained, and the chi-square cut.
  mc_and_decay_card - the exclusive-MC configuration and the EvtGen decay card: the decay chain \
declared, the generator model on each node, and the requested event count.

If a dimension is genuinely not addressed by the REFERENCE at all (for example the reference \
declares no exclusive MC), return null for that dimension rather than a number. Do not return null \
merely because the CANDIDATE omits something the REFERENCE has - that is a score of 0.

Respond with a single JSON object and nothing else. Fill the fields IN ORDER: record the evidence \
first, then the per-dimension scores, then the summary fields.

{
  "reference_summary": "<1-2 sentences: the measurement the REFERENCE specifies>",
  "candidate_summary": "<1-2 sentences: the measurement the CANDIDATE specifies>",
  "substantive_differences": ["<one short bullet per substantive difference you found; empty list if none>"],
  "minor_differences":       ["<one short bullet per minor (sub-20%) threshold difference>"],
  "scores": {
    "decay_chain_and_final_state": <0-100 or null>,
    "dataset_and_energy":          <0-100 or null>,
    "topology_and_tagging":        <0-100 or null>,
    "track_photon_selection":      <0-100 or null>,
    "pid":                         <0-100 or null>,
    "fits_and_constraints":        <0-100 or null>,
    "mc_and_decay_card":           <0-100 or null>
  },
  "fatal_mismatch": true | false,
  "verdict": "<one sentence>"
}

Set fatal_mismatch to true if the two specifications measure fundamentally different quantities - \
a different decay chain or a different parent particle - such that no adjustment of thresholds \
could reconcile them."""


JUDGE_USER = """REFERENCE specification:
```ruby
{reference}
```

CANDIDATE specification:
```ruby
{candidate}
```

Review the CANDIDATE against the REFERENCE. Reply with the JSON object only."""


# ---------------------- scoring ----------------------------------------


def weighted_total(scores: dict) -> float | None:
    """Renormalise WEIGHTS over the dimensions the judge actually scored."""
    num = den = 0.0
    for d, w in WEIGHTS.items():
        v = scores.get(d)
        if v is None:
            continue
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        v = max(0.0, min(100.0, v))
        num += w * v
        den += w
    if den == 0:
        return None
    return num / den


def _parse(text: str) -> dict | None:
    if not text:
        return None
    t = re.sub(r"^```(?:json)?\s*", "", text.strip())
    t = re.sub(r"\s*```$", "", t)
    i, j = t.find("{"), t.rfind("}")
    if i == -1 or j <= i:
        return None
    try:
        obj = json.loads(t[i:j + 1])
    except json.JSONDecodeError:
        return None
    if not isinstance(obj, dict) or not isinstance(obj.get("scores"), dict):
        return None
    tot = weighted_total(obj["scores"])
    if tot is None:
        return None
    obj["total"] = tot
    return obj


def judge(reference: str, candidate: str, *, base: str, key: str,
          alt: bool = False) -> dict:
    user = JUDGE_USER.format(reference=reference, candidate=candidate)
    model = ALT_JUDGE_MODEL if alt else JUDGE_MODEL
    err = None
    for _ in range(2):
        try:
            r = pc.chat(system=JUDGE_SYSTEM, user=user, model=model,
                        base_url=base, api_key=key,
                        max_tokens=JUDGE_MAX_TOKENS, temperature=JUDGE_TEMP,
                        logprobs=False, stream=False, timeout=JUDGE_TIMEOUT,
                        extra_payload=NO_REASONING)
            v = _parse(r.get("answer") or "")
            if v is not None:
                return v
            err = f"unparseable(finish={r.get('finish_reason')})"
        except Exception as e:
            err = repr(e)[:200]
    return {"total": None, "error": err or "unknown"}


# ---------------------- run --------------------------------------------


def _ck(t: dict, alt: bool = False) -> str:
    si = t.get("sample_idx", 0)
    return f"{t['qid']}|{t['kind']}|{t['cand_stem']}|s{si}" + ("|ALT" if alt else "")


def load_cache() -> dict[str, dict]:
    out: dict[str, dict] = {}
    if JUDGMENTS.exists():
        for l in JUDGMENTS.open():
            l = l.strip()
            if not l:
                continue
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                continue
            if r.get("total") is not None:
                out[r["key"]] = r
    return out


def run(workers: int, limit: int | None, cross_judge_check: int,
        cand_file: Path | None = None,
        baseline_file: Path | None = None,
        raw_file: Path | None = None) -> None:
    base, key = pc.load_cfg()
    tasks = build_worklist(limit, cand_file=cand_file,
                           baseline_file=baseline_file, raw_file=raw_file)
    cache = load_cache()

    jobs = [(t, False) for t in tasks if _ck(t) not in cache]
    if cross_judge_check:
        pool = [t for t in tasks if t["kind"] == "cand"]
        random.Random(2).shuffle(pool)
        for t in pool[:cross_judge_check]:
            if _ck(t, True) not in cache:
                jobs.append((t, True))

    print(f"[phys100] tasks={len(tasks)} cached={len(cache)} jobs={len(jobs)} "
          f"judge={JUDGE_MODEL} workers={workers} cross={cross_judge_check}",
          flush=True)
    if not jobs:
        print("[phys100] nothing to do")
        return

    fp = JUDGMENTS.open("a", encoding="utf-8")
    t0 = time.time()
    ok = fail = 0

    def _do(job):
        t, alt = job
        ref = (GOLD / f"{t['qid']}.rb").read_text(encoding="utf-8")
        return t, alt, judge(ref, t["ruby"], base=base, key=key, alt=alt)

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        for i, fut in enumerate(cf.as_completed([ex.submit(_do, j) for j in jobs]), 1):
            try:
                t, alt, v = fut.result()
            except Exception as e:
                fail += 1
                print(f"[{i}/{len(jobs)}] EXC {e!r}", flush=True)
                continue
            row = {"key": _ck(t, alt), "qid": t["qid"], "kind": t["kind"],
                   "cand_stem": t["cand_stem"], "cand_rank": t["cand_rank"],
                   "sample_idx": t.get("sample_idx", 0),
                   "alt_judge": alt, **v}
            with _lock:
                fp.write(json.dumps(row, ensure_ascii=False) + "\n")
                fp.flush()
            ok += v.get("total") is not None
            fail += v.get("total") is None
            if i % 50 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] ok={ok} fail={fail} {el:.0f}s "
                      f"({el/i:.2f}s/job)", flush=True)
    fp.close()
    print(f"[phys100] done ok={ok} fail={fail} wall={time.time()-t0:.0f}s")


# ---------------------- summarize --------------------------------------


def summarize(inspect: int = 0) -> None:
    cache = load_cache()
    if not cache:
        print("[phys100] no judgments yet")
        return
    # Group judgments by (qid, kind, cand_stem); values are per-sample_idx dicts.
    # Historic behaviour (only sample 0 was judged) still yields one entry per group.
    by_group: dict[tuple, dict[int, dict]] = defaultdict(dict)
    for r in cache.values():
        if r.get("alt_judge"):
            continue
        gk = (r["qid"], r["kind"], r["cand_stem"])
        by_group[gk][int(r.get("sample_idx", 0))] = r

    def mean_total(gk):
        rows = list(by_group.get(gk, {}).values())
        vals = [r["total"] for r in rows if r.get("total") is not None]
        return sum(vals) / len(vals) if vals else None

    def fatal_frac(gk):
        rows = list(by_group.get(gk, {}).values())
        if not rows:
            return None
        return sum(bool(r.get("fatal_mismatch")) for r in rows) / len(rows)

    # "pri" preserves the legacy single-total view (mean across whatever samples exist)
    pri = {gk: {"total": mean_total(gk),
                "fatal_mismatch": (fatal_frac(gk) or 0) >= 0.5,
                "scores": next(iter(rows.values()))["scores"] if (rows := by_group[gk]) else {}}
           for gk in by_group}

    cand_rows = [json.loads(l) for l in (CAND_FILE_OVERRIDE or (GATE / "all.jsonl")).open()]
    by_qid: dict[str, list[dict]] = defaultdict(list)
    for r in cand_rows:
        by_qid[r["qid"]].append(r)
    judged_qids = sorted({q for (q, _, _) in pri})

    def tot(qid, kind, cs):
        r = pri.get((qid, kind, cs))
        return r["total"] if r else None

    def stem_for(strat, qid):
        cs = by_qid.get(qid) or []
        real = [c for c in cs if not c.get("is_gold") and c.get("cand_rank") != 0]
        if not cs:
            return None
        if strat == "reranker_top1":
            return min(real or cs, key=lambda c: c.get("cand_rank", 999))["cand_stem"]
        if strat.startswith("entropy"):
            a = float(strat.split("=")[1].rstrip(")"))
            return max(real or cs, key=lambda c: a * (c.get("dH_full") or 0.0)
                       + (1 - a) * (c.get("dH_sel") or 0.0))["cand_stem"]
        if strat == "random":
            return random.Random(42).choice(real or cs)["cand_stem"]
        if strat == "oracle_phys100":
            best, bs = None, -1.0
            for c in (real or cs):
                v = tot(qid, "cand", c["cand_stem"])
                if v is not None and v > bs:
                    bs, best = v, c["cand_stem"]
            return best
        if strat == "oracle_gold":
            for c in cs:
                if c.get("is_gold") or c.get("cand_rank") == 0:
                    return c["cand_stem"]
            return None
        raise ValueError(strat)

    strategies = ["no_rag", "reranker_top1",
                  "entropy(a=0)", "entropy(a=0.3)", "entropy(a=0.5)",
                  "entropy(a=0.7)", "entropy(a=1)",
                  "random", "oracle_phys100", "oracle_gold"]
    lines = ["# Phys-Acc-100 — strict multi-aspect physics equivalence", "",
             f"- Judge `{JUDGE_MODEL}`, temperature {JUDGE_TEMP}, reasoning disabled, "
             f"asymmetric REFERENCE/CANDIDATE presentation.",
             f"- Weighted total computed in code over "
             f"{len(WEIGHTS)} dimensions; weights renormalised over dimensions the "
             f"judge marked applicable.",
             f"- Judged: {len(pri)} unique (qid, kind, candidate) triples over "
             f"{len(judged_qids)} queries.", "",
             "## Policy comparison", "",
             "| policy | mean total (0–100) | median | fatal-mismatch rate | n |",
             "|---|---|---|---|---|"]
    for s in strategies:
        vals, fatals = [], []
        for q in judged_qids:
            v = tot(q, "baseline", None) if s == "no_rag" else \
                tot(q, "cand", stem_for(s, q))
            if v is None:
                continue
            vals.append(v)
            r = pri.get((q, "baseline", None) if s == "no_rag"
                        else (q, "cand", stem_for(s, q)))
            fatals.append(bool(r.get("fatal_mismatch")))
        if not vals:
            continue
        lines.append(f"| {s} | {st.mean(vals):.2f} | {st.median(vals):.2f} | "
                     f"{sum(fatals)/len(fatals)*100:.1f}% | {len(vals)} |")

    # per-dimension means over all candidate verdicts
    lines += ["", "## Per-dimension mean over all candidate verdicts", "",
              "| dimension | weight | mean | n scored | n null |", "|---|---|---|---|---|"]
    for d in DIMS:
        vs = [r["scores"].get(d) for r in pri.values() if r.get("kind") == "cand"]
        nums = [float(v) for v in vs if isinstance(v, (int, float))]
        nulls = sum(1 for v in vs if v is None)
        if nums:
            lines.append(f"| `{d}` | {WEIGHTS[d]} | {st.mean(nums):.1f} | "
                         f"{len(nums)} | {nulls} |")

    # cross-judge agreement
    alt = {(r["qid"], r["cand_stem"]): r for r in cache.values() if r.get("alt_judge")}
    both = [k for k in alt if (k[0], "cand", k[1]) in pri]
    if both:
        d = [pri[(q, "cand", c)]["total"] - alt[(q, c)]["total"] for q, c in both]
        lines += ["", "## Cross-judge agreement", "",
                  f"- n={len(both)}, mean signed difference (primary − alt): {st.mean(d):+.2f} "
                  f"points on the 0–100 scale",
                  f"- mean absolute difference: {st.mean([abs(x) for x in d]):.2f} points"]

    n_fail = sum(1 for l in JUDGMENTS.open()
                 if l.strip() and json.loads(l).get("total") is None)
    lines.append(f"\n- Unparseable / failed judgments: {n_fail}")

    out = RESULTS / "summary.md"
    out.write_text("\n".join(lines))
    print("\n".join(lines))
    print(f"\n[+] {out}")

    # ---- human-readable inspection dump ----
    if inspect:
        rows = [r for r in pri.values() if r.get("kind") == "cand"]
        rows.sort(key=lambda r: r["total"])
        pick = rows[:inspect // 2] + rows[-(inspect - inspect // 2):]
        dump = ["# Phys-Acc-100 inspection sample (lowest and highest scoring)", ""]
        for r in pick:
            dump.append(f"## {r['qid']}  ← cand {r['cand_stem']} (rank {r.get('cand_rank')})  "
                        f"**total {r['total']:.1f}**"
                        + ("  ⚠ FATAL" if r.get("fatal_mismatch") else ""))
            dump.append(f"- REF: {r.get('reference_summary','')}")
            dump.append(f"- CAND: {r.get('candidate_summary','')}")
            sc = r.get("scores", {})
            dump.append("- scores: " + ", ".join(
                f"{d.split('_')[0]}={sc.get(d)}" for d in DIMS))
            sd = r.get("substantive_differences") or []
            md = r.get("minor_differences") or []
            if sd:
                dump.append("- substantive:")
                dump += [f"  - {x}" for x in sd]
            if md:
                dump.append("- minor:")
                dump += [f"  - {x}" for x in md]
            dump.append(f"- verdict: {r.get('verdict','')}")
            dump.append("")
        p = RESULTS / "inspection.md"
        p.write_text("\n".join(dump))
        print(f"[+] {p}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--cross-judge-check", type=int, default=0)
    ap.add_argument("--inspect", type=int, default=12,
                    help="dump this many extreme cases to inspection.md")
    ap.add_argument("--summarize-only", action="store_true")
    ap.add_argument("--cand-file", type=Path, default=None,
                    help="override candidate jsonl (default runs_gate/all.jsonl)")
    ap.add_argument("--baseline-file", type=Path, default=None,
                    help="override baseline jsonl (default runs_gate/all_baseline.jsonl)")
    ap.add_argument("--raw-file", type=Path, default=None,
                    help="if set, judge every (qid, kind, cand_stem, sample_idx) from "
                         "the gzipped raw jsonl instead of only sample_idx=0")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="override output dir (default results/phys_acc100/)")
    args = ap.parse_args()
    global RESULTS, JUDGMENTS, CAND_FILE_OVERRIDE
    if args.out_dir:
        RESULTS = args.out_dir
        RESULTS.mkdir(parents=True, exist_ok=True)
        JUDGMENTS = RESULTS / "judgments.jsonl"
    if args.cand_file:
        CAND_FILE_OVERRIDE = args.cand_file
    if not args.summarize_only:
        run(args.workers, args.limit, args.cross_judge_check,
            cand_file=args.cand_file, baseline_file=args.baseline_file,
            raw_file=args.raw_file)
    summarize(args.inspect)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
