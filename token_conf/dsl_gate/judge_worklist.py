"""Physics-Equivalence Accuracy (Phys-Acc): LLM-judged semantic metric; also
builds the judging worklist reused by judge_phys_acc.py.

The strict (dataset_em / decay_card_em) and structural (method-F1 / slot-F1)
metrics in evaluator.py deliberately ignore concrete numeric values. Phys-Acc
is the third tier: it asks an LLM judge whether the generated program and the
gold program specify *the same physics measurement* — same final state, same
tag/signal topology, same dataset, and selection thresholds within a range a
BESIII analyst would accept. It is the only metric in this paper that is
sensitive to numeric cut values.

Design decisions (all disclosed in the paper):

  * JUDGE MODEL is `deepseek-ai/deepseek-v4.1-flash`, the same model that
    generated the candidates. Self-preference bias is therefore possible, so we
    re-judge a subsample with an independent-vendor model (`zhipu/glm-5.2`) and
    report cross-judge agreement instead of asserting the bias away.
  * PROMPT follows G-Eval-style form filling: the judge must write its
    comparison into an `analysis` field BEFORE emitting the score, which
    recovers most of the chain-of-thought benefit without paying for a
    reasoning stream.
  * TEMPERATURE 0 for reproducibility.
  * POSITION BIAS. The gold is presented as PROGRAM_A or PROGRAM_B at random
    (seeded per pair), and the order is recorded. A `--swap-check N` subsample
    is additionally judged in BOTH orders so we can report verdict-flip rate.
  * SCOPE. We judge every (qid, candidate) pair plus every baseline generation,
    so Phys-Acc is available for *any* selection policy (including an
    oracle-by-Phys-Acc) without further LLM calls.
  * RESUMABLE. Verdicts are cached in judgments.jsonl keyed by
    (qid, kind, cand_stem, order); reruns skip cached keys.

Usage:
    python token_conf/dsl_gate/judge_worklist.py --workers 6
    python token_conf/dsl_gate/judge_worklist.py --limit 20        # smoke
    python token_conf/dsl_gate/judge_worklist.py --summarize-only  # no LLM calls
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _prereq import require  # noqa: E402

import pipeline_common as pc                       # noqa: E402
from pipeline_common import strip_code_fence       # noqa: E402


GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"       # GOLD_A (skill-regenerated corpus)
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results" / "phys_acc"
RESULTS.mkdir(parents=True, exist_ok=True)

JUDGE_MODEL = "deepseek-ai/deepseek-v4.1-flash"   # same family as the generator
ALT_JUDGE_MODEL = "zhipu/glm-5.2"                 # independent judge, for agreement
JUDGE_TEMP = 0.0
JUDGE_MAX_TOKENS = 1500
JUDGE_TIMEOUT = 300.0

# Both candidate judges emit a multi-thousand-token reasoning stream that shares
# `max_tokens`; with a small budget it starves the answer entirely (measured:
# 12/22 empty replies at max_tokens=1200 on glm, and the same failure mode on
# flash in generate_candidates). Disabling it drops per-call latency from ~53s to ~3.4s. To keep
# the chain-of-thought benefit G-Eval relies on, the judge is instead forced to
# emit its comparison *before* the score inside the JSON (form-filling).
NO_REASONING = {"reasoning_effort": "none"}

# NOTE ON BIAS. The primary judge is the SAME model family as the generator, so
# self-enhancement / self-preference bias is a live concern (Panickssery et al.,
# NeurIPS 2024). We therefore re-judge a subsample with ALT_JUDGE_MODEL (a
# different vendor) and report the cross-judge agreement, rather than merely
# asserting that the bias is absent.

JUDGMENTS = RESULTS / "judgments.jsonl"

_lock = threading.Lock()


JUDGE_SYSTEM = """You are a senior BESIII collaboration analyst. You review whether two BOSS \
analysis specifications, written in a Ruby domain-specific language, describe THE SAME PHYSICS \
MEASUREMENT.

You judge physics, not style. Ignore: variable names, comment wording, formatting, ordering of \
independent statements, and free-text `note(...)` blocks. Pay attention to: the parent and daughter \
particles and the full decay chain, the final-state particles reconstructed, the tag/signal \
topology (single-tag vs double-tag, which side is tagged), the dataset / energy point used, and \
whether the selection thresholds (track and photon quality, PID, kinematic-fit chi-square, mass \
windows) are close enough that the two programs would measure the same quantity.

Selection thresholds need not be identical. BESIII analyses legitimately differ in cut values. \
Judge whether the difference is inside the range a competent analyst would accept for the same \
measurement, or whether it changes what is being measured.

Respond with a single JSON object and nothing else. Fill the fields IN ORDER: state your \
comparison first, then the per-dimension findings, and only then the overall score.

{
  "analysis":          "<2-4 sentences: name the decay chain, final state, topology and dataset of each program, then say where they agree and where they differ>",
  "same_decay_chain":  true | false,
  "same_final_state":  true | false,
  "same_topology":     true | false,
  "same_dataset":      true | false,
  "cuts_acceptable":   true | false,
  "score":             1 | 2 | 3 | 4 | 5,
  "reason":            "<one sentence justifying the score>"
}

Score rubric:
  5 — physically equivalent; differences are cosmetic or well inside analyst tolerance.
  4 — same measurement; minor selection differences that shift efficiency but not the measured quantity.
  3 — same physics goal, but materially different selection; the two would not give consistent results.
  2 — related but different measurement (e.g. a different decay mode of the same parent).
  1 — unrelated measurement.
"""

JUDGE_USER = """PROGRAM_A:
```ruby
{prog_a}
```

PROGRAM_B:
```ruby
{prog_b}
```

Do PROGRAM_A and PROGRAM_B specify the same physics measurement? Reply with the JSON object only."""


# ------------------------- judging -------------------------------------


def _parse_verdict(text: str) -> dict | None:
    """Extract the JSON object from the judge's reply, tolerating code fences."""
    if not text:
        return None
    t = text.strip()
    t = re.sub(r"^```(?:json)?\s*", "", t)
    t = re.sub(r"\s*```$", "", t)
    # take the outermost {...}
    i, j = t.find("{"), t.rfind("}")
    if i == -1 or j == -1 or j <= i:
        return None
    try:
        obj = json.loads(t[i:j + 1])
    except json.JSONDecodeError:
        return None
    if not isinstance(obj, dict) or "score" not in obj:
        return None
    try:
        obj["score"] = int(obj["score"])
    except (TypeError, ValueError):
        return None
    if not 1 <= obj["score"] <= 5:
        return None
    return obj


def judge_pair(gold_ruby: str, pred_ruby: str, *, gold_first: bool,
               base: str, key: str, alt_judge: bool = False) -> dict:
    prog_a, prog_b = ((gold_ruby, pred_ruby) if gold_first
                      else (pred_ruby, gold_ruby))
    user = JUDGE_USER.format(prog_a=prog_a, prog_b=prog_b)
    model = ALT_JUDGE_MODEL if alt_judge else JUDGE_MODEL
    last_err = None
    for _ in range(2):
        try:
            r = pc.chat(system=JUDGE_SYSTEM, user=user, model=model,
                        base_url=base, api_key=key,
                        max_tokens=JUDGE_MAX_TOKENS, temperature=JUDGE_TEMP,
                        logprobs=False, stream=False, timeout=JUDGE_TIMEOUT,
                        extra_payload=NO_REASONING)
            v = _parse_verdict(r.get("answer") or "")
            if v is not None:
                v["_raw_len"] = len(r.get("answer") or "")
                return v
            last_err = f"unparseable(finish={r.get('finish_reason')})"
        except Exception as e:  # network / gateway hiccup
            last_err = repr(e)[:200]
    return {"score": None, "error": last_err or "unknown"}


# ------------------------- work list -----------------------------------


def build_worklist(limit: int | None = None,
                   cand_file: Path | None = None,
                   baseline_file: Path | None = None,
                   raw_file: Path | None = None) -> list[dict]:
    """One task per (qid, kind, cand_stem, sample_idx).

    Default behaviour (raw_file=None): one task per (qid, kind, cand_stem)
    keyed on `gen_ruby_first` from the summary jsonl — matches the historic
    main-table protocol.

    raw_file supplied: expand across every sample_idx by reading the gzipped
    per-sample answers written by generate_candidates. This is what the two-level
    (qid x sample_idx) paired bootstrap needs to estimate significance
    without collapsing the sampling variance to zero.
    """
    cf_path = cand_file or (GATE / "all.jsonl")
    bf_path = baseline_file or (GATE / "all_baseline.jsonl")
    cand_rows = [json.loads(l) for l in cf_path.open()]
    base_rows = [json.loads(l) for l in bf_path.open()]

    tasks: list[dict] = []
    seen: set[tuple] = set()

    if raw_file is not None:
        import gzip
        # (qid, kind, cand_stem) -> {sample_idx -> ruby}
        by_key: dict[tuple, dict[int, str]] = {}
        opener = gzip.open if str(raw_file).endswith(".gz") else open
        try:
            for l in opener(raw_file, "rt", encoding="utf-8"):
                r = json.loads(l)
                key = (r["qid"], r["kind"], r.get("cand_stem"))
                by_key.setdefault(key, {})[int(r["sample_idx"])] = \
                    strip_code_fence(r.get("answer", "") or "")
        except (EOFError, OSError) as e:
            # Truncated gzip (generate_candidates killed mid-write) — use what we have. Callers
            # can spot-check by counting rows per qid before scoring.
            print(f"[build_worklist] raw file truncated: {type(e).__name__}: {e}", flush=True)
        # cand_rank / gold flag come from the summary rows
        rank_by_key = {(r["qid"], "cand", r["cand_stem"]): r for r in cand_rows}
        for (qid, kind, cs), samples in by_key.items():
            meta = rank_by_key.get((qid, kind, cs), {})
            for si, ruby in samples.items():
                if not ruby:
                    continue
                tasks.append({
                    "qid": qid, "kind": kind, "cand_stem": cs,
                    "cand_rank": meta.get("cand_rank"),
                    "sample_idx": si, "ruby": ruby,
                })
    else:
        for b in base_rows:
            qid = b["qid"]
            k = (qid, "baseline", None)
            if k in seen:
                continue
            seen.add(k)
            tasks.append({"qid": qid, "kind": "baseline", "cand_stem": None,
                          "cand_rank": None, "sample_idx": 0,
                          "ruby": strip_code_fence(b.get("gen_ruby_first", "") or "")})
        for r in cand_rows:
            qid, cs = r["qid"], r.get("cand_stem")
            k = (qid, "cand", cs)
            if k in seen:
                continue
            seen.add(k)
            tasks.append({"qid": qid, "kind": "cand", "cand_stem": cs,
                          "cand_rank": r.get("cand_rank"), "sample_idx": 0,
                          "ruby": strip_code_fence(r.get("gen_ruby_first", "") or "")})

    # keep only qids that have a gold file
    tasks = [t for t in tasks if (GOLD / f"{t['qid']}.rb").exists()]
    if limit:
        keep_qids = sorted({t["qid"] for t in tasks})[:limit]
        tasks = [t for t in tasks if t["qid"] in keep_qids]
    return tasks


def _cache_key(t: dict, gold_first: bool, alt_judge: bool = False) -> str:
    suffix = "|ALT" if alt_judge else ""
    return (f"{t['qid']}|{t['kind']}|{t['cand_stem']}|"
            f"{'GF' if gold_first else 'PF'}{suffix}")


def load_cache() -> dict[str, dict]:
    cache: dict[str, dict] = {}
    if JUDGMENTS.exists():
        for l in JUDGMENTS.open():
            l = l.strip()
            if not l:
                continue
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                continue
            if r.get("score") is not None:
                cache[r["key"]] = r
    return cache


# ------------------------- run -----------------------------------------


def run(workers: int, limit: int | None, swap_check: int,
        cross_judge_check: int = 0) -> None:
    base, key = pc.load_cfg()
    tasks = build_worklist(limit)
    cache = load_cache()

    rng = random.Random(0)
    # primary pass: one randomized presentation order per pair
    jobs: list[tuple[dict, bool, bool]] = []
    for t in tasks:
        gold_first = rng.random() < 0.5
        if _cache_key(t, gold_first) not in cache:
            jobs.append((t, gold_first, False))

    # swap-consistency subsample: the SAME pairs judged in the opposite order
    swap_pool = [t for t in tasks if t["kind"] == "cand"]
    rng2 = random.Random(1)
    rng2.shuffle(swap_pool)
    for t in swap_pool[:swap_check]:
        for gf in (True, False):
            if _cache_key(t, gf) not in cache:
                jobs.append((t, gf, False))

    # cross-judge subsample: same pairs, gold-first, judged by ALT_JUDGE_MODEL
    if cross_judge_check:
        rng3 = random.Random(2)
        pool = [t for t in tasks if t["kind"] == "cand"]
        rng3.shuffle(pool)
        for t in pool[:cross_judge_check]:
            if _cache_key(t, True, alt_judge=True) not in cache:
                jobs.append((t, True, True))
            if _cache_key(t, True) not in cache:
                jobs.append((t, True, False))

    print(f"[phys-acc] tasks={len(tasks)} cached={len(cache)} jobs={len(jobs)} "
          f"judge={JUDGE_MODEL} alt={ALT_JUDGE_MODEL} workers={workers} "
          f"swap_check={swap_check} cross_judge_check={cross_judge_check}", flush=True)
    if not jobs:
        print("[phys-acc] nothing to do")
        return

    fp = JUDGMENTS.open("a", encoding="utf-8")
    t0 = time.time()
    ok = fail = 0

    def _do(job) -> tuple[dict, bool, bool, dict]:
        t, gold_first, alt_judge = job
        gold = (GOLD / f"{t['qid']}.rb").read_text(encoding="utf-8")
        v = judge_pair(gold, t["ruby"], gold_first=gold_first,
                       base=base, key=key, alt_judge=alt_judge)
        return t, gold_first, alt_judge, v

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_do, j) for j in jobs]
        for i, fut in enumerate(cf.as_completed(futs), 1):
            try:
                t, gold_first, alt_judge, v = fut.result()
            except Exception as e:
                fail += 1
                print(f"[{i}/{len(jobs)}] EXC {e!r}", flush=True)
                continue
            row = {"key": _cache_key(t, gold_first, alt_judge), "qid": t["qid"],
                   "kind": t["kind"], "cand_stem": t["cand_stem"],
                   "cand_rank": t["cand_rank"], "gold_first": gold_first,
                   "alt_judge": alt_judge, **v}
            with _lock:
                fp.write(json.dumps(row, ensure_ascii=False) + "\n")
                fp.flush()
            if v.get("score") is None:
                fail += 1
            else:
                ok += 1
            if i % 100 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] ok={ok} fail={fail} "
                      f"{el:.0f}s ({el/i:.2f}s/job)", flush=True)
    fp.close()
    print(f"[phys-acc] done ok={ok} fail={fail} wall={time.time()-t0:.0f}s")


# ------------------------- summarize ----------------------------------


def _score_map(cache: dict[str, dict]) -> dict[tuple, list[dict]]:
    """(qid, kind, cand_stem) -> list of verdicts across presentation orders.

    Alt-judge rows are held out: they exist only to measure cross-judge
    agreement and must not enter the reported metric."""
    m: dict[tuple, list[dict]] = defaultdict(list)
    for r in cache.values():
        if r.get("alt_judge"):
            continue
        m[(r["qid"], r["kind"], r.get("cand_stem"))].append(r)
    return m


def summarize() -> None:
    cache = load_cache()
    if not cache:
        print("[phys-acc] no judgments yet")
        return
    m = _score_map(cache)

    cand_rows = [json.loads(l) for l in require(GATE / "all.jsonl").open()]
    by_qid: dict[str, list[dict]] = defaultdict(list)
    for r in cand_rows:
        by_qid[r["qid"]].append(r)

    def verdict(qid, kind, cand_stem) -> float | None:
        """Mean score across available orders (position-bias averaged)."""
        vs = m.get((qid, kind, cand_stem))
        if not vs:
            return None
        ss = [v["score"] for v in vs if v.get("score") is not None]
        return st.mean(ss) if ss else None

    # policies
    def pick_stem(strat: str, qid: str, cands: list[dict]) -> str | None:
        if strat == "reranker_top1":
            return min(cands, key=lambda c: c.get("cand_rank", 999)).get("cand_stem")
        if strat.startswith("entropy"):
            a = float(strat.split("=")[1].rstrip(")")) if "=" in strat else 0.5
            def s(c):
                return a * (c.get("dH_full") or 0.0) + (1 - a) * (c.get("dH_sel") or 0.0)
            return max(cands, key=s).get("cand_stem")
        if strat == "random":
            return random.Random(42).choice(cands).get("cand_stem")
        if strat == "oracle_phys":
            best, bs = None, -1.0
            for c in cands:
                v = verdict(qid, "cand", c.get("cand_stem"))
                if v is not None and v > bs:
                    bs, best = v, c.get("cand_stem")
            return best
        raise ValueError(strat)

    strategies = ["no_rag", "reranker_top1", "entropy(a=0)", "entropy(a=0.5)",
                  "entropy(a=1)", "random", "oracle_phys"]
    table = []
    for strat in strategies:
        scores = []
        for qid in sorted(by_qid):
            if strat == "no_rag":
                v = verdict(qid, "baseline", None)
            else:
                cs = pick_stem(strat, qid, by_qid[qid])
                v = verdict(qid, "cand", cs) if cs else None
            if v is not None:
                scores.append(v)
        if not scores:
            continue
        acc = sum(1 for s in scores if s >= 4.0) / len(scores)
        table.append({"strategy": strat, "phys_acc": acc,
                      "mean_score": st.mean(scores), "n": len(scores)})

    # position-bias: pairs judged in both orders
    flips, both = 0, 0
    delta = []
    for k, vs in m.items():
        gf = [v for v in vs if v.get("gold_first") and v.get("score") is not None]
        pf = [v for v in vs if not v.get("gold_first") and v.get("score") is not None]
        if gf and pf:
            both += 1
            a, b = gf[0]["score"], pf[0]["score"]
            delta.append(a - b)
            if (a >= 4) != (b >= 4):
                flips += 1

    # sub-dimension agreement rates (over all cand verdicts)
    dims = ["same_decay_chain", "same_final_state", "same_topology",
            "same_dataset", "cuts_acceptable"]
    dim_rate = {}
    for d in dims:
        vals = [bool(r[d]) for r in cache.values()
                if r.get("kind") == "cand" and isinstance(r.get(d), bool)]
        if vals:
            dim_rate[d] = sum(vals) / len(vals)

    lines = ["# Phys-Acc. — LLM-judged physics equivalence", "",
             f"- Judge: `{JUDGE_MODEL}` (temperature {JUDGE_TEMP}), reasoning stream disabled, "
             f"form-filling prompt (analysis emitted before score).",
             f"- The judge shares the generator's model family, so self-preference bias is "
             f"possible; we quantify it with a cross-judge subsample against "
             f"`{ALT_JUDGE_MODEL}` below rather than assuming it away.",
             f"- Gold: `token_conf/batch/runs/<qid>.rb` (GOLD_A, skill-regenerated corpus).",
             f"- Verdicts cached: {len(cache)} over {len(m)} unique (qid, kind, candidate) pairs.",
             f"- Phys-Acc. := fraction of queries whose selected program scores $\\geq 4$ on the 1–5 rubric.",
             "", "## Policy comparison", "",
             "| policy | Phys-Acc. | mean score (1–5) | n |", "|---|---|---|---|"]
    for r in table:
        lines.append(f"| {r['strategy']} | {r['phys_acc']*100:.1f}% | "
                     f"{r['mean_score']:.3f} | {r['n']} |")

    if table:
        nr = next((r for r in table if r["strategy"] == "no_rag"), None)
        rr = next((r for r in table if r["strategy"] == "reranker_top1"), None)
        orc = next((r for r in table if r["strategy"] == "oracle_phys"), None)
        if nr and rr and orc and (orc["phys_acc"] - nr["phys_acc"]) > 0:
            lines += ["", "## Oracle-gap recovery on Phys-Acc.", "",
                      "| policy | recovery |", "|---|---|"]
            gap = orc["phys_acc"] - nr["phys_acc"]
            for r in table:
                rec = (r["phys_acc"] - nr["phys_acc"]) / gap * 100
                lines.append(f"| {r['strategy']} | {rec:+.1f}% |")

    lines += ["", "## Judge reliability", ""]
    if both:
        lines.append(f"- Pairs judged in BOTH presentation orders: {both}")
        lines.append(f"- Binary verdict flips (score$\\geq$4 changes): {flips}/{both} "
                     f"({flips/both*100:.1f}%)")
        if delta:
            lines.append(f"- Mean signed score difference (gold-first − pred-first): "
                         f"{st.mean(delta):+.3f} (a value near 0 indicates no position bias)")
    else:
        lines.append("- No swap-order pairs judged yet (run with `--swap-check N`).")

    # cross-judge agreement (primary judge vs independent-vendor judge)
    alt = {(r["qid"], r.get("cand_stem")): r for r in cache.values()
           if r.get("alt_judge") and r.get("score") is not None and r.get("gold_first")}
    pri = {(r["qid"], r.get("cand_stem")): r for r in cache.values()
           if not r.get("alt_judge") and r.get("score") is not None and r.get("gold_first")}
    shared = sorted(set(alt) & set(pri))
    if shared:
        agree_bin = sum(1 for k in shared
                        if (alt[k]["score"] >= 4) == (pri[k]["score"] >= 4))
        d = [pri[k]["score"] - alt[k]["score"] for k in shared]
        exact = sum(1 for x in d if x == 0)
        bias = st.mean(d)
        lines += ["",
                  f"- Cross-judge subsample (`{JUDGE_MODEL}` vs `{ALT_JUDGE_MODEL}`): "
                  f"n={len(shared)}",
                  f"  - binary-verdict agreement: {agree_bin}/{len(shared)} "
                  f"({agree_bin/len(shared)*100:.1f}%)",
                  f"  - exact 1–5 score agreement: {exact}/{len(shared)} "
                  f"({exact/len(shared)*100:.1f}%)",
                  f"  - mean signed difference (primary − alt): {bias:+.3f}",
                  "  - A positive difference is the magnitude of the primary judge's",
                  "    self-preference: it scores programs written by its own family higher",
                  "    than an independent judge does. Since the same judge scores every",
                  "    policy, this shifts all rows together and does not by itself explain",
                  "    a difference *between* policies; but the absolute Phys-Acc. level",
                  "    should be read with this offset in mind."]

    n_fail = sum(1 for l in JUDGMENTS.open()
                 if l.strip() and json.loads(l).get("score") is None)
    lines.append(f"- Unparseable / failed judgments: {n_fail}")

    if dim_rate:
        lines += ["", "## Sub-dimension true-rate over all candidate verdicts", "",
                  "| dimension | rate |", "|---|---|"]
        for d, v in dim_rate.items():
            lines.append(f"| `{d}` | {v*100:.1f}% |")

    out = RESULTS / "summary.md"
    out.write_text("\n".join(lines))
    print("\n".join(lines))
    print(f"\n[+] {out}")

    # spot-check file for human review
    rows = [r for r in cache.values() if r.get("kind") == "cand"]
    random.Random(7).shuffle(rows)
    sc = ["# Phys-Acc. spot-check sample (for human review)", ""]
    for r in rows[:30]:
        sc.append(f"## {r['qid']}  ← cand {r['cand_stem']} (rank {r.get('cand_rank')})")
        sc.append(f"- score **{r['score']}**, gold_first={r['gold_first']}")
        for d in dims:
            if d in r:
                sc.append(f"  - {d}: {r[d]}")
        sc.append(f"- reason: {r.get('reason', '')}")
        sc.append("")
    (RESULTS / "spotcheck.md").write_text("\n".join(sc))
    print(f"[+] {RESULTS / 'spotcheck.md'}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None,
                    help="only judge the first N qids (smoke test)")
    ap.add_argument("--swap-check", type=int, default=80,
                    help="how many candidate pairs to judge in BOTH orders")
    ap.add_argument("--cross-judge-check", type=int, default=60,
                    help="how many candidate pairs to additionally judge with "
                         "ALT_JUDGE_MODEL, to measure self-preference bias")
    ap.add_argument("--summarize-only", action="store_true")
    args = ap.parse_args()

    if not args.summarize_only:
        run(args.workers, args.limit, args.swap_check, args.cross_judge_check)
    summarize()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
