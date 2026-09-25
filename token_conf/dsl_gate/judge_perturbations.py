"""Phys-Acc-100 judging for the perturbation variants.

Runs the Phys-Acc judge over the perturbation variants at T=0.1, writing a
separate judgment cache so it does not collide with the main run's.

Usage::
    python token_conf/dsl_gate/judge_perturbations.py --workers 8
    python token_conf/dsl_gate/judge_perturbations.py --summarize-only
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, json, statistics as st
import sys, threading, time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

import pipeline_common as pc
from pipeline_common import strip_code_fence
from variants import VARIANT_NAMES, tag
import judge_phys_acc as judge_phys_acc
from judge_phys_acc import judge
judge_phys_acc.JUDGE_MODEL = pc.model_id("deepseek-flash")
JUDGE_MODEL = judge_phys_acc.JUDGE_MODEL

GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results" / "phys_acc100"
RESULTS.mkdir(parents=True, exist_ok=True)
JUDGMENTS = RESULTS / "perturb_T01_judgments.jsonl"

# Splits to judge; a missing runs_gate tag is skipped with a warning.
SPLIT_PREFIXES = ("pert20", "pert50")

# (variant key, runs_gate tag). The key IS the tag — one name per split × variant.
VARIANTS = [(tag(v, prefix), tag(v, prefix))
            for prefix in SPLIT_PREFIXES for v in VARIANT_NAMES]

_lock = threading.Lock()


def _entropy_score(c): return 0.5 * (c.get("dH_full") or 0.0) + 0.5 * (c.get("dH_sel") or 0.0)


def _ck(t):
    return f"{t['variant']}|{t['qid']}|{t['kind']}|{t.get('cand_stem') or ''}"


def load_cache() -> dict:
    if not JUDGMENTS.exists():
        return {}
    out = {}
    for l in JUDGMENTS.open():
        l = l.strip()
        if not l:
            continue
        r = json.loads(l)
        out[r["key"]] = r
    return out


def build_worklist() -> list[dict]:
    tasks: list[dict] = []
    for vname, tag in VARIANTS:
        cand_p = GATE / f"{tag}.jsonl"
        base_p = GATE / f"{tag}_baseline.jsonl"
        if not cand_p.exists() or not base_p.exists():
            print(f"[skip] {vname}: missing files", flush=True); continue
        cands: dict[str, list[dict]] = defaultdict(list)
        for l in cand_p.open():
            r = json.loads(l); cands[r["qid"]].append(r)
        base = {json.loads(l)["qid"]: json.loads(l) for l in base_p.open()}
        for q in sorted(cands.keys() & base.keys()):
            # baseline (no_rag)
            b = base[q]
            tasks.append({"variant": vname, "qid": q, "kind": "baseline",
                          "cand_stem": None, "cand_rank": None,
                          "ruby": strip_code_fence(b.get("gen_ruby_first", "") or "")})
            cs = cands[q]
            picks = {"rel": min(cs, key=lambda c: c.get("cand_rank", 999)),
                     "ent": max(cs, key=_entropy_score)}
            seen = set()
            for _, r in picks.items():
                stem = r.get("cand_stem")
                if stem in seen:
                    continue
                seen.add(stem)
                tasks.append({"variant": vname, "qid": q, "kind": "cand",
                              "cand_stem": stem, "cand_rank": r.get("cand_rank"),
                              "ruby": strip_code_fence(r.get("gen_ruby_first", "") or "")})
    return tasks


def do_judgments(workers: int = 8) -> None:
    tasks = build_worklist()
    cache = load_cache()
    jobs = [t for t in tasks if _ck(t) not in cache]
    print(f"[judge_perturbations-T01] tasks={len(tasks)} cached={len(cache)} jobs={len(jobs)} "
          f"judge={JUDGE_MODEL} workers={workers}", flush=True)
    if not jobs:
        print("[judge_perturbations-T01] nothing to do"); return
    base, key = pc.load_cfg()

    fp = JUDGMENTS.open("a", encoding="utf-8")
    t0 = time.time(); ok = fail = 0

    def _do(t):
        ref = (GOLD / f"{t['qid']}.rb").read_text(encoding="utf-8")
        return t, judge(ref, t["ruby"], base=base, key=key)

    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_do, j) for j in jobs]
        for i, fut in enumerate(cf.as_completed(futs), 1):
            try:
                t, v = fut.result()
            except Exception as e:
                fail += 1
                print(f"[{i}/{len(jobs)}] EXC {e!r}", flush=True)
                continue
            row = {"key": _ck(t), "variant": t["variant"], "qid": t["qid"],
                   "kind": t["kind"], "cand_stem": t["cand_stem"],
                   "cand_rank": t.get("cand_rank"), **v}
            with _lock:
                fp.write(json.dumps(row, ensure_ascii=False) + "\n")
                fp.flush()
            ok += v.get("total") is not None
            fail += v.get("total") is None
            if i % 25 == 0 or i == len(jobs):
                el = time.time() - t0
                print(f"[{i}/{len(jobs)}] ok={ok} fail={fail} {el:.0f}s "
                      f"({el/i:.2f}s/job)", flush=True)
    fp.close()
    print(f"[judge_perturbations-T01] done ok={ok} fail={fail} wall={time.time()-t0:.0f}s")


def summarize() -> None:
    cache = load_cache()
    if not cache:
        print("[judge_perturbations-T01] no judgments yet"); return
    pri = {(r["variant"], r["qid"], r["kind"], r["cand_stem"]): r
           for r in cache.values()}
    cands: dict[tuple, list[dict]] = defaultdict(list)
    baselines: dict[tuple, dict] = {}
    for vname, tag in VARIANTS:
        cp = GATE / f"{tag}.jsonl"; bp = GATE / f"{tag}_baseline.jsonl"
        if not cp.exists() or not bp.exists():
            continue
        for l in cp.open():
            r = json.loads(l); cands[(vname, r["qid"])].append(r)
        for l in bp.open():
            b = json.loads(l); baselines[(vname, b["qid"])] = b
    print("| variant | n | no_rag | reranker_top1 | entropy(a=.5) |")
    print("|---|---|---|---|---|")
    for vname, _ in VARIANTS:
        qids = sorted({q for (v, q) in cands if v == vname})
        cells = []
        for strat in ("no_rag", "rel", "ent"):
            vals, n = [], 0
            for q in qids:
                if strat == "no_rag":
                    r = pri.get((vname, q, "baseline", None))
                else:
                    cs = cands.get((vname, q), [])
                    if not cs: continue
                    pick = (min(cs, key=lambda c: c.get("cand_rank", 999))
                            if strat == "rel" else max(cs, key=_entropy_score))
                    r = pri.get((vname, q, "cand", pick.get("cand_stem")))
                if r is None: continue
                n += 1; vals.append(r["total"])
            cells.append(f"{st.mean(vals):.2f} (n={n})" if vals else "—")
        print(f"| {vname} | {len(qids)} | {cells[0]} | {cells[1]} | {cells[2]} |")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--summarize-only", action="store_true")
    args = ap.parse_args()
    if not args.summarize_only:
        do_judgments(workers=args.workers)
    summarize()


if __name__ == "__main__":
    main()
