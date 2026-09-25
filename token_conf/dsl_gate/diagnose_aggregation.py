"""Aggregation diagnostic — quick A/B check: does swapping mean(H) for sum(H) / topK(H)
help the entropy_top1 selector beat reranker_top1?

Reruns generate_candidates on 5 dev qids but STORES the full per-sample content_logprobs.
Then computes several aggregation variants on the raw H_tail sequences and
reports which one best correlates ΔH-selected candidate with oracle-selected
candidate on the small probe set.

Usage:
    python token_conf/dsl_gate/diagnose_aggregation.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

import pipeline_common as pc                                    # noqa: E402
from compute_entropy import tag_tokens, token_entropy           # noqa: E402
from segment_modules import tag_modules                         # noqa: E402
from search_dsl import search_dsl                               # noqa: E402
from evaluator import score_pair                             # noqa: E402


PROBE_QIDS = [
    "1506.06018v2",
    "1309.1896v2",
    "2212.13072v2",
    "2206.08554v2",
    "2409.00427v2",
]
N_SAMPLES = 3
TOP_K = 5
DSL_MODEL = pc.model_id("deepseek-ai/deepseek-v4.1-flash")
MAX_TOKENS = 12000
TEMPERATURE = 0.7
EXTRA_PAYLOAD = {"reasoning_effort": "none"}

OUT_DIR = ROOT / "token_conf" / "dsl_gate" / "runs_gate_probe"
OUT_DIR.mkdir(parents=True, exist_ok=True)
RAW_OUT = OUT_DIR / "probe_raw.jsonl"

SYSTEM_TMPL = ("You are executing the description2DSL-BOSS skill. Follow its rules verbatim.\n"
               "Output MUST be a single Ruby code block, no prose before or after.\n"
               "\nThe skill package (SKILL.md + references) follows. Treat it as authoritative:\n\n{skill}")
USER_BASE = "Description:\n{desc}\n\nGenerate the BOSS DSL for the analysis described above.\n"
USER_REF = ("A retrieved reference DSL for a related BESIII analysis is shown below.\n"
            "Use it as a stylistic and structural hint; do NOT copy values that do not match\n"
            "the target description. Adapt cuts, particle modes, and dataset choices to the\n"
            "target description.\n\n--- reference DSL (from arXiv:{ref_stem}) ---\n{ref_rb}\n--- end reference ---\n\n"
            "Description:\n{desc}\n\nGenerate the BOSS DSL for the target description.\n")


# ---------------- collect raw (LLM heavy) ---------------------------------


def collect_raw() -> list[dict]:
    if RAW_OUT.exists():
        print(f"[+] loading cached probe raw from {RAW_OUT}")
        return [json.loads(l) for l in RAW_OUT.open()]

    base, key = pc.load_cfg()
    system = SYSTEM_TMPL.format(skill=pc.assemble_skill_context())
    print(f"[probe] system={len(system):,} chars; model={DSL_MODEL}; qids={len(PROBE_QIDS)}")

    def _call(user: str) -> dict:
        for _ in range(2):
            r = pc.chat(system=system, user=user, model=DSL_MODEL, base_url=base,
                        api_key=key, max_tokens=MAX_TOKENS, temperature=TEMPERATURE,
                        logprobs=True, top_logprobs=20, timeout=600.0,
                        extra_payload=EXTRA_PAYLOAD)
            if r.get("content_logprobs"):
                return r
        return r

    all_rows: list[dict] = []
    t0 = time.time()
    for qi, qid in enumerate(PROBE_QIDS, 1):
        desc = (ROOT / "data" / "corpus" / "descriptions" / f"{qid}.txt").read_text().strip()
        hits = search_dsl(desc, top_k=TOP_K + 1)
        hits = [h for h in hits if h["stem"] != qid][:TOP_K]

        # baseline (shared across candidates in prod, but we still run N samples)
        base_samples = []
        for si in range(N_SAMPLES):
            r = _call(USER_BASE.format(desc=desc))
            base_samples.append({
                "answer": r["answer"],
                "content_logprobs": r["content_logprobs"],
            })
        all_rows.append({"qid": qid, "kind": "baseline",
                         "samples": base_samples})

        for h in hits:
            cand_samples = []
            for si in range(N_SAMPLES):
                r = _call(USER_REF.format(ref_stem=h["stem"], ref_rb=h["rb_text"], desc=desc))
                cand_samples.append({
                    "answer": r["answer"],
                    "content_logprobs": r["content_logprobs"],
                })
            all_rows.append({"qid": qid, "kind": "cand",
                             "cand_stem": h["stem"], "cand_rank": hits.index(h) + 1,
                             "samples": cand_samples})
        print(f"[{qi}/{len(PROBE_QIDS)}] {qid} done  wall={time.time()-t0:.0f}s")

    with RAW_OUT.open("w") as f:
        for r in all_rows:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"[+] saved {len(all_rows)} raw rows -> {RAW_OUT}  ({RAW_OUT.stat().st_size/1e6:.1f} MB)")
    return all_rows


# ---------------- aggregation variants ------------------------------------


def per_sample_H_by_module(sample: dict) -> tuple[list[float], list[float]]:
    """Return (all_code_H, selection_only_H) lists for one sample."""
    toks = [dict(t) for t in sample["content_logprobs"]]
    if not toks:
        return [], []
    tag_tokens(toks)
    tag_modules(toks)
    Hs_all, Hs_sel = [], []
    for t in toks:
        if not t.get("kept"):
            continue
        H = token_entropy(t)["H_tail_nats"]
        Hs_all.append(H)
        if t.get("module") == "selection":
            Hs_sel.append(H)
    return Hs_all, Hs_sel


def aggregate(Hs: list[float], variant: str) -> float | None:
    if not Hs:
        return None
    if variant == "mean":
        return sum(Hs) / len(Hs)
    if variant == "sum":
        return sum(Hs)
    if variant.startswith("topk_mean_"):
        k = int(variant.split("_")[-1])
        top = sorted(Hs, reverse=True)[:k]
        return sum(top) / len(top) if top else None
    if variant.startswith("topk_sum_"):
        k = int(variant.split("_")[-1])
        top = sorted(Hs, reverse=True)[:k]
        return sum(top) if top else None
    if variant.startswith("frac_ge_"):
        thr = float(variant.split("_")[-1])
        return sum(1 for h in Hs if h >= thr) / len(Hs)
    if variant == "count_hot":         # count H > 0.5 nats
        return sum(1 for h in Hs if h > 0.5)
    if variant == "max":
        return max(Hs)
    raise ValueError(variant)


VARIANTS = [
    "mean", "sum",
    "topk_mean_10", "topk_mean_20", "topk_mean_50",
    "topk_sum_10", "topk_sum_20", "topk_sum_50",
    "frac_ge_0.5", "count_hot", "max",
]


# ---------------- probe experiment ----------------------------------------


def compute_dH(rows: list[dict], variant: str) -> dict[str, dict[str, dict]]:
    """qid -> {cand_stem or "_baseline_": {"full": val, "sel": val}}"""
    # 1) per-row aggregate over N samples: mean-of-per-sample-aggregate
    per_row: dict[tuple[str, str], dict[str, float]] = {}
    for r in rows:
        Hs_all_samples, Hs_sel_samples = [], []
        for s in r["samples"]:
            Hs_all, Hs_sel = per_sample_H_by_module(s)
            va = aggregate(Hs_all, variant)
            vs = aggregate(Hs_sel, variant)
            if va is not None:
                Hs_all_samples.append(va)
            if vs is not None:
                Hs_sel_samples.append(vs)
        key = (r["qid"], r.get("cand_stem", "_baseline_"))
        per_row[key] = {
            "full": (sum(Hs_all_samples) / len(Hs_all_samples)) if Hs_all_samples else None,
            "sel": (sum(Hs_sel_samples) / len(Hs_sel_samples)) if Hs_sel_samples else None,
        }
    # 2) dH per (qid, cand)
    dH: dict[str, dict[str, dict]] = {}
    for (qid, key), vals in per_row.items():
        if key == "_baseline_":
            continue
        base = per_row.get((qid, "_baseline_"))
        if base is None:
            continue
        entry = dH.setdefault(qid, {})
        entry[key] = {
            "dH_full": (base["full"] - vals["full"])
                        if (base["full"] is not None and vals["full"] is not None) else None,
            "dH_sel":  (base["sel"] - vals["sel"])
                        if (base["sel"] is not None and vals["sel"] is not None) else None,
        }
    return dH


def eval_variant(rows: list[dict], variant: str, alpha: float) -> dict:
    """For each qid pick the cand with max(alpha*dH_full + (1-alpha)*dH_sel);
    score against gold; report mean struct_mean and how often it agrees with
    oracle (highest struct_mean among the 5)."""
    dH = compute_dH(rows, variant)
    # first pass: for each qid, gather (cand, first-sample ruby, gold-score)
    per_q_cands: dict[str, list[tuple[str, float]]] = {}
    for r in rows:
        if r["kind"] != "cand":
            continue
        qid = r["qid"]
        ruby = r["samples"][0]["answer"] if r["samples"] else ""
        gold = (ROOT / "data" / "corpus" / "generated_dsl" / f"{qid}.rb").read_text()
        # strip fence
        from pipeline_common import strip_code_fence
        m = score_pair(strip_code_fence(ruby), gold)
        per_q_cands.setdefault(qid, []).append((r["cand_stem"], float(m["struct_mean"])))

    struct_means = []
    oracle_gap = []
    oracle_hits = 0
    for qid, cand_scores in per_q_cands.items():
        oracle_best = max(cand_scores, key=lambda x: x[1])[1]
        entries = dH.get(qid, {})
        def _s(cand_stem: str) -> float:
            v = entries.get(cand_stem)
            if v is None:
                return -1e9
            f = v["dH_full"] if v["dH_full"] is not None else 0.0
            s = v["dH_sel"]  if v["dH_sel"]  is not None else 0.0
            return alpha * f + (1 - alpha) * s
        chosen_stem, _ = max(cand_scores, key=lambda x: _s(x[0]))
        chosen_score = dict(cand_scores)[chosen_stem]
        struct_means.append(chosen_score)
        oracle_gap.append(oracle_best - chosen_score)
        oracle_hits += (chosen_score == oracle_best)
    return {
        "variant": variant, "alpha": alpha,
        "mean_struct": sum(struct_means) / len(struct_means),
        "mean_gap_to_oracle": sum(oracle_gap) / len(oracle_gap),
        "oracle_hit_rate": oracle_hits / len(struct_means),
        "n": len(struct_means),
    }


def eval_reranker_top1(rows: list[dict]) -> dict:
    struct_means, oracle_gap, hits = [], [], 0
    from pipeline_common import strip_code_fence
    per_q: dict[str, list[tuple[int, str, float]]] = {}
    for r in rows:
        if r["kind"] != "cand":
            continue
        qid = r["qid"]
        ruby = r["samples"][0]["answer"] if r["samples"] else ""
        gold = (ROOT / "data" / "corpus" / "generated_dsl" / f"{qid}.rb").read_text()
        m = score_pair(strip_code_fence(ruby), gold)
        per_q.setdefault(qid, []).append((r["cand_rank"], r["cand_stem"], float(m["struct_mean"])))
    for qid, xs in per_q.items():
        xs.sort()
        chosen = xs[0][2]
        oracle_best = max(x[2] for x in xs)
        struct_means.append(chosen)
        oracle_gap.append(oracle_best - chosen)
        hits += (chosen == oracle_best)
    return {"variant": "reranker_top1", "alpha": None,
            "mean_struct": sum(struct_means) / len(struct_means),
            "mean_gap_to_oracle": sum(oracle_gap) / len(oracle_gap),
            "oracle_hit_rate": hits / len(struct_means),
            "n": len(struct_means)}


def main() -> int:
    rows = collect_raw()
    results = [eval_reranker_top1(rows)]
    for v in VARIANTS:
        for a in (0.0, 0.5, 1.0):
            results.append(eval_variant(rows, v, a))

    # sort by mean_struct desc
    results.sort(key=lambda r: -r["mean_struct"])
    print(f"\n{'variant':<20} {'a':>4}  {'struct_mean':>11}  {'gap_to_oracle':>13}  {'oracle_hit':>10}  n")
    print("-" * 75)
    for r in results:
        a_s = f"{r['alpha']:.1f}" if r["alpha"] is not None else "  --"
        print(f"{r['variant']:<20} {a_s:>4}  {r['mean_struct']:>11.4f}  "
              f"{r['mean_gap_to_oracle']:>13.4f}  {r['oracle_hit_rate']:>10.2f}  {r['n']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
