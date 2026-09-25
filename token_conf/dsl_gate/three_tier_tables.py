"""Assemble the three-tier (strict / structural / Phys-Acc) tables.

Two outputs, both in the paper's 4-column layout
(Strict Match | Structural Match | Phys-Acc | Match Recovery):

  1. GATE table (267 held-out targets, all.jsonl / all_baseline.jsonl):
     strategies as rows.
  2. PERTURBATION table (20 dev targets, pert_* jsonl): variants as rows.

Strict Match  = strict_joint = 0.5*dataset_em + 0.5*decay_card_em
                (values-must-match tier: exact IDs + exact decay cards)
Structural    = struct_mean = 0.5*method-F1 + 0.5*slot-F1  (the paper's primary)
Phys-Acc      = Phys-Acc-100 / 100  (7-dim weighted LLM judge, judge_phys_acc / judge_perturbations)
Match Recovery= oracle-gap recovery on struct_mean, normalized to [No-RAG=0, Oracle=100].

Offline for strict/structural; Phys-Acc read from the judge_phys_acc/judge_perturbations judgment caches.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _prereq import require  # noqa: E402

from evaluator import score_pair               # noqa: E402
from pipeline_common import strip_code_fence      # noqa: E402
from variants import tags                        # noqa: E402

GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"
PA100 = ROOT / "token_conf" / "dsl_gate" / "results" / "phys_acc100"
PA100_PERT = PA100 / "perturb_judgments.jsonl"
PA100_MAIN = PA100 / "judgments.jsonl"

# (display label, judge variant key, runs_gate tag)
PERT_PREFIX = "pert20"
VARIANTS = [(label, t_, t_) for _, label, t_ in tags(PERT_PREFIX)]


def ruby(c: dict) -> str:
    return strip_code_fence(c.get("gen_ruby_first", "") or "")


def ent_score(c: dict) -> float:
    return 0.5 * (c.get("dH_full") or 0.0) + 0.5 * (c.get("dH_sel") or 0.0)


def load_phys100_main() -> dict[tuple, dict]:
    out: dict[tuple, dict] = {}
    if PA100_MAIN.exists():
        for l in PA100_MAIN.open():
            l = l.strip()
            if not l:
                continue
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                continue
            if r.get("total") is not None and not r.get("alt_judge"):
                out[(r["qid"], r["kind"], r.get("cand_stem"))] = r
    return out


def load_phys100_pert() -> dict[tuple, dict]:
    out: dict[tuple, dict] = {}
    if PA100_PERT.exists():
        for l in PA100_PERT.open():
            l = l.strip()
            if not l:
                continue
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                continue
            if r.get("total") is not None:
                out[(r["variant"], r["qid"], r["kind"], r.get("cand_stem"))] = r
    return out


def strategy_ruby(strat: str, qid: str, cands: list[dict], base: dict,
                  phys: dict | None = None) -> str:
    if strat == "no_rag":
        return ruby(base[qid])
    if strat == "rel":
        return ruby(min(cands, key=lambda c: c.get("cand_rank", 999)))
    if strat.startswith("ent"):
        a = float(strat.split("=")[1])
        return ruby(max(cands, key=lambda c: a * (c.get("dH_full") or 0.0)
                        + (1 - a) * (c.get("dH_sel") or 0.0)))
    if strat == "oracle_struct":
        def sm(c):
            return score_pair(ruby(c), (GOLD / f"{qid}.rb").read_text())["struct_mean"]
        return ruby(max(cands, key=sm))
    if strat == "oracle_phys":
        def pa(c):
            r = (phys or {}).get((qid, "cand", c.get("cand_stem")))
            return r["total"] if r else -1
        return ruby(max(cands, key=pa))
    raise ValueError(strat)


def pick_stem(strat: str, qid: str, cands: list[dict],
              phys: dict | None = None) -> str | None:
    """The cand_stem a strategy selects (for Phys-Acc lookup)."""
    if strat == "no_rag":
        return None
    if strat == "rel":
        return min(cands, key=lambda c: c.get("cand_rank", 999)).get("cand_stem")
    if strat.startswith("ent"):
        a = float(strat.split("=")[1])
        return max(cands, key=lambda c: a * (c.get("dH_full") or 0.0)
                   + (1 - a) * (c.get("dH_sel") or 0.0)).get("cand_stem")
    if strat == "oracle_struct":
        def sm(c):
            return score_pair(ruby(c), (GOLD / f"{qid}.rb").read_text())["struct_mean"]
        return max(cands, key=sm).get("cand_stem")
    if strat == "oracle_phys":
        def pa(c):
            r = (phys or {}).get((qid, "cand", c.get("cand_stem")))
            return r["total"] if r else -1
        return max(cands, key=pa).get("cand_stem")
    raise ValueError(strat)


def gate_table() -> None:
    by_qid: dict[str, list[dict]] = defaultdict(list)
    base: dict[str, dict] = {}
    for l in require(GATE / "all.jsonl").open():
        r = json.loads(l)
        by_qid[r["qid"]].append(r)
    for l in require(GATE / "all_baseline.jsonl").open():
        b = json.loads(l)
        base[b["qid"]] = b
    qids = sorted({q for q in by_qid if (GOLD / f"{q}.rb").exists()})
    phys = load_phys100_main()

    strats = ["no_rag", "rel", "ent=0", "ent=0.5", "ent=1",
              "oracle_struct", "oracle_phys"]
    labels = {"no_rag": "No-RAG", "rel": "Relevance Top-1",
              "ent=0": r"Confidence Top-1 ($\alpha=0.0$)",
              "ent=0.5": r"Confidence Top-1 ($\alpha=0.5$)",
              "ent=1": r"Confidence Top-1 ($\alpha=1.0$)",
              "oracle_struct": "Oracle (struct upper)",
              "oracle_phys": "Oracle (Phys-Acc upper)"}

    rows = {}
    for st in strats:
        strict, struct, pa = [], [], []
        for q in qids:
            m = score_pair(strategy_ruby(st, q, by_qid[q], base, phys),
                           (GOLD / f"{q}.rb").read_text())
            # NEW: strict_joint = 0.5*dataset_em + 0.5*decay_card_em
            strict.append(0.5 * m["dataset_em"] + 0.5 * m["decay_card_em"])
            struct.append(m["struct_mean"])
            stem = pick_stem(st, q, by_qid[q], phys)
            r = phys.get((q, "baseline" if st == "no_rag" else "cand", stem))
            if r is not None:
                pa.append(r["total"])
        rows[st] = {"strict": sum(strict) / len(strict),
                    "struct": sum(struct) / len(struct),
                    "phys": sum(pa) / len(pa) / 100.0 if pa else float("nan")}
    nr, orc = rows["no_rag"]["struct"], rows["oracle_struct"]["struct"]
    gap = orc - nr

    print("\n### GATE table (n=267)\n")
    print("| strategy | Strict Match | Structural Match | Phys-Acc | Match Recovery |")
    print("|---|---|---|---|---|")
    for st in strats:
        r = rows[st]
        rec = 0.0 if st == "no_rag" else (100.0 if st == "oracle_struct"
                                          else (r["struct"] - nr) / gap * 100)
        print(f"| {labels[st]} | {r['strict']:.3f} | {r['struct']:.3f} | "
              f"{r['phys']:.3f} | {rec:.1f}% |")


def perturb_table() -> None:
    phys = load_phys100_pert()
    print("\n### PERTURBATION table (n=20 per variant)\n")
    print("| variant | Strict (NR/Rel/Ent) | Struct (NR/Rel/Ent) | Phys-Acc (NR/Rel/Ent) |")
    print("|---|---|---|---|")
    for label, vkey, tag in VARIANTS:
        by_qid: dict[str, list[dict]] = defaultdict(list)
        base: dict[str, dict] = {}
        for l in (GATE / f"{tag}.jsonl").open():
            r = json.loads(l)
            by_qid[r["qid"]].append(r)
        for l in (GATE / f"{tag}_baseline.jsonl").open():
            b = json.loads(l)
            base[b["qid"]] = b
        qids = sorted(by_qid)
        agg = {"no_rag": {"strict": [], "struct": [], "pa": []},
               "rel": {"strict": [], "struct": [], "pa": []},
               "ent": {"strict": [], "struct": [], "pa": []}}
        for q in qids:
            g = (GOLD / f"{q}.rb").read_text()
            cands = by_qid[q]
            picks = {"no_rag": ruby(base[q]),
                     "rel": ruby(min(cands, key=lambda c: c.get("cand_rank", 999))),
                     "ent": ruby(max(cands, key=ent_score))}
            stems = {"rel": min(cands, key=lambda c: c.get("cand_rank", 999)).get("cand_stem"),
                     "ent": max(cands, key=ent_score).get("cand_stem")}
            for st in ("no_rag", "rel", "ent"):
                m = score_pair(picks[st], g)
                agg[st]["strict"].append(0.5 * m["dataset_em"] + 0.5 * m["decay_card_em"])
                agg[st]["struct"].append(m["struct_mean"])
                stem = None if st == "no_rag" else stems[st]
                r = phys.get((vkey, q, "baseline" if st == "no_rag" else "cand", stem))
                if r is not None:
                    agg[st]["pa"].append(r["total"])
        def f(st, k):
            v = agg[st][k]
            return sum(v) / len(v) if v else float("nan")
        def fmt(st, k):
            return f"{f(st, k):.3f}"
        print(f"| {label} | {fmt('no_rag','strict')}/{fmt('rel','strict')}/{fmt('ent','strict')} "
              f"| {fmt('no_rag','struct')}/{fmt('rel','struct')}/{fmt('ent','struct')} "
              f"| {f('no_rag','pa')/100:.2f}/{f('rel','pa')/100:.2f}/{f('ent','pa')/100:.2f} |")


if __name__ == "__main__":
    gate_table()
    perturb_table()
