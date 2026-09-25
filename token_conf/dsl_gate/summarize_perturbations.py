"""Struct + Strict summary of the perturbation runs.

Offline — uses evaluator.score_pair only, no LLM judge. Prints No-RAG /
relevance / entropy struct_mean and strict_joint per variant, plus the
per-variant oracle (max struct over that target's candidates).

Usage::
    python token_conf/dsl_gate/summarize_perturbations.py
    python token_conf/dsl_gate/summarize_perturbations.py --splits pert20
"""
from __future__ import annotations

import argparse
import json, sys, statistics as st
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
sys.path.insert(0, str(ROOT / "token_conf"))

from evaluator import score_pair
from pipeline_common import strip_code_fence
from variants import tags

GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"


def ruby(c): return strip_code_fence(c.get("gen_ruby_first", "") or "")
def ent(c):  return 0.5 * (c.get("dH_full") or 0.0) + 0.5 * (c.get("dH_sel") or 0.0)


def summarize(prefix: str, variants: list[tuple[str, str]]) -> None:
    print(f"\n### {prefix}")
    print("| variant | n | Strict NR/Rel/Ent | Struct NR/Rel/Ent | Oracle Struct | Recovery Rel%/Ent% (per-var) |")
    print("|---|---|---|---|---|---|")
    for label, tag in variants:
        cand_p = GATE / f"{tag}.jsonl"
        base_p = GATE / f"{tag}_baseline.jsonl"
        if not cand_p.exists() or not base_p.exists():
            print(f"| {label} | — | missing | | | |")
            continue
        by_qid = defaultdict(list)
        for l in cand_p.open():
            r = json.loads(l); by_qid[r["qid"]].append(r)
        base = {json.loads(l)["qid"]: json.loads(l) for l in base_p.open()}
        qids = sorted(by_qid.keys() & base.keys())
        nr_s, rel_s, ent_s, orc_s = [], [], [], []
        nr_x, rel_x, ent_x = [], [], []
        for q in qids:
            g = (GOLD / f"{q}.rb").read_text()
            m_nr = score_pair(ruby(base[q]), g)
            nr_s.append(m_nr["struct_mean"])
            nr_x.append(0.5 * m_nr["dataset_em"] + 0.5 * m_nr["decay_card_em"])
            cs = by_qid[q]
            top_rel = min(cs, key=lambda c: c.get("cand_rank", 999))
            top_ent = max(cs, key=ent)
            m_rel = score_pair(ruby(top_rel), g)
            m_ent = score_pair(ruby(top_ent), g)
            rel_s.append(m_rel["struct_mean"]); ent_s.append(m_ent["struct_mean"])
            rel_x.append(0.5 * m_rel["dataset_em"] + 0.5 * m_rel["decay_card_em"])
            ent_x.append(0.5 * m_ent["dataset_em"] + 0.5 * m_ent["decay_card_em"])
            # Oracle = max struct over all candidates
            orc_s.append(max(score_pair(ruby(c), g)["struct_mean"] for c in cs))
        n = len(qids)
        NR, R, E, O = st.mean(nr_s), st.mean(rel_s), st.mean(ent_s), st.mean(orc_s)
        SNR, SR, SE = st.mean(nr_x), st.mean(rel_x), st.mean(ent_x)
        gap = O - NR
        rec_r = (R - NR) / gap * 100 if gap > 0 else float("nan")
        rec_e = (E - NR) / gap * 100 if gap > 0 else float("nan")
        print(f"| {label} | {n} | {SNR:.3f}/{SR:.3f}/{SE:.3f} | "
              f"{NR:.3f}/{R:.3f}/{E:.3f} | {O:.3f} | {rec_r:+5.1f}%/{rec_e:+5.1f}% |")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--splits", nargs="+", default=["pert20", "pert50"],
                    help="output tag prefixes to summarize (default: pert20 pert50)")
    args = ap.parse_args()
    print("# Perturbation summary")
    for prefix in args.splits:
        summarize(prefix, [(label, t_) for _, label, t_ in tags(prefix)])


if __name__ == "__main__":
    main()
