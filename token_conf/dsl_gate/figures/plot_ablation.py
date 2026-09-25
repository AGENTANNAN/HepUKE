"""Ablation figures over the merged dev+test split.

Reads from runs_gate/all_{summary,perqid}.jsonl + all.jsonl +
all_raw.jsonl.gz. Writes to results/figures/:

  1. fig_alpha_sensitivity.png   struct_mean vs alpha (all/dev/test)
  2. fig_module_dH_violin.png    per-module ΔH violin over all candidate rows
  3. fig_oracle_recovery.png     bar chart of norm-% oracle recovery
  4. fig_calibration.png         dH bin -> mean struct_mean (per bin)

  5. ablation_module_dH.csv      module -> mean|median|std of ΔH across cands
  6. ablation_calibration.csv    (bin, mean_dH_sel, mean_struct_mean, n)
  7. ablation_cost.csv           strategy -> wall_ms mean, calls_per_query, tokens

CLI:
    python token_conf/dsl_gate/plot_ablation.py
"""
from __future__ import annotations

import csv
import gzip
import json
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "token_conf"))

from compute_entropy import tag_tokens, token_entropy   # noqa: E402
from segment_modules import tag_modules, MODULES        # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _prereq import require                             # noqa: E402


GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results"
FIGS = RESULTS / "figures"
FIGS.mkdir(parents=True, exist_ok=True)


# --------------- helpers ---------------

def _load_jsonl(p: Path) -> list[dict]:
    return [json.loads(l) for l in p.open()]


def _load_summary(split: str) -> dict[str, dict]:
    path = require(RESULTS / f"{split}_summary.jsonl",
                   f"python token_conf/dsl_gate/strategy_matrix.py --split {split} "
                   f"--alpha 0 0.3 0.5 0.7 1.0 --per-qid-out token_conf/dsl_gate/results/")
    return {r["strategy"]: r for r in _load_jsonl(path)}


# --------------- 1. alpha sensitivity fig ---------------


def fig_alpha_sensitivity() -> None:
    alphas = [0.0, 0.3, 0.5, 0.7, 1.0]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    colors = {"all": "#8e44ad", "dev": "#2980b9", "test": "#c0392b"}
    for split in ("all", "dev", "test"):
        d = _load_summary(split)
        ys = [d[f"entropy_top1(a={a:g})"]["struct_mean"] for a in alphas]
        n = int(d["oracle(struct_mean)"]["n_eval"])
        ax.plot(alphas, ys, "-o", color=colors[split], label=f"{split} (n={n})", lw=1.8, ms=6)
        # reranker baseline as horizontal dashed
        rr = d["reranker_top1"]["struct_mean"]
        ax.axhline(rr, ls="--", color=colors[split], alpha=0.35, lw=1)
        oracle = d["oracle(struct_mean)"]["struct_mean"]
        ax.axhline(oracle, ls=":", color=colors[split], alpha=0.35, lw=1)
    ax.set_xlabel(r"$\alpha$  (weight on $\Delta H_{\mathrm{full}}$; $1-\alpha$ on $\Delta H_{\mathrm{sel}}$)")
    ax.set_ylabel("struct_mean")
    ax.set_title("Entropy gating: $\\alpha$ sensitivity\n(dashed = reranker top-1; dotted = oracle)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_alpha_sensitivity.png", dpi=150)
    fig.savefig(FIGS / "fig_alpha_sensitivity.pdf")
    plt.close(fig)
    print(f"[+] {FIGS / 'fig_alpha_sensitivity.png'}")


# --------------- 2. per-module ΔH violin ---------------


def compute_per_module_dH() -> tuple[dict[str, list[float]], list[dict]]:
    """For each (qid, cand) sample-average per-module H; compare to baseline
    sample-average per-module H → per-module ΔH. Returns:
      per_module -> list of ΔH values (one per (qid, cand))
      raw_rows_list (for csv)
    """
    print("  loading all_raw.jsonl.gz ...")
    # per (qid, kind, cand_stem) -> mean_per_module_H
    per_row_H: dict[tuple[str, str, str | None], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    with gzip.open(GATE / "all_raw.jsonl.gz", "rt") as f:
        for i, line in enumerate(f):
            r = json.loads(line)
            toks = [dict(t) for t in (r.get("content_logprobs") or [])]
            if not toks:
                continue
            tag_tokens(toks)
            tag_modules(toks)
            # per-module H mean for THIS sample
            mod_H: dict[str, list[float]] = defaultdict(list)
            for t in toks:
                if not t.get("kept"):
                    continue
                H = token_entropy(t)["H_tail_nats"]
                mod_H[t.get("module") or "preamble"].append(H)
            key = (r["qid"], r["kind"], r.get("cand_stem"))
            for m, xs in mod_H.items():
                per_row_H[key][m].append(sum(xs) / len(xs))
            if (i + 1) % 500 == 0:
                print(f"    processed {i+1} raw rows")
    # sample-average -> single value per module per row
    per_row_mean: dict[tuple[str, str, str | None], dict[str, float]] = {}
    for key, mod_map in per_row_H.items():
        per_row_mean[key] = {m: (sum(vs) / len(vs)) if vs else None
                             for m, vs in mod_map.items()}

    # for each qid, baseline vs each cand: ΔH per module
    per_mod_dH: dict[str, list[float]] = {m: [] for m in MODULES}
    per_qid_bases: dict[str, dict[str, float]] = {}
    for (qid, kind, _), mm in per_row_mean.items():
        if kind == "baseline":
            per_qid_bases[qid] = mm
    rows_out: list[dict] = []
    for (qid, kind, cand_stem), mm in per_row_mean.items():
        if kind != "cand":
            continue
        base = per_qid_bases.get(qid)
        if not base:
            continue
        for m in MODULES:
            bv, cv = base.get(m), mm.get(m)
            if bv is None or cv is None:
                continue
            dH = bv - cv
            per_mod_dH[m].append(dH)
            rows_out.append({"qid": qid, "cand": cand_stem, "module": m, "dH": dH})
    return per_mod_dH, rows_out


def fig_module_violin(per_mod_dH: dict[str, list[float]]) -> None:
    mods_no_note = [m for m in MODULES if m != "note"]
    data = [per_mod_dH.get(m, []) for m in mods_no_note]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    parts = ax.violinplot(data, showmedians=True, widths=0.85)
    for pc, m in zip(parts["bodies"], mods_no_note):
        pc.set_facecolor("#8e44ad" if m == "selection" else "#7f8c8d")
        pc.set_alpha(0.55 if m == "selection" else 0.35)
        pc.set_edgecolor("black")
    ax.set_xticks(range(1, len(mods_no_note) + 1))
    ax.set_xticklabels(mods_no_note, rotation=25, ha="right")
    ax.axhline(0, ls="--", color="black", alpha=0.4)
    ax.set_ylabel(r"$\Delta H = H_{\mathrm{base}} - H_{\mathrm{cand}}$  (nats)")
    ax.set_title(f"Per-module ΔH across all (qid, candidate) pairs (n={sum(len(v) for v in per_mod_dH.values())} obs)\n"
                 "Positive = reference helped the module get more confident. `selection` is the widest-spread signal.")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_module_dH_violin.png", dpi=150)
    fig.savefig(FIGS / "fig_module_dH_violin.pdf")
    plt.close(fig)
    print(f"[+] {FIGS / 'fig_module_dH_violin.png'}")


def csv_module_dH(per_mod_dH: dict[str, list[float]]) -> None:
    with (RESULTS / "ablation_module_dH.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["module", "n", "mean", "median", "std", "abs_mean", "pct_positive"])
        for m in MODULES:
            xs = per_mod_dH.get(m, [])
            if not xs:
                w.writerow([m, 0, "", "", "", "", ""])
                continue
            w.writerow([m, len(xs), f"{st.mean(xs):.5f}", f"{st.median(xs):.5f}",
                        f"{st.pstdev(xs):.5f}", f"{sum(abs(x) for x in xs)/len(xs):.5f}",
                        f"{sum(1 for x in xs if x > 0)/len(xs):.3f}"])
    print(f"[+] {RESULTS / 'ablation_module_dH.csv'}")


# --------------- 3. oracle recovery bar ---------------


def fig_oracle_recovery() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    splits = ("all", "dev", "test")
    strategies = ["no_rag", "reranker_top1", "random",
                  "entropy_top1(a=0)", "entropy_top1(a=0.5)", "entropy_top1(a=1)",
                  "oracle(struct_mean)"]
    labels = ["no-RAG", "reranker top-1", "random",
              "entropy α=0", "entropy α=0.5", "entropy α=1", "oracle"]
    x = np.arange(len(strategies))
    width = 0.26
    colors = {"all": "#8e44ad", "dev": "#2980b9", "test": "#c0392b"}
    for i, split in enumerate(splits):
        d = _load_summary(split)
        no_rag = d["no_rag"]["struct_mean"]
        oracle = d["oracle(struct_mean)"]["struct_mean"]
        gap = oracle - no_rag
        norm = []
        for k in strategies:
            v = d[k]["struct_mean"]
            norm.append((v - no_rag) / gap * 100 if gap > 0 else 0)
        ax.bar(x + (i - 1) * width, norm, width=width, color=colors[split],
               label=f"{split} (n={int(d['no_rag']['n_eval'])})", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Oracle-gap recovery (%)")
    ax.set_title("Fraction of (oracle − no-RAG) headroom recovered by each strategy")
    ax.axhline(0, color="black", lw=0.5)
    ax.axhline(100, ls=":", color="black", lw=0.5)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_oracle_recovery.png", dpi=150)
    fig.savefig(FIGS / "fig_oracle_recovery.pdf")
    plt.close(fig)
    print(f"[+] {FIGS / 'fig_oracle_recovery.png'}")


# --------------- 4. calibration fig ---------------


def _score_ruby(ruby: str, gold: str) -> float:
    """Ad-hoc struct_mean via evaluator (cached module import)."""
    from evaluator import score_pair
    if not ruby or not gold:
        return 0.0
    return float(score_pair(ruby, gold)["struct_mean"])


def calibration() -> list[dict]:
    """For each candidate row: bin its dH_sel, score its gen_ruby against gold,
    and average struct_mean per bin."""
    sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
    from pipeline_common import strip_code_fence
    cand_rows = _load_jsonl(GATE / "all.jsonl")
    # bin edges (nats). ΔH_sel range roughly [-0.05, +0.08] from earlier peek.
    edges = [-0.05, -0.01, 0.0, 0.005, 0.01, 0.02, 0.04, 0.08]
    labels = [f"[{edges[i]:+.3g},{edges[i+1]:+.3g})" for i in range(len(edges)-1)]
    bins: dict[str, list[float]] = {l: [] for l in labels}
    for r in cand_rows:
        dHs = r.get("dH_sel")
        if dHs is None:
            continue
        gold_path = ROOT / "data" / "corpus" / "generated_dsl" / f"{r['qid']}.rb"
        if not gold_path.exists():
            continue
        gold = gold_path.read_text()
        ruby = strip_code_fence(r.get("gen_ruby_first", "") or "")
        s = _score_ruby(ruby, gold)
        # find bin
        placed = False
        for i in range(len(edges) - 1):
            if edges[i] <= dHs < edges[i+1]:
                bins[labels[i]].append(s)
                placed = True
                break
        if not placed:
            # out-of-range → left/right tail
            if dHs < edges[0]:
                bins[labels[0]].append(s)
            else:
                bins[labels[-1]].append(s)
    rows_out = []
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    xs = np.arange(len(labels))
    means = []
    ns = []
    for i, lb in enumerate(labels):
        xs_bin = bins[lb]
        if xs_bin:
            m = st.mean(xs_bin)
            se = st.pstdev(xs_bin) / (len(xs_bin) ** 0.5) if len(xs_bin) > 1 else 0
        else:
            m = 0.0
            se = 0.0
        means.append(m)
        ns.append(len(xs_bin))
        rows_out.append({"bin": lb, "n": len(xs_bin), "mean_struct": m, "sem": se})
        ax.errorbar([i], [m], yerr=[se], fmt="o", color="#8e44ad", ms=6, capsize=3, lw=1.2)
    ax.plot(xs, means, "-", color="#8e44ad", alpha=0.5, lw=1.2)
    for i, n in enumerate(ns):
        ax.annotate(f"n={n}", (i, means[i]), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=7, color="grey")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_xlabel(r"$\Delta H_{\mathrm{sel}}$ bin (nats)")
    ax.set_ylabel("mean struct_mean of candidate's ruby")
    ax.set_title("Calibration: does larger $\\Delta H_{\\mathrm{sel}}$ imply a better candidate?\n"
                 "Monotone-rising = ΔH is a good proxy for downstream match")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_calibration.png", dpi=150)
    fig.savefig(FIGS / "fig_calibration.pdf")
    plt.close(fig)
    print(f"[+] {FIGS / 'fig_calibration.png'}")
    # csv
    with (RESULTS / "ablation_calibration.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["bin", "n", "mean_struct", "sem"])
        for r in rows_out:
            w.writerow([r["bin"], r["n"], f"{r['mean_struct']:.4f}", f"{r['sem']:.4f}"])
    print(f"[+] {RESULTS / 'ablation_calibration.csv'}")
    return rows_out


# --------------- 5. cost table ---------------


def cost_table() -> None:
    cand_rows = _load_jsonl(GATE / "all.jsonl")
    base_rows = _load_jsonl(GATE / "all_baseline.jsonl")
    n_qid = len(base_rows)
    wall_base_ms = [b["wall_ms_base"] for b in base_rows if b.get("wall_ms_base")]
    wall_cand_ms = [r["wall_ms_cand"] for r in cand_rows if r.get("wall_ms_cand")]
    n_tok_base = [b.get("n_tok_base_mean") for b in base_rows if b.get("n_tok_base_mean")]
    n_tok_cand = [r.get("n_tok_cand_mean") for r in cand_rows if r.get("n_tok_cand_mean")]
    # strategies
    # no_rag: only baseline calls (3 per query)
    # reranker/entropy/random/oracle: baseline + 5*3 cand = 18 calls per query
    rows = []
    def _row(name, calls_per_q, wall_ms_per_q, tok_per_q):
        rows.append({"strategy": name, "calls_per_query": calls_per_q,
                     "wall_ms_per_query": wall_ms_per_q, "tokens_per_query": tok_per_q})
    # per-q wall = baseline_wall + 5*mean_cand_wall (approx; qid-wise below)
    per_qid_wall: dict[str, list[float]] = defaultdict(list)
    per_qid_tok:  dict[str, list[float]] = defaultdict(list)
    for r in cand_rows:
        per_qid_wall[r["qid"]].append(r.get("wall_ms_cand") or 0)
        per_qid_tok[r["qid"]].append(r.get("n_tok_cand_mean") or 0)
    q_walls = []
    q_toks = []
    for b in base_rows:
        q = b["qid"]
        w = (b.get("wall_ms_base") or 0) + sum(per_qid_wall.get(q, []))
        t = 3 * (b.get("n_tok_base_mean") or 0) + 3 * sum(per_qid_tok.get(q, []))
        q_walls.append(w)
        q_toks.append(t)
    _row("no_rag", 3, st.mean(wall_base_ms) if wall_base_ms else 0,
         3 * st.mean(n_tok_base) if n_tok_base else 0)
    _row("+entropy/reranker/random/oracle (share cache)", 18,
         st.mean(q_walls) if q_walls else 0,
         st.mean(q_toks) if q_toks else 0)
    print("\nstrategy                                       calls/q  wall_ms/q  tokens/q")
    for r in rows:
        print(f"  {r['strategy']:<45}  {r['calls_per_query']:>7}  {r['wall_ms_per_query']:>9.0f}  {r['tokens_per_query']:>8.0f}")
    with (RESULTS / "ablation_cost.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["strategy", "calls_per_query", "wall_ms_per_query", "tokens_per_query"])
        for r in rows:
            w.writerow([r["strategy"], r["calls_per_query"],
                        f"{r['wall_ms_per_query']:.0f}",
                        f"{r['tokens_per_query']:.0f}"])
    print(f"[+] {RESULTS / 'ablation_cost.csv'}")


def main() -> int:
    print("=== ablation figures ===")
    print("[1/5] alpha sensitivity fig ...")
    fig_alpha_sensitivity()
    print("[2/5] per-module ΔH (heavy: loads all_raw.jsonl.gz) ...")
    per_mod_dH, _ = compute_per_module_dH()
    fig_module_violin(per_mod_dH)
    csv_module_dH(per_mod_dH)
    print("[3/5] oracle recovery bar ...")
    fig_oracle_recovery()
    print("[4/5] calibration ...")
    calibration()
    print("[5/5] cost table ...")
    cost_table()
    print("\n=== done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
