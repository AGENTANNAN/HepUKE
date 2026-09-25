"""plot_calibration.py — write-back calibration in confidence (c) space.

Standalone (no source modification). Draws the write-back gate's calibration
with x = c = σ(50·(αΔH_full+(1−α)ΔH_sel)), y = mean struct_mean, over the 267
entropy-top1 write-back candidates. Overlays τ_write=0.5 and τ_high=0.80,
shading the admitted window [0.5, 0.80]. Left panel = quality-vs-c (quantile
bins); right panel = confidence-mass histogram.
"""
from __future__ import annotations
import json, sys, math
from collections import defaultdict
from pathlib import Path
import statistics as st

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _prereq import require  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from evaluator import score_pair
from pipeline_common import strip_code_fence
from splits import load_splits
from kb_evolution import confidence_from_dH


def boot_ci(vals, B=2000, alpha=0.05, seed=0):
    arr = np.asarray(vals, dtype=float)
    if len(arr) < 2:
        m = float(arr.mean()) if len(arr) else float("nan")
        return m, m
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(arr), size=(B, len(arr)))
    means = arr[idx].mean(axis=1)
    return float(np.quantile(means, alpha / 2)), float(np.quantile(means, 1 - alpha / 2))

GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"
OUT = ROOT / "token_conf" / "dsl_gate" / "results" / "kb_evolution"

ALPHA = 0.5
TAU_WRITE, TAU_HIGH = 0.5, 0.80

splits = load_splits()
query = list(splits["dev"]) + list(splits["test"])
cand = [json.loads(l) for l in require(GATE / "all.jsonl").open()]
by_qid = defaultdict(list)
for r in cand:
    by_qid[r["qid"]].append(r)

rows = []  # (c, struct_mean)
for qid in query:
    cs = by_qid.get(qid)
    if not cs:
        continue
    top1 = max(cs, key=lambda c: ALPHA * (c.get("dH_full") or 0) + (1 - ALPHA) * (c.get("dH_sel") or 0))
    ruby = strip_code_fence(top1.get("gen_ruby_first", "") or "")
    g = GOLD / f"{qid}.rb"
    if not g.exists():
        continue
    sm = float(score_pair(ruby, g.read_text())["struct_mean"])
    c = confidence_from_dH(top1.get("dH_full"), top1.get("dH_sel"), ALPHA)
    rows.append((c, sm))

rows.sort()
cs = [r[0] for r in rows]
sms = [r[1] for r in rows]
print(f"n = {len(rows)} write-back candidates")
print(f"c:   min={cs[0]:.4f}  median={st.median(cs):.4f}  max={cs[-1]:.4f}")
print(f"     p10={cs[int(len(cs)*0.10)]:.4f}  p90={cs[int(len(cs)*0.90)]:.4f}")

# quantile bins (equal-count)
NB = 10
edges = [cs[int(i * len(cs) / NB)] for i in range(NB)] + [cs[-1] + 1e-9]
bins = []
for i in range(NB):
    lo, hi = edges[i], edges[i + 1]
    inb = [(c, s) for (c, s) in rows if lo <= c < hi or (i == NB - 1 and lo <= c <= hi)]
    grp = [s for (c, s) in inb]
    ctr = st.mean([c for (c, s) in inb]) if inb else (lo + hi) / 2
    mean = st.mean(grp) if grp else float("nan")
    ci_lo, ci_hi = boot_ci(grp, B=2000, seed=1234 + i)
    fp = sum(1 for s in grp if s < 0.5)
    bins.append((ctr, lo, hi, mean, len(grp), fp, ci_lo, ci_hi))

print("\nquantile bins (x=c bin center, struct_mean, n, fp<0.5, 95% CI):")
for ctr, lo, hi, m, n, fp, cl, ch in bins:
    print(f"  c∈[{lo:.3f},{hi:.3f})  ctr={ctr:.3f}  struct={m:.4f}  "
          f"CI=[{cl:.3f},{ch:.3f}]  n={n:3d}  fp={fp}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

# left: quality vs c (with bootstrap 95% CI + n-based alpha shading)
xs = [b[0] for b in bins]
ys = [b[3] for b in bins]
ns = [b[4] for b in bins]
lo_err = [b[3] - b[6] for b in bins]
hi_err = [b[7] - b[3] for b in bins]

ax[0].axvspan(TAU_WRITE, TAU_HIGH, color="tab:green", alpha=0.10)
ax[0].axvline(TAU_WRITE, color="tab:red", ls="--", lw=1.2, label=r"$\tau_{\mathrm{write}}=0.5$")
ax[0].axvline(TAU_HIGH, color="tab:red", ls="--", lw=1.2, label=r"$\tau_{\mathrm{high}}=0.80$")
ax[0].plot(xs, ys, "-", color="tab:blue", lw=1.2, alpha=0.6, zorder=1)
for x, y, n, le, he in zip(xs, ys, ns, lo_err, hi_err):
    a = min(1.0, max(0.25, n / 30.0))
    color = "tab:blue" if n >= 20 else "gray"
    ax[0].errorbar([x], [y], yerr=[[le], [he]], fmt="o", color=color,
                   ecolor=color, alpha=a, capsize=3, ms=5, zorder=3)
    ax[0].annotate(f"n={n}", (x, y), textcoords="offset points", xytext=(0, 8),
                   ha="center", fontsize=7,
                   color="black" if n >= 20 else "gray")
ax[0].set_xlabel(r"confidence $c=\sigma(50(\alpha\Delta H_{\mathrm{full}}+(1-\alpha)\Delta H_{\mathrm{sel}}))$, $\alpha=0.5$")
ax[0].set_ylabel("mean struct_mean (bootstrap 95% CI)")
ax[0].set_title("Write-back calibration (n<20 in gray)")
ax[0].set_ylim(0.0, 1.05)
ax[0].legend(fontsize=8, loc="lower left")
ax[0].grid(alpha=0.25)

# right: confidence-mass histogram (equal-width)
ax[1].hist(cs, bins=30, range=(0.45, 1.0), color="tab:blue", alpha=0.6, edgecolor="white")
ax[1].axvspan(TAU_WRITE, TAU_HIGH, color="tab:green", alpha=0.10)
ax[1].axvline(TAU_WRITE, color="tab:red", ls="--", lw=1.2)
ax[1].axvline(TAU_HIGH, color="tab:red", ls="--", lw=1.2)
ax[1].set_xlabel(r"confidence $c$")
ax[1].set_ylabel("count")
ax[1].set_title("Confidence mass (admitted window shaded)")
ax[1].grid(alpha=0.25)

plt.tight_layout()
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    plt.savefig(OUT / f"plot_calibration.{ext}", dpi=150, bbox_inches="tight")
print(f"\n[+] saved {OUT/'plot_calibration.png'} / .pdf")
