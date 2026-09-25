"""Compute per-token predictive entropy from a *_out.json produced by
run_dsl_skill.py / get_token_logprobs.py, and emit it as the GOLD signal.

Each content token carries the model's top-K (=20) candidate distribution as
natural-log probabilities. Predictive entropy at position t is

    H_t = - Σ_i p_i · ln p_i        (nats;  /ln2 -> bits)

Two subtleties handled explicitly:

  * TRUNCATION. top-K is a truncated distribution; the tail carries
    m = 1 - Σ_{i<=K} p_i  of the mass. We report both:
      - H_renorm : renormalize the visible top-K to sum 1 (ignores the tail)
      - H_tail   : treat the whole tail as ONE extra bucket of mass m
                   (a lower bound on the true entropy; the true tail, being
                   spread over many tokens, has >= this much entropy)
    H_tail is the more honest per-token uncertainty and is used as the
    default gold value.

Output (JSON, --out):
  {
    "source": <input json path>,
    "model": ...,
    "n_tokens": N,
    "gold": "H_tail_nats",                 # which field is the gold scalar/seq
    "per_token": [ {i, token, logprob, H_renorm_nats, H_tail_nats,
                    H_tail_bits, tail_mass, top2_gap}, ... ],
    "summary": { mean/sum/min/max over H_tail_nats, ... }
  }

The full per-token H_tail sequence IS the gold distribution.

Usage:
    python token_conf/compute_entropy.py dsl_out_x3872.json
    python token_conf/compute_entropy.py dsl_out_x3872.json --out gold_x3872.json --plot gold_x3872.png
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
LN2 = math.log(2.0)

_HEREDOC_OPEN = re.compile(r"<<[~-]?[\"']?([A-Za-z_]\w*)")


def load_tokens(path: Path) -> tuple[list[dict], dict]:
    d = json.loads(path.read_text())
    toks = d.get("content_logprobs") or []
    if not toks:
        raise SystemExit(f"[!] {path} has no content_logprobs (max_tokens too small?)")
    return toks, d


def classify_chars(text: str) -> list[str]:
    """Char-level Ruby lexer: label each char 'c' (code), '#' (comment), or 'w'
    (whitespace). String and heredoc bodies count as CODE (they carry semantics
    like "==2" or the decay-card particle names); '#' is a comment ONLY in
    normal code state, never inside a string or heredoc.
    """
    out: list[str] = []
    state = "normal"          # normal | sstr | dstr | comment | heredoc
    esc = False               # previous char was a backslash (in a string)
    heredoc_term = ""         # active heredoc terminator (leading text stripped)
    pending_heredoc = ""      # terminator captured; body starts after next '\n'
    at_line_start = True      # cursor sits at first non-consumed char of a line

    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        is_ws = ch in " \t\r\n"

        # ---- heredoc body: consume whole lines until the terminator line ----
        if state == "heredoc":
            nl = text.find("\n", i)
            end = nl if nl != -1 else n
            line = text[i:end]
            if line.strip() == heredoc_term:
                state = "normal"          # terminator line is code
                for c in line:
                    out.append("w" if c in " \t\r" else "c")
            else:
                for c in line:
                    out.append("w" if c in " \t\r" else "c")
            if nl != -1:
                out.append("w")           # the newline
                i = end + 1
            else:
                i = end
            continue

        if state == "comment":
            out.append("w" if is_ws else "#")
            if ch == "\n":
                state = "normal"
                at_line_start = True
            i += 1
            continue

        if state in ("sstr", "dstr"):
            out.append("w" if is_ws else "c")
            quote = "'" if state == "sstr" else '"'
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == quote:
                state = "normal"
            i += 1
            continue

        # ---- normal code ----
        if ch == "\n":
            out.append("w")
            i += 1
            if pending_heredoc:
                state, heredoc_term, pending_heredoc = "heredoc", pending_heredoc, ""
            at_line_start = True
            continue
        if is_ws:
            out.append("w")
            i += 1
            continue
        if ch == "#":
            out.append("#")
            state = "comment"
            i += 1
            continue
        if ch == "'":
            out.append("c")
            state = "sstr"
            i += 1
            at_line_start = False
            continue
        if ch == '"':
            out.append("c")
            state = "dstr"
            i += 1
            at_line_start = False
            continue
        if text[i:i + 2] == "<<":
            m = _HEREDOC_OPEN.match(text, i)
            if m:
                pending_heredoc = m.group(1)
        out.append("c")
        at_line_start = False
        i += 1

    if len(out) != len(text):
        out = (out + ["c"] * len(text))[:len(text)]
    return out


def tag_tokens(toks: list[dict]) -> None:
    """Annotate each token in-place with 'kind' (code|comment|ws) and 'kept'.
    A token is KEPT iff it contains at least one code (non-ws, non-comment) char.
    """
    text = "".join(t["token"] for t in toks)
    cls = classify_chars(text)
    pos = 0
    for t in toks:
        s = t["token"]
        seg = cls[pos:pos + len(s)]
        pos += len(s)
        has_code = "c" in seg
        has_comment = "#" in seg
        if has_code:
            kind = "code"
        elif has_comment:
            kind = "comment"
        else:
            kind = "ws"
        t["kind"] = kind
        t["kept"] = has_code


def token_entropy(tok: dict) -> dict:
    """Entropy of one position from its truncated top-K logprob list."""
    tl = tok.get("top_logprobs") or []
    # visible probabilities (natural-log -> linear)
    ps = [math.exp(c["logprob"]) for c in tl]
    vis_mass = sum(ps)
    tail_mass = max(0.0, 1.0 - vis_mass)

    # (1) renormalized-over-visible entropy
    if vis_mass > 0:
        H_renorm = -sum((p / vis_mass) * math.log(p / vis_mass) for p in ps if p > 0)
    else:
        H_renorm = 0.0

    # (2) tail-as-single-bucket entropy (lower bound on true entropy)
    H_tail = -sum(p * math.log(p) for p in ps if p > 0)
    if tail_mass > 1e-12:
        H_tail += -tail_mass * math.log(tail_mass)

    gap = (tl[0]["logprob"] - tl[1]["logprob"]) if len(tl) >= 2 else float("inf")

    return {
        "H_renorm_nats": H_renorm,
        "H_tail_nats": H_tail,
        "H_tail_bits": H_tail / LN2,
        "tail_mass": tail_mass,
        "top2_gap": gap,
    }


def build_gold(toks: list[dict], meta: dict, src: Path, *, keep_all: bool = False) -> dict:
    tag_tokens(toks)  # annotate kind / kept
    per_token = []
    for i, t in enumerate(toks):
        e = token_entropy(t)
        per_token.append({
            "i": i,
            "token": t["token"],
            "logprob": t["logprob"],
            "kind": t["kind"],
            "kept": t["kept"],
            **e,
        })

    # GOLD is computed on CODE tokens only (comment + whitespace excluded),
    # unless keep_all is set. The full per-token list is retained either way.
    kept = per_token if keep_all else [p for p in per_token if p["kept"]]
    Hs = [p["H_tail_nats"] for p in kept]
    n = len(Hs)
    if n == 0:
        raise SystemExit("[!] no code tokens survived filtering")
    finite_gaps = [p["top2_gap"] for p in kept if math.isfinite(p["top2_gap"])]

    from collections import Counter
    kinds = Counter(p["kind"] for p in per_token)

    summary = {
        "n_tokens_total": len(per_token),
        "n_tokens_code": kinds["code"],
        "n_tokens_comment": kinds["comment"],
        "n_tokens_ws": kinds["ws"],
        "n_tokens_gold": n,                    # tokens the gold is computed over
        "H_tail_nats_mean": sum(Hs) / n,
        "H_tail_nats_sum": sum(Hs),
        "H_tail_nats_min": min(Hs),
        "H_tail_nats_max": max(Hs),
        "H_tail_bits_mean": (sum(Hs) / n) / LN2,
        "mean_tail_mass": sum(p["tail_mass"] for p in kept) / n,
        "mean_top2_gap": (sum(finite_gaps) / len(finite_gaps)) if finite_gaps else None,
    }
    return {
        "source": str(src),
        "model": meta.get("model"),
        "gold": "H_tail_nats",
        "gold_scope": "all_tokens" if keep_all else "code_tokens_only",
        "per_token": per_token,               # full sequence, each with kind/kept
        "summary": summary,
    }


def top_uncertain(gold: dict, n: int) -> None:
    # rank only over the gold scope (kept tokens) unless gold spans everything
    scope_all = gold["gold_scope"] == "all_tokens"
    pool = [p for p in gold["per_token"] if scope_all or p["kept"]]
    ranked = sorted(pool, key=lambda p: p["H_tail_nats"], reverse=True)[:n]
    print(f"\n# TOP {n} highest-entropy CODE tokens (semantic decision hotspots)")
    print(f"{'rank':>4} {'idx':>5} {'token':<20} {'H(nats)':>9} {'H(bits)':>9} {'tail':>7}   context")
    toks_by_i = {p["i"]: p for p in gold["per_token"]}
    m = len(gold["per_token"])
    for r, p in enumerate(ranked, 1):
        i = p["i"]
        lo, hi = max(0, i - 5), min(m, i + 6)
        ctx = "".join(
            (f"[{toks_by_i[k]['token']}]" if k == i else toks_by_i[k]["token"])
            for k in range(lo, hi)
        ).replace("\n", "⏎")
        print(f"{r:>4} {i:>5} {p['token']!r:<20} {p['H_tail_nats']:9.4f} "
              f"{p['H_tail_bits']:9.4f} {p['tail_mass']:7.3f}   …{ctx}…")


def plot(gold: dict, out_png: Path, mark_top_n: int = 12) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as _fm

    for _p in ("/usr/share/fonts/google-droid-sans-fonts/DroidSansFallbackFull.ttf",
               "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
               "/usr/share/fonts/truetype/arphic/uming.ttc"):
        if Path(_p).exists():
            _fm.fontManager.addfont(_p)
            plt.rcParams["font.sans-serif"] = ["DejaVu Sans", _fm.FontProperties(fname=_p).get_name()]
            plt.rcParams["axes.unicode_minus"] = False
            break

    scope_all = gold["gold_scope"] == "all_tokens"
    all_pts = gold["per_token"]
    pts = all_pts if scope_all else [p for p in all_pts if p["kept"]]
    xs = [p["i"] for p in pts]
    H = [p["H_tail_nats"] for p in pts]
    order = sorted(range(len(pts)), key=lambda k: pts[k]["H_tail_nats"], reverse=True)[:mark_top_n]

    fig, axes = plt.subplots(2, 1, figsize=(14, 7))
    ax = axes[0]
    ax.plot(xs, H, lw=0.8, color="#8e44ad", marker="." if not scope_all else None, ms=3)
    ax.scatter([xs[k] for k in order], [H[k] for k in order],
               color="crimson", s=28, zorder=3, label=f"top-{mark_top_n} highest entropy")
    for k in order:
        ax.annotate(pts[k]["token"].replace("\n", "⏎"), (xs[k], H[k]),
                    fontsize=7, xytext=(2, 4), textcoords="offset points", color="crimson")
    s = gold["summary"]
    ax.set_ylabel("entropy H_tail (nats)")
    ax.set_title(f"per-CODE-token predictive entropy  (gold n={s['n_tokens_gold']} of "
                 f"{s['n_tokens_total']} tokens; comment+ws excluded)  "
                 f"mean={s['H_tail_nats_mean']:.3f} nats = {s['H_tail_bits_mean']:.3f} bits")
    ax.set_xlabel("token position (original index; gaps = removed comment/ws)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[1]
    ax.hist(H, bins=60, color="#8e44ad", edgecolor="white")
    ax.set_yscale("log")
    ax.set_xlabel("per-token entropy (nats)")
    ax.set_ylabel("count (log)")
    ax.set_title("entropy distribution over code tokens  (mass near 0 = confident; right tail = decision hotspots)")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    print(f"\n[+] plot saved to {out_png}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("json_path", nargs="?", default=str(HERE / "dsl_out_x3872.json"))
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--out", default=None, help="write gold JSON here (default: gold_<stem>.json)")
    ap.add_argument("--plot", default=None, help="write PNG here (default: gold_<stem>.png)")
    ap.add_argument("--keep-all", action="store_true",
                    help="compute gold over ALL tokens (do not drop comment/whitespace)")
    args = ap.parse_args()

    src = Path(args.json_path)
    toks, meta = load_tokens(src)
    gold = build_gold(toks, meta, src, keep_all=args.keep_all)

    s = gold["summary"]
    print(f"# loaded {s['n_tokens_total']} content tokens from {src}")
    print(f"  kinds: code={s['n_tokens_code']}  comment={s['n_tokens_comment']}  ws={s['n_tokens_ws']}")
    print(f"  gold scope: {gold['gold_scope']}  ->  {s['n_tokens_gold']} tokens scored")
    print(f"  mean entropy  = {s['H_tail_nats_mean']:.4f} nats  = {s['H_tail_bits_mean']:.4f} bits")
    print(f"  sum  entropy  = {s['H_tail_nats_sum']:.4f} nats")
    print(f"  min/max       = {s['H_tail_nats_min']:.4f} / {s['H_tail_nats_max']:.4f} nats")
    print(f"  mean tail mass= {s['mean_tail_mass']:.4g}  (top-20 truncation leakage)")

    top_uncertain(gold, args.top)

    out = Path(args.out) if args.out else HERE / f"gold_{src.stem}.json"
    out.write_text(json.dumps(gold, ensure_ascii=False, indent=2))
    print(f"\n[+] gold written to {out}  (per_token[].H_tail_nats over code tokens = gold distribution)")

    png = Path(args.plot) if args.plot else HERE / f"gold_{src.stem}.png"
    plot(gold, png)


if __name__ == "__main__":
    main()
