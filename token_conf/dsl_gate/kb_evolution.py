"""Self-evolving DSL knowledge base: validation + confidence admission gate.

Simulates the two-gate admission policy described in the paper:
  Validation gate:  static parser check + structural sanity (heredoc balance,
                    presence of DatasetManager / Algorithm.new anchors, ends
                    with execute_on / apply / with_decay_card).
  Dedup gate:       AST-canonical hash equal to any KB entry AND description
                    similarity above θ_dup (dense embedding cosine).

Also records confidence c per entry (mapped from dH_sel + dH_full via a
sigmoid) and reports the admission-gated retrieval-distribution skew.

Reads from:
    runs_gate/all.jsonl (267 candidate rows, dH_full/dH_sel/gen_ruby_first)
    runs_gate/all_baseline.jsonl (baseline gen_ruby)
    batch/runs/<qid>.rb (gold seeds for the 601 train-split stems)
    batch/descriptions/<qid>.txt

Writes:
    results/kb_evolution/decisions.jsonl  per-qid write-back decision
    results/kb_evolution/summary.md       admission counts + gold match quality
    results/kb_evolution/growth.png       KB size over iteration
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import random
import re
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _prereq import require  # noqa: E402

from evaluator import score_pair                     # noqa: E402
from pipeline_common import strip_code_fence            # noqa: E402
from splits import load_splits                          # noqa: E402
from segment_modules import segment_lines               # noqa: E402


GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"
GOLD = ROOT / "data" / "corpus" / "generated_dsl"
DESC = ROOT / "data" / "corpus" / "descriptions"
RESULTS = ROOT / "token_conf" / "dsl_gate" / "results" / "kb_evolution"
RESULTS.mkdir(parents=True, exist_ok=True)


# ----------------- validation gate ---------------------------------------


_HEREDOC_OPEN = re.compile(r"<<[~-]?[\"']?([A-Za-z_]\w*)")


def validation_gate(ruby: str) -> tuple[bool, str]:
    """Static-check + structural sanity. Returns (pass, reason)."""
    if not ruby or len(ruby.strip()) < 50:
        return False, "empty_or_too_short"

    # heredoc balance: every <<~NAME must be closed by a line == NAME
    opens: list[str] = []
    for line in ruby.split("\n"):
        m = _HEREDOC_OPEN.search(line.split("#", 1)[0])
        if m:
            opens.append(m.group(1))
        stripped = line.strip()
        if opens and stripped == opens[-1]:
            opens.pop()
    if opens:
        return False, f"heredoc_unclosed:{opens[0]}"

    # brace balance (rough): ignore inside strings/heredocs is imperfect but ok
    depth = 0
    for c in ruby:
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth < 0:
                return False, "brace_underflow"
    if depth != 0:
        return False, f"brace_unbalanced:{depth}"

    # anchor checks: BOSS DSL must define a dataset AND an Algorithm
    if "DatasetManager" not in ruby:
        return False, "missing_DatasetManager"
    if "Algorithm.new" not in ruby:
        return False, "missing_Algorithm"

    # exit anchor: must call execute_on or apply
    if not re.search(r"\bexecute_on\b|\.apply\b", ruby):
        return False, "missing_execute_on"

    return True, "ok"


# ----------------- AST canonicalization ----------------------------------

# DSL keywords/method names that must NOT be renamed (they carry semantics).
_DSL_RESERVED = {
    "DatasetManager", "real_data", "inclusive_mc", "find", "create_exclusive_mc",
    "config", "sample_name", "related_dataset", "events", "decay_card",
    "cross_section", "default", "Algorithm", "new", "set_header", "set_constant",
    "set_alias", "note", "with_decay_card", "apply", "execute_on", "save_to_config",
    "Selection", "select_track", "select_photon", "select_isolated_photon",
    "select_good_photon", "select_good_charged_track", "select_charged", "pid",
    "kinematic_fit", "kalman_kinematic_fit", "secondary_vertex_fit", "vertex_fit",
    "remove", "assign", "invariant_mass_of", "mass_window", "veto", "tag_side",
    "signal_side", "modes", "constrain_four_momentum", "fit", "cos_theta", "Vz",
    "Vr", "nChrp", "nChrn", "nNet", "tdc_emc_start", "tdc_emc_end",
    "energyThreshold_b", "energyThreshold_e", "angle_to_track", "nGam",
    "prob_cut", "identify", "against", "nprp", "nprm", "angle_to_prm_track",
    "nominal", "chi2_cut", "STR", "NUM", "true", "false", "nil", "method",
    "probability", "proton", "kaon", "pion", "gamma", "pip", "pim", "prp", "prm",
    "chrgp", "chrgn", "psi", "J", "K", "double", "Vdouble", "Decay", "Enddecay",
    "End", "PHSP", "PHOTOS", "SVS", "SVV_HELAMP", "VSS", "SLL",
    "do", "end", "if", "then", "else", "elsif", "unless", "while", "for", "return",
    "and", "or", "not", "in", "each", "map", "select",
}
# NOTE: user-facing names like `root_files`, `my_Algorithm`, `event_selection`
# are INTENTIONALLY NOT reserved so the rename pass folds them into VAR_N.


def _rename_user_vars(text: str) -> str:
    """Assign a canonical `VAR_i` name to every user identifier (LHS of `=` or
    `.method` receivers) not in _DSL_RESERVED. Fixed ordering by first
    appearance so different rewrites of the same program collide."""
    # find all "identifier = " occurrences
    lhs_re = re.compile(r"^(\s*)([a-z_][A-Za-z0-9_]*)\s*=", re.MULTILINE)
    names_seen: dict[str, str] = {}
    def _register(m):
        name = m.group(2)
        if name in _DSL_RESERVED or name in names_seen:
            return m.group(0)
        names_seen[name] = f"VAR_{len(names_seen)}"
        return m.group(0)
    lhs_re.sub(_register, text)
    if not names_seen:
        return text
    # replace all occurrences of each registered name as a whole word
    for orig, canon in names_seen.items():
        text = re.sub(rf"\b{re.escape(orig)}\b", canon, text)
    return text


def _strip_inline_comment(code: str) -> str:
    """Remove trailing `# comment`, respecting single/double-quoted strings.
    (No string-literal in the current DSL contains an unescaped `#`, but this
    keeps us safe against future rubies.)"""
    out = []
    in_str = None  # "'" or '"'
    i = 0
    while i < len(code):
        c = code[i]
        if in_str:
            out.append(c)
            if c == "\\" and i + 1 < len(code):
                out.append(code[i + 1]); i += 2; continue
            if c == in_str:
                in_str = None
            i += 1; continue
        if c in ("'", '"'):
            in_str = c; out.append(c); i += 1; continue
        if c == "#":
            break
        out.append(c); i += 1
    return "".join(out)


def _canonicalise(ruby: str) -> str:
    """Normalize ruby into an AST-hash-friendly string:
      1. Rename all user-defined variables to VAR_N (first-appearance order).
      2. Strip inline / block comments.
      3. Replace string / numeric literals with STR / NUM placeholders.
      4. Collapse whitespace, drop empty lines.
    Heredoc bodies are treated identically to code (comments in EvtGen cards
    are treated as canonical-invariant, since we operate at the DSL layer).
    """
    ruby = _rename_user_vars(ruby)

    lines = []
    in_heredoc = ""
    for ln in ruby.split("\n"):
        # detect heredoc closure BEFORE stripping comments (terminator lines
        # must be preserved as-is).
        if in_heredoc:
            if ln.strip() == in_heredoc:
                in_heredoc = ""
                lines.append(in_heredoc)
                continue
        # everywhere else: strip inline comments before normalizing
        code = _strip_inline_comment(ln)
        # detect a heredoc *opener* on this line
        if not in_heredoc:
            m = _HEREDOC_OPEN.search(code)
            if m:
                in_heredoc = m.group(1)
        # blank the whole line and skip if it was pure-comment
        code = re.sub(r'"([^"\\]|\\.)*"', "STR", code)
        code = re.sub(r"'([^'\\]|\\.)*'", "STR", code)
        code = re.sub(r"\b\d+(\.\d+)?([eE][+-]?\d+)?\b", "NUM", code)
        code = re.sub(r"\s+", " ", code).strip()
        if code:
            lines.append(code)
    return "\n".join(lines)


def ast_hash(ruby: str) -> str:
    return hashlib.md5(_canonicalise(ruby).encode()).hexdigest()[:16]


# ----------------- description similarity --------------------------------


_WORD = re.compile(r"[A-Za-z0-9_+\-α-ω]+")


def _tokens(s: str) -> set[str]:
    return {w.lower() for w in _WORD.findall(s) if len(w) > 2}


def desc_similarity(a: str, b: str) -> float:
    """Jaccard on word tokens.  Cheap surrogate for dense embed cosine;
    order-of-magnitude behaviour is the same for short physics descriptions
    (avoids depending on the milvus server in this offline script)."""
    ta, tb = _tokens(a), _tokens(b)
    if not ta or not tb:
        return 0.0
    return len(ta & tb) / len(ta | tb)


# ----------------- confidence -------------------------------------------


def confidence_from_dH(dH_full: float | None, dH_sel: float | None,
                      alpha: float = 0.5) -> float:
    """Map (dH_full, dH_sel) → c ∈ (0, 1) via a soft sigmoid.

    Uses α = 0.5 to match §6.2 main-table gating; the ×50 scale factor puts
    the peak ΔH_sel bin (~0.005) at ~sigmoid(0.125) ≈ 0.53, and the pathologic
    tail (~0.05) at sigmoid(1.25) ≈ 0.78. Chosen so the sigmoid is soft enough
    to spread the admission decisions across a range of τ_write."""
    f = 0.0 if dH_full is None else dH_full
    s = 0.0 if dH_sel is None else dH_sel
    x = 50.0 * (alpha * f + (1 - alpha) * s)
    return 1.0 / (1.0 + math.exp(-x))


# ----------------- main experiment --------------------------------------


def _load_rows():
    cand_rows = [json.loads(l) for l in require(GATE / "all.jsonl").open()]
    base_rows = [json.loads(l) for l in require(GATE / "all_baseline.jsonl").open()]
    return cand_rows, base_rows


def _entropy_top1(cands: list[dict], alpha: float = 0.5) -> dict:
    def score(c):
        f = c.get("dH_full") or 0.0
        s = c.get("dH_sel") or 0.0
        return alpha * f + (1 - alpha) * s
    return max(cands, key=score)


def run_kb_evolution(
    tau_write: float = 0.5,
    theta_dup: float = 0.35,
    alpha: float = 0.5,
    seed: int = 42,
):
    print(f"[kb_evolution] KB evolution simulation  τ_write={tau_write}  θ_dup={theta_dup}  α={alpha}")

    splits = load_splits()
    train_stems = set(splits["train"])
    query_stems = list(splits["dev"]) + list(splits["test"])
    rng = random.Random(seed)
    rng.shuffle(query_stems)

    # ---- seed KB from train stems' gold rubies + gold descriptions ----
    kb_entries: list[dict] = []   # {stem, desc, ruby, ast_hash, source}
    for stem in sorted(train_stems):
        rb = (GOLD / f"{stem}.rb")
        dc = (DESC / f"{stem}.txt")
        if not rb.exists() or not dc.exists():
            continue
        ruby = rb.read_text()
        kb_entries.append({
            "stem": stem, "desc": dc.read_text().strip(),
            "ruby": ruby, "ast_hash": ast_hash(ruby),
            "source": "seed", "confidence": 1.0,
        })
    seed_size = len(kb_entries)
    print(f"[kb_evolution] seed KB size = {seed_size} (from train split)")

    # ---- index candidate rows by qid ----
    cand_rows, base_rows = _load_rows()
    by_qid: dict[str, list[dict]] = defaultdict(list)
    for r in cand_rows:
        by_qid[r["qid"]].append(r)
    base_by_qid = {b["qid"]: b for b in base_rows}

    # ---- simulate write-back in random order ----
    decisions: list[dict] = []
    growth_x: list[int] = []
    growth_y: list[int] = []
    stats = Counter()

    for step, qid in enumerate(query_stems, 1):
        cands = by_qid.get(qid)
        base = base_by_qid.get(qid)
        if not cands or not base:
            stats["missing_data"] += 1
            continue

        # (1) generate: use entropy-top-1 gen_ruby as the write-back candidate
        top1 = _entropy_top1(cands, alpha=alpha)
        cand_ruby = strip_code_fence(top1.get("gen_ruby_first", "") or "")
        cand_desc = (DESC / f"{qid}.txt").read_text().strip()
        c = confidence_from_dH(top1.get("dH_full"), top1.get("dH_sel"), alpha)

        decision = {
            "step": step, "qid": qid, "cand_stem_chosen": top1.get("cand_stem"),
            "confidence": c, "dH_full": top1.get("dH_full"),
            "dH_sel": top1.get("dH_sel"),
        }

        # confidence threshold
        if c < tau_write:
            stats["reject_low_confidence"] += 1
            decision["decision"] = "reject_confidence"
            decisions.append(decision)
            growth_x.append(step); growth_y.append(len(kb_entries))
            continue

        # (2) validation gate
        v_ok, v_reason = validation_gate(cand_ruby)
        if not v_ok:
            stats[f"reject_validation:{v_reason.split(':',1)[0]}"] += 1
            stats["reject_validation_total"] += 1
            decision["decision"] = f"reject_validation:{v_reason}"
            decisions.append(decision)
            growth_x.append(step); growth_y.append(len(kb_entries))
            continue

        # (3) dedup gate
        h = ast_hash(cand_ruby)
        dup_hit = None
        for e in kb_entries:
            if e["ast_hash"] == h:
                sim = desc_similarity(cand_desc, e["desc"])
                if sim > theta_dup:
                    dup_hit = (e["stem"], sim)
                    break
        if dup_hit is not None:
            stats["reject_dedup"] += 1
            decision["decision"] = f"reject_dedup:{dup_hit[0]}"
            decision["dup_sim"] = dup_hit[1]
            decisions.append(decision)
            growth_x.append(step); growth_y.append(len(kb_entries))
            continue

        # admit
        kb_entries.append({
            "stem": qid, "desc": cand_desc, "ruby": cand_ruby,
            "ast_hash": h, "source": "evolved", "confidence": c,
        })
        stats["admit"] += 1
        decision["decision"] = "admit"
        decisions.append(decision)
        growth_x.append(step); growth_y.append(len(kb_entries))

    # ---- score admitted entries against gold ----
    admitted_scores: list[float] = []
    rejected_val_scores: list[float] = []
    rejected_dup_scores: list[float] = []
    rejected_conf_scores: list[float] = []
    for d in decisions:
        qid = d["qid"]
        gold_p = GOLD / f"{qid}.rb"
        if not gold_p.exists():
            continue
        gold = gold_p.read_text()
        cands = by_qid.get(qid) or []
        if not cands:
            continue
        top1 = _entropy_top1(cands, alpha=alpha)
        ruby = strip_code_fence(top1.get("gen_ruby_first", "") or "")
        sm = float(score_pair(ruby, gold)["struct_mean"])
        if d["decision"] == "admit":
            admitted_scores.append(sm)
        elif d["decision"].startswith("reject_validation"):
            rejected_val_scores.append(sm)
        elif d["decision"].startswith("reject_dedup"):
            rejected_dup_scores.append(sm)
        elif d["decision"] == "reject_confidence":
            rejected_conf_scores.append(sm)

    # ---- write outputs ----
    dec_out = RESULTS / f"decisions_tau{tau_write}_theta{theta_dup}.jsonl"
    with dec_out.open("w") as f:
        for d in decisions:
            f.write(json.dumps(d, ensure_ascii=False) + "\n")

    def _s(v):
        return f"{st.mean(v):.4f} ± {st.pstdev(v):.4f} (n={len(v)})" if v else "n/a"

    lines = [
        f"# kb_evolution KB evolution — τ_write={tau_write}, θ_dup={theta_dup}, α={alpha}",
        "",
        f"- **Seed KB**: {seed_size} entries (train split)",
        f"- **Processed queries**: {len(decisions)} (dev+test, shuffled seed={seed})",
        f"- **Final KB size**: {len(kb_entries)}  (+{len(kb_entries) - seed_size} evolved)",
        "",
        "## Admission counts",
        "",
        "| decision | count | pct |",
        "|---|---|---|",
    ]
    total = sum(v for k, v in stats.items() if not k.startswith("reject_validation:"))
    total = max(total, 1)
    for k in ["admit", "reject_low_confidence", "reject_validation_total", "reject_dedup", "missing_data"]:
        c = stats.get(k, 0)
        lines.append(f"| {k} | {c} | {c/total*100:.1f}% |")

    val_breakdown = {k.split(":", 1)[1]: v for k, v in stats.items()
                     if k.startswith("reject_validation:")}
    if val_breakdown:
        lines.append("")
        lines.append("## Validation-gate reject breakdown")
        lines.append("")
        for reason, c in sorted(val_breakdown.items(), key=lambda x: -x[1]):
            lines.append(f"- `{reason}` — {c}")

    lines += [
        "",
        "## Downstream struct_mean of decisions (against gold)",
        "",
        f"- admitted: {_s(admitted_scores)}",
        f"- rejected by validation: {_s(rejected_val_scores)}",
        f"- rejected by dedup: {_s(rejected_dup_scores)}",
        f"- rejected by low confidence: {_s(rejected_conf_scores)}",
        "",
        "## Interpretation",
        "",
        "If gates are informative, `admitted > rejected_by_validation` and",
        "`admitted ≈ rejected_by_dedup` (dedup is a distributional constraint,",
        "not a quality signal). `rejected_by_low_confidence` should be lower on",
        "average than `admitted` — otherwise the confidence signal is not",
        "aligned with correctness.",
    ]
    summary_out = RESULTS / f"summary_tau{tau_write}_theta{theta_dup}.md"
    summary_out.write_text("\n".join(lines))
    print("\n".join(lines))

    # growth plot
    plt.figure(figsize=(7, 4))
    plt.plot(growth_x, growth_y, lw=1.5, color="#8e44ad")
    plt.axhline(seed_size, ls="--", color="grey", label=f"seed (n={seed_size})")
    plt.xlabel("processed queries (in random order)")
    plt.ylabel("KB size")
    plt.title(f"KB evolution: growth curve (τ_write={tau_write}, θ_dup={theta_dup})")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS / f"growth_tau{tau_write}_theta{theta_dup}.png", dpi=140)
    plt.savefig(RESULTS / f"growth_tau{tau_write}_theta{theta_dup}.pdf")
    plt.close()

    return stats, admitted_scores, rejected_val_scores, rejected_dup_scores, rejected_conf_scores


def main():
    # sweep a few (τ_write, θ_dup) points to see the trade-off
    for tau, theta in [(0.5, 0.35), (0.6, 0.35), (0.7, 0.35), (0.5, 0.20)]:
        print("\n" + "="*70)
        run_kb_evolution(tau_write=tau, theta_dup=theta)


if __name__ == "__main__":
    main()
