"""Evaluator: compare a generated Ruby DSL to its reference counterpart.

Three tiers, from strict to loose:

  strict (numbers matter):
    * dataset_em     — set of DatasetManager.*.find("<id>") argument strings
                       (order-independent EM on the ids).
    * decay_card_em  — normalized EvtGen decay card text EM (whitespace collapsed,
                       trailing colons/semicolons kept as separators).
                       Averaged over cards when both sides have >=1.

  structural (methods & knob NAMES matter, values do NOT):
    * sel_method_f1  — set F1 over the top-level chain methods on
                       `Selection.new` (select_track/select_photon/pid/
                       select_isolated_photon/kinematic_fit/remove/assign/…).
    * sel_slot_f1    — for each shared method, the F1 over the SET of keys
                       inside its `{ … }` body (cos_theta, Vz, nGam, prob_cut,
                       chi2_cut, angle_to_track, energyThreshold_b, …).
                       Slot F1 = mean over intersected methods.
    * struct_mean    — 0.5 * sel_method_f1 + 0.5 * sel_slot_f1  (the main metric).

  soft (placeholder):
    * llm_judge      — reserved hook for a small LLM-as-judge pass; not implemented.

CLI:
    python token_conf/dsl_gate/evaluator.py \\
        --pairs '{"pred_stem":"1506.06018v2","gold_stem":"1506.06018v2"}' \\
        --gold-dir token_conf/batch/runs \\
        --pred-dir token_conf/batch/runs

Or import score_pair(pred_text, gold_text) -> dict.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


# ------------- helpers ---------------------------------------------------


_FIND_ARG = re.compile(r'DatasetManager\.(?:real_data|inclusive_mc)\.find\(\s*"([^"]+)"\s*\)')


def _dataset_ids(text: str) -> set[str]:
    return set(_FIND_ARG.findall(text))


_HEREDOC_OPEN = re.compile(r"<<[~-]?[\"']?([A-Za-z_]\w*)")


def _decay_cards(text: str) -> list[str]:
    """Extract each heredoc body whose terminator name contains 'DECAY' (case-insensitive).
    Returns the raw inner text of every card, in order."""
    lines = text.split("\n")
    cards: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = _HEREDOC_OPEN.search(line)
        if m and "DECAY" in m.group(1).upper():
            term = m.group(1)
            body: list[str] = []
            j = i + 1
            while j < len(lines):
                if lines[j].strip() == term:
                    break
                body.append(lines[j])
                j += 1
            cards.append("\n".join(body))
            i = j + 1
            continue
        i += 1
    return cards


def _normalize_decay(card: str) -> str:
    """Collapse whitespace and strip line-terminating comments so equivalent
    EvtGen cards compare equal.  Preserves semicolons/colons as tokens."""
    out_lines: list[str] = []
    for ln in card.split("\n"):
        # strip # comments (uncommon in cards but safe)
        code = ln.split("#", 1)[0].strip()
        if not code:
            continue
        # collapse internal whitespace
        code = re.sub(r"\s+", " ", code)
        out_lines.append(code)
    return "\n".join(out_lines).strip()


# ------------- selection segmentation -----------------------------------
#
# We reuse hepuke's segmenter idea but keep this file self-contained. We
# tokenize the `Selection.new … chain by scanning the file at brace depth 0
# for the top-level chain methods; then for each such method we extract the
# `{ … }` body (matched by depth) and read its knob KEY set.


_CHAIN_METHOD_RE = re.compile(
    r"\.(select_track|select_photon|select_good_photon|select_isolated_photon|"
    r"pid|kinematic_fit|kalman_kinematic_fit|secondary_vertex_fit|vertex_fit|"
    r"remove|assign|invariant_mass_of|mass_window|veto|tag_side|signal_side|"
    r"modes|constrain_four_momentum|fit|select_charged|select_good_charged_track)"
    r"(?:\s*\(([^)]*)\))?"                  # optional (…) args
    r"(?:\s*\{)?",                          # optional opening brace
)

_KNOB_LINE_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\b")


def _selection_region(text: str) -> str:
    """Return the substring of `text` between the first `Selection.new` occurrence
    and the end of its chain (approximated by end-of-file — the chain runs to
    the next Selection or execute_on, whichever comes first)."""
    m = re.search(r"\bSelection\.new\b", text)
    if not m:
        return ""
    start = m.end()
    # end at next `execute_on`, `my_Algorithm.with_decay_card`, or another Selection.new
    stops = []
    for pat in (r"\bexecute_on\b", r"\.with_decay_card\b", r"\bSelection\.new\b"):
        m2 = re.search(pat, text[start + 1:])
        if m2:
            stops.append(start + 1 + m2.start())
    end = min(stops) if stops else len(text)
    return text[start:end]


def _find_matching_brace(text: str, open_idx: int) -> int:
    """Return index of `}` that matches the `{` at open_idx (assumed). Naive:
    tracks brace depth ignoring strings and comments.  Good enough for our DSL."""
    depth = 0
    i = open_idx
    n = len(text)
    while i < n:
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def _chain_methods_with_bodies(sel_text: str) -> list[tuple[str, str]]:
    """Return [(method_name, body_text)] for each top-level `.method { … }`
    encountered in the Selection chain. Methods without a body get body=''."""
    out: list[tuple[str, str]] = []
    depth = 0
    i = 0
    n = len(sel_text)
    while i < n:
        c = sel_text[i]
        if c == "{":
            depth += 1
            i += 1
            continue
        if c == "}":
            depth -= 1
            i += 1
            continue
        if depth == 0 and c == ".":
            m = _CHAIN_METHOD_RE.match(sel_text, i)
            if m:
                method = m.group(1)
                body = ""
                after = m.end()
                # find optional brace body immediately after
                # (either matched by the regex '{' or a separate '{' with ws)
                brace_i = after - 1 if sel_text[after - 1:after] == "{" else -1
                if brace_i == -1:
                    # search for `{` in the immediate whitespace
                    k = after
                    while k < n and sel_text[k] in " \t":
                        k += 1
                    if k < n and sel_text[k] == "{":
                        brace_i = k
                if brace_i != -1:
                    end_brace = _find_matching_brace(sel_text, brace_i)
                    if end_brace != -1:
                        body = sel_text[brace_i + 1:end_brace]
                        out.append((method, body))
                        i = end_brace + 1
                        continue
                out.append((method, body))
                i = m.end()
                continue
        i += 1
    return out


def _knobs_of(body: str) -> set[str]:
    """Extract KEY names from a chain-method body. One key per non-blank line's
    first identifier."""
    keys: set[str] = set()
    for ln in body.split("\n"):
        s = ln.split("#", 1)[0].strip()
        if not s:
            continue
        m = _KNOB_LINE_RE.match(s)
        if m:
            keys.add(m.group(1))
    return keys


def _knob_pairs(body: str) -> set[tuple[str, str]]:
    """Extract (key, normalized_value) pairs from a chain-method body.

    Strict tier: values matter. One pair per non-blank line's first identifier
    plus its trailing value (whitespace-collapsed). Numeric forms are NOT
    normalised (0.9 != 0.90), so this is deliberately harsher than slot-F1."""
    pairs: set[tuple[str, str]] = set()
    for ln in body.split("\n"):
        s = ln.split("#", 1)[0].strip()
        if not s:
            continue
        m = _KNOB_LINE_RE.match(s)
        if not m:
            continue
        key = m.group(1)
        val = re.sub(r"\s+", " ", s[m.end():].strip())
        pairs.add((key, val))
    return pairs


def _chain_method_sigs(sel_text: str) -> set[tuple[str, str]]:
    """Top-level chain methods with normalised args: {(method, args)}.

    Mirrors _chain_methods_with_bodies but keeps the optional (…) arguments so
    the strict method tier distinguishes `invariant_mass_of(:pipi)` from
    `invariant_mass_of(:pim)."""
    out: set[tuple[str, str]] = set()
    depth = 0
    i = 0
    n = len(sel_text)
    while i < n:
        c = sel_text[i]
        if c == "{":
            depth += 1
            i += 1
            continue
        if c == "}":
            depth -= 1
            i += 1
            continue
        if depth == 0 and c == ".":
            m = _CHAIN_METHOD_RE.match(sel_text, i)
            if m:
                args = re.sub(r"\s+", " ", (m.group(2) or "").strip())
                out.add((m.group(1), args))
                after = m.end()
                # _CHAIN_METHOD_RE consumes an optional `{`; if it did, skip to
                # the matching `}` exactly like _chain_methods_with_bodies so we
                # do not rescan the brace body at depth 0.
                if after >= 1 and sel_text[after - 1:after] == "{":
                    end_brace = _find_matching_brace(sel_text, after - 1)
                    i = (end_brace + 1) if end_brace != -1 else after
                else:
                    i = after
                continue
        i += 1
    return out


# ------------- F1 helpers -----------------------------------------------


def _f1(pred: set, gold: set) -> float:
    if not pred and not gold:
        return 1.0
    if not pred or not gold:
        return 0.0
    tp = len(pred & gold)
    if tp == 0:
        return 0.0
    p = tp / len(pred)
    r = tp / len(gold)
    return 2 * p * r / (p + r)


# ------------- top-level scoring -----------------------------------------


def score_pair(pred_text: str, gold_text: str) -> dict:
    """Compute all metrics for one (pred, gold) DSL pair."""
    # ---- strict ----
    pred_ds = _dataset_ids(pred_text)
    gold_ds = _dataset_ids(gold_text)
    if not gold_ds:
        dataset_em = 1.0 if not pred_ds else 0.0
    else:
        dataset_em = 1.0 if pred_ds == gold_ds else 0.0

    pred_cards = [_normalize_decay(c) for c in _decay_cards(pred_text)]
    gold_cards = [_normalize_decay(c) for c in _decay_cards(gold_text)]
    if not gold_cards:
        decay_card_em = 1.0 if not pred_cards else 0.0
    else:
        # bipartite-ish: for each gold, best-matching pred is a set match
        pred_set, gold_set = set(pred_cards), set(gold_cards)
        matched = len(pred_set & gold_set)
        decay_card_em = matched / max(len(gold_cards), len(pred_cards))

    # ---- structural ----
    pred_chain = _chain_methods_with_bodies(_selection_region(pred_text))
    gold_chain = _chain_methods_with_bodies(_selection_region(gold_text))
    pred_methods = {m for m, _ in pred_chain}
    gold_methods = {m for m, _ in gold_chain}
    sel_method_f1 = _f1(pred_methods, gold_methods)

    # slot F1: for methods present in BOTH, compare the key sets of the FIRST
    # occurrence of each method (chain usually has one occurrence per method;
    # if multiple, we union across occurrences per side).
    def _union_bodies(chain: list[tuple[str, str]]) -> dict[str, set[str]]:
        out: dict[str, set[str]] = {}
        for m, b in chain:
            out.setdefault(m, set()).update(_knobs_of(b))
        return out

    pred_slots = _union_bodies(pred_chain)
    gold_slots = _union_bodies(gold_chain)
    common = pred_methods & gold_methods
    slot_f1s = [_f1(pred_slots.get(m, set()), gold_slots.get(m, set())) for m in common]
    sel_slot_f1 = (sum(slot_f1s) / len(slot_f1s)) if slot_f1s else 0.0

    struct_mean = 0.5 * sel_method_f1 + 0.5 * sel_slot_f1

    # ---- strict structural (names AND values) ----
    strict_method_f1 = _f1(_chain_method_sigs(_selection_region(pred_text)),
                           _chain_method_sigs(_selection_region(gold_text)))

    def _union_pairs(chain: list[tuple[str, str]]) -> dict[str, set[tuple[str, str]]]:
        out: dict[str, set[tuple[str, str]]] = {}
        for m, b in chain:
            out.setdefault(m, set()).update(_knob_pairs(b))
        return out

    pred_pairs = _union_pairs(pred_chain)
    gold_pairs = _union_pairs(gold_chain)
    strict_slot_f1s = [_f1(pred_pairs.get(m, set()), gold_pairs.get(m, set()))
                       for m in common]
    strict_slot_f1 = (sum(strict_slot_f1s) / len(strict_slot_f1s)) if strict_slot_f1s else 0.0

    strict_struct_mean = 0.5 * strict_method_f1 + 0.5 * strict_slot_f1

    return {
        "dataset_em": dataset_em,
        "decay_card_em": decay_card_em,
        "sel_method_f1": sel_method_f1,
        "sel_slot_f1": sel_slot_f1,
        "struct_mean": struct_mean,
        "strict_method_f1": strict_method_f1,
        "strict_slot_f1": strict_slot_f1,
        "strict_struct_mean": strict_struct_mean,
        "n_pred_ds": len(pred_ds),
        "n_gold_ds": len(gold_ds),
        "n_pred_cards": len(pred_cards),
        "n_gold_cards": len(gold_cards),
        "n_pred_methods": len(pred_methods),
        "n_gold_methods": len(gold_methods),
    }


# ------------- CLI -------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", type=Path, required=False,
                    help="one predicted .rb (single-pair mode)")
    ap.add_argument("--gold", type=Path, required=False,
                    help="one gold .rb (single-pair mode)")
    ap.add_argument("--sanity", action="store_true",
                    help="score a gold .rb against itself as sanity check")
    ap.add_argument("--stem", default="1001.5328v1")
    args = ap.parse_args()

    if args.sanity:
        gold = ROOT / "data" / "corpus" / "generated_dsl" / f"{args.stem}.rb"
        text = gold.read_text(encoding="utf-8")
        print(json.dumps(score_pair(text, text), indent=2))
        return 0

    if not (args.pred and args.gold):
        ap.error("supply --pred and --gold, or --sanity")
    p = args.pred.read_text(encoding="utf-8")
    g = args.gold.read_text(encoding="utf-8")
    print(json.dumps(score_pair(p, g), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
