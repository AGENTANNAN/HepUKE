"""Segment a generated BOSS DSL into TEMPLATE MODULES and tag every token.

The corpus follows one skeleton (verified over 476 generated DSLs). We split each
file into 8 ordered modules and assign every content token to exactly one:

    1 preamble      loader/section-header comments before the first dataset line
    2 dataset       DatasetManager.real_data / .inclusive_mc .find(...)
    3 decay_card    heredoc EvtGen decay card(s)  (<<~DECAYCARD ... Decay ...)
    4 exclusive_mc  create_exclusive_mc do |config| ... end  (+ config.* lines)
    5 algorithm     Algorithm.new + set_header/set_constant/set_alias
    6 selection     Selection.new -> select_track/select_photon/pid/kinematic_fit chain
    7 note          .note(...) blocks (constraints not expressible in the DSL)
    8 execute       with_decay_card/apply/execute_on/root_files/save_to_config

Method: a LINE-level state machine over the reconstructed source. Heredoc bodies are
tracked explicitly (reusing the same rule as the code/comment lexer). A line that
matches a module anchor switches the current module; a continuation line with no
anchor inherits the current module (so `config.*`, cut lines, chained `.set_*`
stay with their opener). Each token inherits the module of its first CODE char
(comment/ws tokens inherit the surrounding module but are excluded from gold).
"""
from __future__ import annotations

import re

MODULES = ["preamble", "dataset", "decay_card", "exclusive_mc",
           "algorithm", "selection", "note", "execute"]

_HEREDOC_OPEN = re.compile(r"<<[~-]?[\"']?([A-Za-z_]\w*)")

# line anchors, checked in order; first hit sets the module for that line
_ANCHORS = [
    ("execute",      re.compile(r"\b(execute_on|with_decay_card|\.apply\b|save_to_config|root_files\s*=)")),
    ("note",         re.compile(r"\.note\b")),
    ("exclusive_mc", re.compile(r"\bcreate_exclusive_mc(_for)?\b")),
    ("dataset",      re.compile(r"DatasetManager\.(real_data|inclusive_mc)")),
    ("algorithm",    re.compile(r"\bAlgorithm\.new\b|\bset_header\b|\bset_constant\b|\bset_alias\b")),
    ("selection",    re.compile(r"\bSelection\.new\b|\bevent_selection\b|\bselect_\w+|\bpid\b|"
                                r"\bkinematic_fit\b|\.remove\b|\.assign\b|\bveto\b|\bmass_window\b")),
]


def segment_lines(text: str) -> list[str]:
    """Return a per-LINE module label list (len == number of lines in text)."""
    lines = text.split("\n")
    labels: list[str] = []
    cur = "preamble"
    seen_code = False
    heredoc_term = ""      # active heredoc terminator, or "" if not in a heredoc
    for ln in lines:
        stripped = ln.strip()

        # --- inside a heredoc body: it belongs to whatever opened it (decay_card) ---
        if heredoc_term:
            labels.append(cur)
            if stripped == heredoc_term:
                heredoc_term = ""      # terminator line closes the heredoc
            continue

        # blank / comment-only line: inherit current module (preamble if nothing yet)
        code_part = stripped.split("#", 1)[0].strip()
        if code_part:
            seen_code = True
            matched = None
            for mod, rx in _ANCHORS:
                if rx.search(code_part):
                    matched = mod
                    break
            # a bare heredoc opener line for a decay card
            m = _HEREDOC_OPEN.search(code_part)
            if m and ("DECAY" in m.group(1).upper() or "decay" in code_part.lower()):
                matched = "decay_card"
            if matched:
                cur = matched
        elif not seen_code:
            cur = "preamble"

        labels.append(cur)

        # open a heredoc if this line starts one (body begins next line)
        m = _HEREDOC_OPEN.search(code_part) if code_part else None
        if m:
            heredoc_term = m.group(1)
    return labels


def tag_modules(toks: list[dict]) -> None:
    """Annotate each token in-place with 'module' (one of MODULES).

    Assumes tokens already carry char offsets implicitly via concatenation order.
    A token's module = the module of the LINE its first character sits on.
    """
    text = "".join(t["token"] for t in toks)
    line_of_char: list[int] = []
    ln = 0
    for ch in text:
        line_of_char.append(ln)
        if ch == "\n":
            ln += 1
    line_labels = segment_lines(text)
    nlines = len(line_labels)

    pos = 0
    for t in toks:
        s = t["token"]
        li = line_of_char[pos] if pos < len(line_of_char) else nlines - 1
        li = min(max(li, 0), nlines - 1)
        t["module"] = line_labels[li]
        pos += len(s)


# ---- selection-chain SUB-modules ------------------------------------------------
# within the `selection` module, classify each line by the nearest chain opener.
SEL_SUBS = ["track", "photon", "pid", "kinematic_fit", "vertex_fit", "other_sel"]

_SEL_ANCHORS = [
    ("track",         re.compile(r"\bselect_track\b|\bselect_good_charged_track\b|\bselect_charged\b")),
    ("photon",        re.compile(r"\bselect_photon\b|\bselect_isolated_photon\b|\bselect_good_photon\b")),
    ("pid",           re.compile(r"\bpid\b")),
    ("vertex_fit",    re.compile(r"\bsecondary_vertex_fit\b|\bvertex_fit\b")),
    ("kinematic_fit", re.compile(r"\bkinematic_fit\b|\bkalman_kinematic_fit\b|\.fit\b")),
    # everything else in the chain: remove/assign/constrain/mass windows/veto/...
    ("other_sel",     re.compile(r"\.remove\b|\.assign\b|\bconstrain\w*\b|\binvariant_mass_of\b|"
                                 r"\bmass_window\b|\bveto\b|\btag_side\b|\bsignal_side\b|\bmodes\b")),
]


def segment_selection_subs(text: str, line_labels: list[str]) -> list[str | None]:
    """Per-LINE sub-label for lines whose module == 'selection'; None elsewhere.

    Only CHAIN-LEVEL lines (outside any open `{ }` block) may switch the sub-block;
    lines inside a `{ }` body (the cut values, e.g. `chi2_cut 60`,
    `constrain_four_momentum`) inherit the sub-block of the opener. This prevents a
    body keyword like `constrain_four_momentum` from being mistaken for a new
    chain method. Brace depth is tracked across selection lines.
    """
    lines = text.split("\n")
    subs: list[str | None] = []
    cur = "other_sel"
    depth = 0            # net `{` minus `}` seen so far within the selection region
    in_sel = False
    for lab, ln in zip(line_labels, lines):
        if lab != "selection":
            subs.append(None)
            in_sel = False
            depth = 0
            continue
        if not in_sel:
            cur, depth, in_sel = "other_sel", 0, True
        code_part = ln.split("#", 1)[0]
        # only switch when at chain level (no enclosing block open)
        if depth <= 0:
            for sub, rx in _SEL_ANCHORS:
                if rx.search(code_part):
                    cur = sub
                    break
        subs.append(cur)
        depth += code_part.count("{") - code_part.count("}")
        if depth < 0:
            depth = 0
    return subs


def tag_selection_subs(toks: list[dict]) -> None:
    """Annotate each token with 'sel_sub' (one of SEL_SUBS) when module=='selection',
    else None. Requires tag_modules() to have run (uses each token's line)."""
    text = "".join(t["token"] for t in toks)
    line_of_char: list[int] = []
    ln = 0
    for ch in text:
        line_of_char.append(ln)
        if ch == "\n":
            ln += 1
    line_labels = segment_lines(text)
    subs = segment_selection_subs(text, line_labels)
    nlines = len(subs)
    pos = 0
    for t in toks:
        li = line_of_char[pos] if pos < len(line_of_char) else nlines - 1
        li = min(max(li, 0), nlines - 1)
        t["sel_sub"] = subs[li]
        pos += len(t["token"])
