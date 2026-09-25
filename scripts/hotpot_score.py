"""HotpotQA official Ans/Sup/Joint EM+F1 scorer (A3).

Consumes the JSONL written by ``scripts/eval_answer_hotpotqa.py`` — gold
answer and gold supporting_facts are inlined per row, so this script has
no dependency on the raw HotpotQA JSON at scoring time.

Metric implementations are direct ports of the official
``hotpot_evaluate_v1.py`` (Yang et al. 2018, MIT). Every rule mirrors the
reference so the ICLR submission's Ans/Sup/Joint numbers are apples-to-
apples with prior work.

Usage
-----
    python scripts/hotpot_score.py --pred storage/eval_ctrl.jsonl
    python scripts/hotpot_score.py --pred storage/eval_no_ctrl.jsonl storage/eval_ctrl.jsonl
    # compare two runs side by side
    python scripts/hotpot_score.py --pred a.jsonl b.jsonl --labels no_ctrl ctrl
"""

from __future__ import annotations

import argparse
import json
import re
import string
import sys
from collections import Counter
from pathlib import Path


# --------------------------------------------------------------------------- #
# Official metric primitives (verbatim port from hotpot_evaluate_v1.py)       #
# --------------------------------------------------------------------------- #
def normalize_answer(s: str) -> str:
    def remove_articles(text: str) -> str:
        return re.sub(r"\b(a|an|the)\b", " ", text)

    def white_space_fix(text: str) -> str:
        return " ".join(text.split())

    def remove_punc(text: str) -> str:
        exclude = set(string.punctuation)
        return "".join(ch for ch in text if ch not in exclude)

    def lower(text: str) -> str:
        return text.lower()

    return white_space_fix(remove_articles(remove_punc(lower(s or ""))))


def f1_score(prediction: str, ground_truth: str) -> tuple[float, float, float]:
    """Token-overlap F1 with yes/no/noanswer short-circuit.

    Verbatim port of the official rule: if either side is one of the
    "special" answers and they disagree, F1 is zero regardless of token
    overlap. Prevents ``"yes"`` from getting partial credit against
    ``"yesterday"`` and similar false positives.
    """
    npred = normalize_answer(prediction)
    ngold = normalize_answer(ground_truth)
    ZERO = (0.0, 0.0, 0.0)
    if npred in ("yes", "no", "noanswer") and npred != ngold:
        return ZERO
    if ngold in ("yes", "no", "noanswer") and npred != ngold:
        return ZERO
    pred_toks = npred.split()
    gold_toks = ngold.split()
    common = Counter(pred_toks) & Counter(gold_toks)
    num_same = sum(common.values())
    if num_same == 0:
        return ZERO
    precision = num_same / len(pred_toks)
    recall = num_same / len(gold_toks)
    f1 = 2 * precision * recall / (precision + recall)
    return f1, precision, recall


def exact_match_score(prediction: str, ground_truth: str) -> bool:
    return normalize_answer(prediction) == normalize_answer(ground_truth)


def update_answer(
    metrics: dict, prediction: str, gold: str,
) -> tuple[float, float, float]:
    em = exact_match_score(prediction, gold)
    f1, prec, recall = f1_score(prediction, gold)
    metrics["em"] += float(em)
    metrics["f1"] += f1
    metrics["prec"] += prec
    metrics["recall"] += recall
    return float(em), prec, recall


def update_sp(
    metrics: dict, prediction: list, gold: list,
) -> tuple[float, float, float]:
    """Set-comparison over (title, sent_idx) tuples.

    Sp EM = 1 iff prediction set == gold set (no missing, no extra).
    Sp F1 = harmonic mean of prec/recall over the SAME tuples.
    """
    cur = set(map(tuple, prediction or []))
    au = set(map(tuple, gold or []))
    tp = fp = fn = 0
    for e in cur:
        if e in au:
            tp += 1
        else:
            fp += 1
    for e in au:
        if e not in cur:
            fn += 1
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * prec * recall / (prec + recall) if (prec + recall) > 0 else 0.0
    em = 1.0 if (fp + fn == 0) else 0.0
    metrics["sp_em"] += em
    metrics["sp_f1"] += f1
    metrics["sp_prec"] += prec
    metrics["sp_recall"] += recall
    return em, prec, recall


def _empty_metrics() -> dict:
    return {
        "em": 0.0, "f1": 0.0, "prec": 0.0, "recall": 0.0,
        "sp_em": 0.0, "sp_f1": 0.0, "sp_prec": 0.0, "sp_recall": 0.0,
        "joint_em": 0.0, "joint_f1": 0.0, "joint_prec": 0.0, "joint_recall": 0.0,
    }


# --------------------------------------------------------------------------- #
# Score one JSONL                                                             #
# --------------------------------------------------------------------------- #
def score_file(path: Path) -> dict:
    """Read the A2 JSONL and return normalised metrics + counts.

    Rows with ``error`` set (agent crashed after retry) are counted as
    n_error and excluded from the metric denominator — the official
    evaluator penalises missing predictions, but here we want to isolate
    the *controller effect* from upstream API flakes. n_error is reported
    so the reader can see how much of the sample was dropped.

    Rows with a valid ``pred_ans`` but empty ``pred_sp`` (sp picker
    failure) still count toward Ans metrics; their Sp is a legitimate
    zero — the pipeline didn't produce sp facts.
    """
    metrics = _empty_metrics()
    n_ok = n_err = n_joint_ok = 0
    n_sp_err = 0
    rounds_sum = 0
    latency_sum = 0.0
    ctrl_on_flags: set[bool] = set()

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("error"):
                n_err += 1
                continue

            ctrl_on_flags.add(bool(row.get("ctrl_on")))
            # Prefer the compressed short answer written by the A2 short-ans
            # extractor; fall back to the raw long answer for jsonl produced
            # before that field existed. HotpotQA gold is 1-6 tokens, so the
            # short field is what the official metric expects.
            pred_ans = row.get("short_ans") or row.get("pred_ans") or ""
            gold_ans = row.get("gold_ans") or ""
            pred_sp = row.get("pred_sp") or []
            gold_sp = row.get("gold_sp") or []

            em, prec, recall = update_answer(metrics, pred_ans, gold_ans)
            sp_em, sp_prec, sp_recall = update_sp(metrics, pred_sp, gold_sp)

            # Joint = product of ans-side * sp-side (official rule).
            joint_prec = prec * sp_prec
            joint_recall = recall * sp_recall
            joint_f1 = (
                2 * joint_prec * joint_recall / (joint_prec + joint_recall)
                if (joint_prec + joint_recall) > 0 else 0.0
            )
            joint_em = em * sp_em
            metrics["joint_em"] += joint_em
            metrics["joint_f1"] += joint_f1
            metrics["joint_prec"] += joint_prec
            metrics["joint_recall"] += joint_recall

            n_ok += 1
            n_joint_ok += 1
            if row.get("sp_error"):
                n_sp_err += 1
            rounds_sum += int(row.get("rounds") or 0)
            latency_sum += float(row.get("latency_s") or 0.0)

    if n_ok == 0:
        # Return zeros with counts so caller sees the empty run rather than
        # a divide-by-zero crash.
        return {
            "n_ok": 0, "n_error": n_err, "n_sp_error": n_sp_err,
            "ctrl_on": None, "metrics": metrics,
            "avg_rounds": 0.0, "avg_latency_s": 0.0,
        }
    for k in metrics:
        metrics[k] /= n_ok
    return {
        "n_ok": n_ok,
        "n_error": n_err,
        "n_sp_error": n_sp_err,
        "ctrl_on": (next(iter(ctrl_on_flags)) if len(ctrl_on_flags) == 1 else "mixed"),
        "metrics": metrics,
        "avg_rounds": rounds_sum / n_ok,
        "avg_latency_s": latency_sum / n_ok,
    }


# --------------------------------------------------------------------------- #
# Presentation                                                                #
# --------------------------------------------------------------------------- #
_METRIC_ROWS = [
    ("Ans EM",   "em"),
    ("Ans F1",   "f1"),
    ("Ans Prec", "prec"),
    ("Ans Rec",  "recall"),
    ("Sup EM",   "sp_em"),
    ("Sup F1",   "sp_f1"),
    ("Sup Prec", "sp_prec"),
    ("Sup Rec",  "sp_recall"),
    ("Joint EM", "joint_em"),
    ("Joint F1", "joint_f1"),
    ("Joint Prec","joint_prec"),
    ("Joint Rec", "joint_recall"),
]


def _print_report(paths: list[Path], labels: list[str], results: list[dict]) -> None:
    # Header row: metric | label1 | label2 | ...
    col_w = max(10, max((len(l) for l in labels), default=10))
    print()
    print(f"{'':<12}" + "".join(f"{l:>{col_w+2}}" for l in labels))
    print("-" * (12 + (col_w + 2) * len(labels)))
    for pretty, key in _METRIC_ROWS:
        cells = []
        for r in results:
            v = r["metrics"].get(key, 0.0) * 100.0
            cells.append(f"{v:>{col_w+2}.2f}")
        print(f"{pretty:<12}" + "".join(cells))
    print("-" * (12 + (col_w + 2) * len(labels)))
    for pretty, key in [
        ("n_ok",         "n_ok"),
        ("n_error",      "n_error"),
        ("n_sp_error",   "n_sp_error"),
        ("avg_rounds",   "avg_rounds"),
        ("avg_latency",  "avg_latency_s"),
        ("ctrl_on",      "ctrl_on"),
    ]:
        cells = []
        for r in results:
            v = r.get(key)
            if isinstance(v, float):
                cells.append(f"{v:>{col_w+2}.2f}")
            else:
                cells.append(f"{str(v):>{col_w+2}}")
        print(f"{pretty:<12}" + "".join(cells))
    print()
    for p, r in zip(paths, results):
        print(f"  {p}  ({r['n_ok']} scored, {r['n_error']} errored)")
    print()


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pred", type=Path, nargs="+", required=True,
                   help="One or more prediction JSONLs (A2 output). "
                        "Multiple files are scored and printed side-by-side.")
    p.add_argument("--labels", nargs="+", default=None,
                   help="Column labels matching --pred order. Defaults "
                        "to the file stems.")
    p.add_argument("--json", action="store_true",
                   help="Also emit the full per-file metric dict as JSON "
                        "to stdout after the table (for downstream tools).")
    return p.parse_args()


def main() -> int:
    args = _parse_args()

    labels = args.labels or [p.stem for p in args.pred]
    if len(labels) != len(args.pred):
        print(f"error: --labels ({len(labels)}) must match --pred ({len(args.pred)})",
              file=sys.stderr)
        return 2

    results = []
    for path in args.pred:
        if not path.exists():
            print(f"error: {path} not found", file=sys.stderr)
            return 2
        results.append(score_file(path))

    _print_report(args.pred, labels, results)

    if args.json:
        payload = {
            label: {"path": str(path), **res}
            for label, path, res in zip(labels, args.pred, results)
        }
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
