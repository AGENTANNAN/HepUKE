"""Consume trace JSONL + gold cases → RetrievalReport / AnswerReport.

The runner (`benchmarks/runner.py`) is responsible for:
  1. `configure_trace(enabled=True, file=..., include_context=True)`
  2. For each case, `set_trace_context(case_id=..., suite=..., config=...)`
     so `emit()` records automatically carry the case_id.
This module reads the resulting JSONL back, groups events by `case_id`, and
scores against a gold-case JSONL.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from hepuke.observability.metrics import (
    aggregate_answer,
    aggregate_retrieval,
    exact_match,
    mrr as _mrr,
    ndcg_at_k,
    recall_at_k,
    token_f1,
)
from hepuke.observability.reports import AnswerReport, RetrievalReport
from hepuke.observability.sinks import is_sqlite_path, iter_sqlite_events


# ---------------------------------------------------------------------------
# Case
# ---------------------------------------------------------------------------

@dataclass
class Case:
    case_id: str
    question: str
    answer: list[str] = field(default_factory=list)
    relevant_doc_ids: list[str] = field(default_factory=list)
    suite: str = ""
    meta: dict = field(default_factory=dict)


def _coerce_str_list(v) -> list[str]:
    if v is None:
        return []
    if isinstance(v, str):
        return [v]
    if isinstance(v, (list, tuple)):
        return [str(x) for x in v if x is not None]
    return [str(v)]


def load_cases(path: str | Path) -> list[Case]:
    """Read one Case per JSONL line. Extra fields land in `meta`."""
    p = Path(path)
    cases: list[Case] = []
    with p.open("r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            cid = row.pop("case_id", None) or f"case_{i:06d}"
            q = row.pop("question", "") or ""
            answer = _coerce_str_list(row.pop("answer", None))
            relevant = _coerce_str_list(row.pop("relevant_doc_ids", None))
            suite = row.pop("suite", "") or ""
            meta = row.pop("meta", {}) or {}
            for k, v in row.items():
                meta.setdefault(k, v)
            cases.append(Case(case_id=cid, question=q, answer=answer,
                              relevant_doc_ids=relevant, suite=suite, meta=meta))
    return cases


# ---------------------------------------------------------------------------
# retrieval eval
# ---------------------------------------------------------------------------

def _iter_trace_records(
    path: str | Path,
    *,
    run_id: str | None = None,
) -> Iterable[dict]:
    p = Path(path)
    if is_sqlite_path(p):
        yield from iter_sqlite_events(p, run_id=run_id)
        return
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if run_id is not None and rec.get("run_id") != run_id:
                continue
            yield rec


def _last_retrieve_per_case(
    trace_file: str | Path,
    *,
    run_id: str | None = None,
) -> dict[str, dict]:
    """Return {case_id: last `retrieve` record for that case}."""
    latest: dict[str, dict] = {}
    for rec in _iter_trace_records(trace_file, run_id=run_id):
        if rec.get("kind") != "retrieve":
            continue
        cid = rec.get("case_id")
        if not cid:
            continue
        prev = latest.get(cid)
        if prev is None or rec.get("ts", 0) >= prev.get("ts", 0):
            latest[cid] = rec
    return latest


def _retrieved_doc_ids(record: dict) -> list[str]:
    """Prefer the doc_id list; fall back to parent_ids (extract doc_id prefix)."""
    top = record.get("top_doc_ids")
    if top:
        return [str(d) for d in top if d]
    parents = record.get("top_parent_ids") or []
    out: list[str] = []
    for pid in parents:
        if not pid:
            continue
        s = str(pid)
        # parent_id shape: "{doc_id}::{section_type}::{section_number}"
        out.append(s.split("::", 1)[0])
    return out


def eval_retrieval(
    cases: list[Case],
    trace_file: str | Path,
    *,
    ks: tuple[int, ...] | list[int] = (1, 3, 5, 10),
    suite: str = "",
    config: str = "",
    run_id: str | None = None,
) -> RetrievalReport:
    ks_list = list(ks)
    latest = _last_retrieve_per_case(trace_file, run_id=run_id)
    per_case: list[dict] = []
    for case in cases:
        rec = latest.get(case.case_id)
        retrieved = _retrieved_doc_ids(rec) if rec else []
        rel = case.relevant_doc_ids
        row = {
            "case_id": case.case_id,
            "n_retrieved": len(retrieved),
            "n_relevant": len(rel),
            "recall": {k: recall_at_k(retrieved, rel, k) for k in ks_list},
            "mrr": _mrr(retrieved, rel),
            "ndcg": {k: ndcg_at_k(retrieved, rel, k) for k in ks_list},
            "retrieved": retrieved[:max(ks_list) if ks_list else 10],
            "relevant": rel,
        }
        per_case.append(row)
    agg = aggregate_retrieval(per_case, ks_list)
    return RetrievalReport(
        suite=suite or (cases[0].suite if cases else ""),
        config=config,
        n_cases=len(cases),
        ks=ks_list,
        recall=agg["recall"],
        mrr=agg["mrr"],
        ndcg=agg["ndcg"],
        trace_file=str(trace_file),
        per_case=per_case,
    )


# ---------------------------------------------------------------------------
# answer eval
# ---------------------------------------------------------------------------

def eval_answer(
    cases: list[Case],
    answers: dict[str, str],
    *,
    suite: str = "",
    config: str = "",
) -> AnswerReport:
    per_case: list[dict] = []
    for case in cases:
        pred = answers.get(case.case_id, "") or ""
        em = exact_match(pred, case.answer)
        f1 = token_f1(pred, case.answer)
        per_case.append({
            "case_id": case.case_id,
            "pred": pred,
            "gold": case.answer,
            "em": em,
            "f1": f1,
        })
    agg = aggregate_answer(per_case)
    return AnswerReport(
        suite=suite or (cases[0].suite if cases else ""),
        config=config,
        n_cases=len(cases),
        em=agg["em"],
        f1=agg["f1"],
        per_case=per_case,
    )
