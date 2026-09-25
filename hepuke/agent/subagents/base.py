"""Wrap a nested `RagAgent` as a single `Tool` for the parent agent.

The nested agent runs its own tool-use loop with its own toolset + prompt;
only its final structured return string reaches the parent. This is the
subagent-as-tool pattern: the parent's context stays clean of the child's
tool-call trace, and the child's specialization stays out of the parent's
prompt.
"""
from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Iterable, Optional

from hepuke.agent.loop import RagAgent
from hepuke.agent.tools import Tool, ToolRegistry
from hepuke.observability.trace import emit


_ANSWER_JSON_RE = re.compile(r"\{.*\}", re.DOTALL)

# Doc-id shapes we accept as "grounded" identifiers. Extend as new corpora
# come online. Kept intentionally strict — a hit means the LLM referenced a
# specific document, not a general term.
_DOC_ID_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\bBAM-\d{3,5}\b", re.IGNORECASE),
    re.compile(r"\barXiv:\s*\d{4}\.\d{4,5}\b", re.IGNORECASE),
)

# Tools whose result bodies are trusted sources of doc_ids. A doc_id that
# shows up ONLY in the LLM's answer prose but not in any of these tool
# results is treated as an ungrounded (potentially hallucinated) claim.
_EVIDENCE_TOOL_NAMES: frozenset[str] = frozenset({
    "search_index", "read_summary",
    "grep", "list_memos", "preview_memo", "read_memo",
    "retrieve_vec", "filter_rows", "get_full_doc",
    "milvus_search",
})


def _extract_doc_ids(text: str) -> list[str]:
    """Return de-duplicated doc_ids referenced in `text`, preserving order.

    Case-normalised to upper for BAM ids so they match tool return bodies
    regardless of how the LLM cased them in prose.
    """
    if not text:
        return []
    seen: set[str] = set()
    out: list[str] = []
    for pat in _DOC_ID_PATTERNS:
        for m in pat.finditer(text):
            canon = re.sub(r"\s+", "", m.group(0)).upper()
            if canon not in seen:
                seen.add(canon)
                out.append(canon)
    return out


def _evidence_from_trace(
    trace: list[dict], mentioned_ids: list[str],
) -> list[dict]:
    """Reconstruct evidence records from the tool trace.

    For each doc_id the LLM mentioned in its answer, look for it inside
    the trace's tool `result_preview` blobs. The intersection is what we
    treat as grounded. This is the hard rule: no matter what the LLM
    filled into `evidence`, an id that never appeared in a real tool
    return does not count as evidence.

    Each reconstructed record carries `source: "trace-fallback"` and the
    tool name that produced the hit, so downstream code can tell these
    apart from LLM-authored evidence records.
    """
    if not trace or not mentioned_ids:
        return []
    mentioned_upper = [d.upper() for d in mentioned_ids]
    hits: list[dict] = []
    seen: set[str] = set()
    for entry in trace:
        tool = entry.get("tool")
        if tool not in _EVIDENCE_TOOL_NAMES:
            continue
        preview = (entry.get("result_preview") or "").upper()
        if not preview:
            continue
        for doc_id in mentioned_upper:
            if doc_id in seen:
                continue
            if doc_id in preview:
                hits.append({
                    "source": "trace-fallback",
                    "doc_id": doc_id,
                    "tool": tool,
                    "step": entry.get("step"),
                })
                seen.add(doc_id)
    return hits


_SALVAGE_TOOL_NAMES: frozenset[str] = frozenset({
    # tools whose non-empty return means the corpus actually contained
    # something for this question. Used by the max_steps salvage path.
    "sql_query", "filter_rows", "retrieve_vec", "grep",
    "search_index", "read_summary", "read_memo", "preview_memo",
    "get_full_doc",
})


def _evidence_from_observed(observed: list[dict]) -> list[dict]:
    """Convert loop.py's `observed_hits` into subagent-contract evidence.

    Preferred over `_salvage_from_trace` on the max_steps path because
    the loop extracts hits from the RAW tool result (not the 400-char
    `result_preview`), so a `sql_query` that returned 19 rows contributes
    all 19 doc_ids instead of just the ones that fit in the preview.
    De-dupes by doc_id; keeps a single `row_count`-signal record per
    tool+step so hypernews-style rowset returns without doc_ids still
    surface as a non-dry signal.
    """
    if not observed:
        return []
    out: list[dict] = []
    seen_docs: set[str] = set()
    seen_signals: set[tuple[str, int]] = set()
    for h in observed:
        if not isinstance(h, dict):
            continue
        tool = h.get("tool")
        step = h.get("step")
        doc_id = h.get("doc_id")
        if isinstance(doc_id, str) and doc_id.strip():
            key = doc_id.strip()
            if key in seen_docs:
                continue
            seen_docs.add(key)
            rec = {
                "source": "observed",
                "doc_id": key,
                "tool": tool,
                "step": step,
            }
            score = h.get("score")
            if score is not None:
                rec["score"] = score
            out.append(rec)
            continue
        # No doc_id: row-count signal (hypernews / count queries).
        row_count = h.get("row_count")
        if isinstance(row_count, int) and row_count > 0:
            sig_key = (str(tool), int(step) if isinstance(step, int) else -1)
            if sig_key in seen_signals:
                continue
            seen_signals.add(sig_key)
            out.append({
                "source": "observed",
                "doc_id": None,
                "tool": tool,
                "step": step,
                "row_count": row_count,
                "note": "tool returned rows but agent ran out of steps",
            })
    return out


def _salvage_from_trace(trace: list[dict]) -> list[dict]:
    """When an inner subagent runs out of steps before it can write its
    contract JSON, we still want the parent to know retrieval actually hit
    something. Scan the recorded tool trace for signs of a non-empty
    result — doc_ids inside `result_preview`, or textual hints like
    `row_count`/`n_hits`/`n_matches` above zero — and mint one evidence
    row per hit. `source: "trace-salvage"` marks these so callers can
    tell them apart from genuine LLM-authored evidence.

    Best-effort only: `result_preview` is capped at 400 chars, so this
    reflects "at least one hit exists", not the full result set.
    """
    if not trace:
        return []
    hits: list[dict] = []
    seen_docs: set[str] = set()
    for entry in trace:
        tool = entry.get("tool")
        if tool not in _SALVAGE_TOOL_NAMES:
            continue
        if entry.get("error"):
            continue
        preview = entry.get("result_preview") or ""
        if not preview:
            continue
        for doc_id in _extract_doc_ids(preview):
            if doc_id in seen_docs:
                continue
            seen_docs.add(doc_id)
            hits.append({
                "source": "trace-salvage",
                "doc_id": doc_id,
                "tool": tool,
                "step": entry.get("step"),
            })
        # Even without a parsable doc_id (e.g. hypernews_db count/author
        # rows), signal that the tool returned SOMETHING so the parent
        # doesn't treat the subagent as fully dry.
        preview_lower = preview.lower()
        row_hit = False
        for probe in ("row_count", "n_hits", "n_matches"):
            if probe in preview_lower and probe + '": 0' not in preview_lower \
               and probe + '":0' not in preview_lower:
                row_hit = True
                break
        if row_hit:
            hits.append({
                "source": "trace-salvage",
                "doc_id": None,
                "tool": tool,
                "step": entry.get("step"),
                "note": "tool returned rows but agent ran out of steps",
            })
    return hits


def _coerce_return(raw: dict, *, name: str) -> dict:
    """Normalise the inner agent's `.run()` return into the subagent contract.

    Inner agents are told to return a JSON object matching the contract; if
    they emit prose or a partial object we still surface something the parent
    can act on rather than a black-box error.
    """
    answer_text = (raw.get("answer") or "").strip()
    steps = int(raw.get("steps", 0) or 0)
    trace = raw.get("trace") or []
    terminated = raw.get("terminated")
    observed = raw.get("observed_hits") or []

    parsed: Optional[dict] = None
    if answer_text:
        m = _ANSWER_JSON_RE.search(answer_text)
        if m:
            try:
                candidate = json.loads(m.group(0))
                if isinstance(candidate, dict):
                    parsed = candidate
            except json.JSONDecodeError:
                parsed = None

    if parsed is not None:
        answer = parsed.get("answer")
        if isinstance(answer, str):
            answer = answer.strip() or None
        evidence = parsed.get("evidence") or []
        if not isinstance(evidence, list):
            evidence = []
        # Hard rule: dry ↔ empty evidence. We deliberately ignore whatever
        # `dry` the inner LLM self-reported — it's the same LLM that just
        # produced the (possibly hallucinated) answer, so its self-audit is
        # unreliable. Evidence is the observable ground truth.
        #
        # Fallback: when the LLM forgot to populate `evidence` but its
        # answer names doc_ids we can find in the tool trace, reconstruct
        # evidence from the trace. The intersection of (mentioned in
        # answer) × (returned by a real tool) is deterministic and does
        # not trust the LLM's self-report.
        if not evidence and isinstance(answer, str):
            mentioned = _extract_doc_ids(answer)
            recovered = _evidence_from_trace(trace, mentioned)
            if recovered:
                evidence = recovered
        # Second fallback: max_steps salvage. When the inner loop ran out
        # of steps we treat any real tool hit as at least a "signal"
        # evidence, so the parent sees "we did find something but timed
        # out" instead of a false-negative dry report.
        #
        # Prefer `observed_hits` (accumulated in loop.py from RAW tool
        # results, not truncated previews) over the legacy trace scan.
        if not evidence and terminated == "max_steps":
            salvaged = _evidence_from_observed(observed)
            if not salvaged:
                salvaged = _salvage_from_trace(trace)
            if salvaged:
                evidence = salvaged
        dry = not evidence
        return {
            "answer": answer,
            "evidence": evidence,
            "dry": dry,
            "steps": steps,
            "subagent": name,
        }

    # No JSON — try the same trace-fallback against the prose answer. A
    # subagent that emits prose but names doc_ids the tools actually
    # returned is still grounded; only unsourced prose stays dry.
    if answer_text:
        mentioned = _extract_doc_ids(answer_text)
        recovered = _evidence_from_trace(trace, mentioned)
        if recovered:
            return {
                "answer": answer_text,
                "evidence": recovered,
                "dry": False,
                "steps": steps,
                "subagent": name,
                "format_warning": "inner agent did not return contract JSON",
            }

    # Third fallback: max_steps salvage on the raw-prose path too. Same
    # rationale as above — retrieval did fire, we should not report dry.
    if terminated == "max_steps":
        salvaged = _evidence_from_observed(observed)
        if not salvaged:
            salvaged = _salvage_from_trace(trace)
        if salvaged:
            return {
                "answer": answer_text or None,
                "evidence": salvaged,
                "dry": False,
                "steps": steps,
                "subagent": name,
                "format_warning": (
                    "inner agent hit max_steps before writing an answer; "
                    "evidence recovered from tool trace"
                ),
            }

    return {
        "answer": answer_text or None,
        "evidence": [],
        "dry": True,
        "steps": steps,
        "subagent": name,
        "format_warning": "inner agent did not return contract JSON",
    }


def _merge_candidate_runs(results: list[dict], *, name: str) -> dict:
    """Fold N candidate inner runs into ONE subagent return.

    Rules:
      * `evidence` — union across runs, deduped by doc_id (first seen
        wins so richer LLM-authored records outrank later trace-fallback
        ones).
      * `dry` — True only when EVERY candidate was dry. One grounded
        candidate suffices to call the whole fan-out grounded.
      * `answer` — the first non-dry candidate's answer; if all dry,
        the one with the most steps (proxy for "tried hardest").
      * `steps` — sum across candidates (transparent about total cost).
    """
    if not results:
        return {
            "answer": None, "evidence": [], "dry": True,
            "steps": 0, "subagent": name,
        }

    merged_evidence: list[dict] = []
    seen_doc_ids: set[str] = set()
    total_steps = 0
    any_grounded = False
    first_grounded_answer: Optional[str] = None
    longest_dry_answer: Optional[str] = None
    longest_dry_steps = -1

    for r in results:
        total_steps += int(r.get("steps", 0) or 0)
        r_evidence = r.get("evidence") or []
        for e in r_evidence:
            if not isinstance(e, dict):
                continue
            key = e.get("doc_id")
            if key and key in seen_doc_ids:
                continue
            if key:
                seen_doc_ids.add(key)
            merged_evidence.append(e)
        if not r.get("dry", True):
            any_grounded = True
            if first_grounded_answer is None:
                first_grounded_answer = r.get("answer")
        else:
            steps = int(r.get("steps", 0) or 0)
            if steps > longest_dry_steps:
                longest_dry_steps = steps
                longest_dry_answer = r.get("answer")

    return {
        "answer": first_grounded_answer if any_grounded else longest_dry_answer,
        "evidence": merged_evidence,
        "dry": not any_grounded,
        "steps": total_steps,
        "subagent": name,
    }


def make_subagent_tool(
    *,
    name: str,
    description: str,
    llm: Any,
    tools: Iterable[Tool],
    system_prompt: str,
    max_steps: int = 6,
    memory: Any = None,
    rewriter: Any = None,
    max_candidates: int = 3,
) -> Tool:
    """Package a nested `RagAgent` as a Tool the parent can call.

    Parameters
    ----------
    name, description : identity fields exposed to the parent LLM. `description`
        is the ONLY place the parent learns when to route here — write it as
        one line WHAT + one line WHEN vs siblings.
    llm : same interface as the outer agent expects (must have
        `.chat_with_tools`).
    tools : narrow toolset the subagent is allowed to call. The parent tools
        are intentionally invisible to the subagent — that's the isolation.
    system_prompt : subagent-specific rules. Should end with the contract
        reminder so the inner LLM emits the expected JSON.
    max_steps : inner loop cap. Overridable per-call by the parent.
    memory : optional `SessionMemory`. When provided, a dry return records
        a dead_end file so a sibling / retry can short-circuit. Kept as
        `Any` (duck-typed) so this module doesn't import memory.py — the
        Tool wrapper stays a leaf.
    rewriter : optional multi-candidate query rewriter (LLMQueryRewriter
        or duck-typed equivalent with a `.rewrite_multi(q)` method
        returning objects that expose `.dense_text`). When set, the
        parent's `question` is fanned out into up to `max_candidates`
        English rewrites and each rewrite runs its own inner
        RagAgent.run() concurrently — evidence is merged, dry AND-
        reduces. This is the root fix for "父 LLM 用中文 / 语料英文"
        and "同一 subagent 该并行多个改写" both at the entry point.
        When None, the legacy single-run path is used byte-for-byte.
    max_candidates : cap on how many rewrites we actually launch — even
        if the rewriter produces more.
    """

    from hepuke.observability.trace import (
        _ctx as _trace_ctx,
        set_trace_context,
        clear_trace_context,
    )

    def _call(
        question: str,
        max_steps: Optional[int] = None,
        doc_id_hint: Optional[str] = None,
    ) -> dict:
        # `max_steps` kwarg from the parent is honoured by rebuilding
        # the inner cap inside _run_one — but keeping the closure form
        # was noisy, so we just clamp here and re-close via a nested
        # helper if the parent overrode.
        effective_max_steps = int(max_steps) if max_steps else 6

        def _one_with_cap(prompt: str) -> dict:
            registry = ToolRegistry()
            registry.register_all(list(tools))
            inner = RagAgent(
                llm=llm,
                tools=registry,
                max_steps=effective_max_steps,
                system_prompt=system_prompt,
            )
            parent_ctx = dict(_trace_ctx())
            set_trace_context(**{**parent_ctx, "subagent_name": name})
            try:
                raw = inner.run(prompt)
            finally:
                clear_trace_context()
                if parent_ctx:
                    set_trace_context(**parent_ctx)
            return _coerce_return(raw, name=name)

        # Build the prompt suffix once — the hint is invariant across
        # rewrites.
        def _wrap(q: str) -> str:
            return q if not doc_id_hint else (
                f"{q}\n\n[hint: focus on doc_id `{doc_id_hint}`]"
            )

        # Rewriter fan-out: N candidates, N concurrent inner loops.
        # Each candidate carries its own English `dense_text`; we pass
        # that as the inner user prompt (the sparse tokens don't help
        # the LLM read — they exist for BM25). Concurrency is bounded
        # by max_candidates.
        candidates: list = []
        if rewriter is not None:
            try:
                cands = rewriter.rewrite_multi(question)
                candidates = list(cands[: max(1, int(max_candidates))])
            except Exception:
                candidates = []

        emit(
            "subagent_start", name=name,
            has_hint=bool(doc_id_hint),
            question=question,
            doc_id_hint=doc_id_hint,
            n_candidates=len(candidates) if candidates else 1,
        )

        if not candidates:
            out = _one_with_cap(_wrap(question))
        else:
            # Snapshot the outer trace context so worker threads inherit
            # agent_id / depth. Same pattern as loop._dispatch_all.
            outer_ctx = dict(_trace_ctx())

            def _worker(prompt: str) -> dict:
                if outer_ctx:
                    set_trace_context(**outer_ctx)
                try:
                    return _one_with_cap(prompt)
                finally:
                    clear_trace_context()

            prompts = [_wrap(c.dense_text) for c in candidates]
            results: list[dict] = [{}] * len(prompts)
            with ThreadPoolExecutor(
                max_workers=len(prompts), thread_name_prefix=f"sub-{name}",
            ) as pool:
                futs = {
                    pool.submit(_worker, p): i for i, p in enumerate(prompts)
                }
                for fut in futs:
                    results[futs[fut]] = fut.result()
            out = _merge_candidate_runs(results, name=name)
            # Attach per-candidate audit trail so the parent (and the
            # step-review UI) can see what each rewrite did.
            out["candidates"] = [
                {
                    "label": getattr(c, "label", "") or f"cand{i}",
                    "dense_text": getattr(c, "dense_text", ""),
                    "dry": bool(results[i].get("dry", True)),
                    "evidence_n": len(results[i].get("evidence") or []),
                    "steps": int(results[i].get("steps", 0) or 0),
                }
                for i, c in enumerate(candidates)
            ]

        emit(
            "subagent_end", name=name,
            dry=out["dry"], steps=out["steps"],
            evidence_n=len(out.get("evidence") or []),
            n_candidates=len(candidates) if candidates else 1,
        )
        # Dead-end recording. A dry run means the subagent explored this
        # (subagent, question) pair and produced nothing groundable —
        # future siblings and retries can skip. Guarded on the memory
        # kwarg so tests / callers without a session_root keep the
        # historical behaviour.
        if memory is not None and out.get("dry"):
            try:
                memory.record_dead_end(
                    subagent=name,
                    question=question,
                    reason=f"dry after {out.get('steps', 0)} steps",
                    steps=int(out.get("steps", 0) or 0),
                )
            except Exception:
                # Memory is best-effort; a broken disk write must never
                # take down the agent loop.
                pass
        return out

    # Cap max_steps aggressively — a subagent that wanders 16 steps has
    # already failed. 8 is a soft ceiling; parent can lower it.
    parameters = {
        "type": "object",
        "properties": {
            "question": {
                "type": "string",
                "description": (
                    "Natural-language sub-question rewritten for this "
                    "corpus. Include named entities verbatim."
                ),
            },
            "max_steps": {
                "type": "integer",
                "default": max_steps,
                "minimum": 1,
                "maximum": max(max_steps, 16),
            },
            "doc_id_hint": {
                "type": "string",
                "description": (
                    "Optional doc_id to anchor the search (e.g. 'BAM-00713'). "
                    "Pass when a sibling subagent already found the id."
                ),
            },
        },
        "required": ["question"],
    }
    return Tool(
        name=name,
        description=description,
        parameters=parameters,
        call=_call,
    )


__all__ = [
    "make_subagent_tool",
    "_extract_doc_ids",
    "_evidence_from_trace",
    "_evidence_from_observed",
    "_salvage_from_trace",
    "_merge_candidate_runs",
]
