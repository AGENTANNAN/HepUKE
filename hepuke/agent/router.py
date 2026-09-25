"""Entry-level router: classify a question and (optionally) decompose it.

One LLM call before the tool-use loop starts. Three outcomes:

  * `direct`     — the question doesn't need the BESIII knowledge base
                   (small-talk, general knowledge, math, greetings, meta).
                   Skip the subagent fan-out entirely; answer with a plain
                   chat call.
  * `rag`        — single-topic domain question. Run the current subagent
                   loop unchanged.
  * `decompose`  — the question is a bundle. Split into 2-4 independent
                   sub-questions, run each through its own subagent loop
                   (in parallel), then synthesize.

Design notes:
  * The router prompt is tiny (~600 tokens) so its LLM call cost is close
    to noise once prefix caching kicks in.
  * Failure modes fall through to `rag` — a broken router must never lose
    the user's question, only miss a shortcut.
  * The router is a pure classifier; the actual `direct` reply and
    `decompose` synthesis live in `RagAgent.run()`.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Any, Literal, Optional

from hepuke.agent import hep_entities
from hepuke.observability.trace import emit


Mode = Literal["direct", "rag", "decompose"]


@dataclass
class RoutingDecision:
    mode: Mode
    sub_questions: list[str]
    rationale: str = ""

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "sub_questions": list(self.sub_questions),
            "rationale": self.rationale,
        }


ROUTER_SYSTEM_PROMPT = """\
You classify the user's question for a BESIII physics research assistant.

Choose EXACTLY ONE mode:

  "direct"     — the question can be answered without the BESIII knowledge
                 base. Small-talk, greetings ("你好"), simple math ("1+1"),
                 general programming, generic factual (capitals, dates),
                 meta questions about the assistant itself. If in doubt,
                 do NOT pick direct — pick rag.

  "rag"        — a single-topic BESIII question. Anything naming a paper
                 id (BAM-xxxxx), an author, a referee, a memo, a physics
                 process measured at BESIII, a HyperNews thread, or a
                 methodology used in a specific memo. Also: broad
                 "what has X done" over the corpus.

  "decompose"  — the question bundles two or more INDEPENDENT sub-topics
                 that each need their own retrieval. Typical signals:
                   • compares two or more paper ids / authors / topics
                   • lists multiple named entities to fetch metadata for
                   • asks two disconnected things joined by "and" / "以及".
                 Decompose ONLY when the parts can be answered separately;
                 do NOT decompose "who reviewed BAM-00713?" — that's one
                 topic.

Contract — return JSON, no prose, no fences:
  {
    "mode": "direct" | "rag" | "decompose",
    "sub_questions": [ "...", ... ],   // 2-4 items iff mode="decompose"; else []
    "rationale": "<one short sentence>"
  }

For "decompose", each sub_question MUST be self-contained (a downstream
agent will see it in isolation) and MUST preserve every named entity from
the original question.
"""


_JSON_RE = re.compile(r"\{.*\}", re.DOTALL)


def _parse_router_output(raw: str) -> Optional[RoutingDecision]:
    """Try hard to pull a valid RoutingDecision out of the model's reply.

    Accepts fenced JSON, unfenced JSON, or JSON embedded in prose. On any
    schema violation returns None so the caller can fall back to `rag`.
    """
    if not raw:
        return None
    m = _JSON_RE.search(raw)
    if not m:
        return None
    try:
        obj = json.loads(m.group(0))
    except json.JSONDecodeError:
        return None
    if not isinstance(obj, dict):
        return None
    mode = obj.get("mode")
    if mode not in ("direct", "rag", "decompose"):
        return None
    subs = obj.get("sub_questions") or []
    if not isinstance(subs, list):
        subs = []
    subs = [s.strip() for s in subs if isinstance(s, str) and s.strip()]
    if mode == "decompose":
        if not (2 <= len(subs) <= 4):
            # Decompose without a valid split is a routing failure; caller
            # will fall back to rag rather than dispatch nothing.
            return None
    else:
        subs = []
    rationale = obj.get("rationale") or ""
    if not isinstance(rationale, str):
        rationale = ""
    return RoutingDecision(mode=mode, sub_questions=subs, rationale=rationale)


def classify(
    question: str,
    *,
    llm: Any,
    system_prompt: str = ROUTER_SYSTEM_PROMPT,
    fallback: Mode = "rag",
) -> RoutingDecision:
    """Run one LLM call and return a `RoutingDecision`.

    `llm` must expose either `chat_with_tools(messages, tools=[])` (the
    project's LLMClient) — the router doesn't want tool calling, it just
    needs a chat response. We pass an empty `tools` list; the response's
    `.content` is what we parse.

    Any exception, empty response, or schema violation is caught and turned
    into `RoutingDecision(mode=fallback, ...)` — the router must never
    break the pipeline.
    """
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": question},
    ]
    try:
        msg = llm.chat_with_tools(messages, tools=[], tool_choice="none")
    except TypeError:
        # Some fakes don't accept tool_choice; retry without.
        msg = llm.chat_with_tools(messages, tools=[])
    except Exception as e:
        emit("router_error", error=f"{type(e).__name__}: {e}")
        return RoutingDecision(mode=fallback, sub_questions=[], rationale="llm error")

    content = getattr(msg, "content", None) or ""
    decision = _parse_router_output(content)
    if decision is None:
        emit(
            "router_fallback",
            reason="unparseable",
            content_preview=content[:200],
        )
        decision = RoutingDecision(
            mode=fallback, sub_questions=[], rationale="parse fail",
        )

    # HEP-entity override: when the question mentions ≥2 specific states /
    # paper ids / kinematically-tagged processes, force decompose regardless
    # of what the LLM said. This is a hard mechanical rule — the LLM
    # frequently under-decomposes "Zc(3900) and Zc(4020) and Zc(4200)" as a
    # single rag topic even though each state has its own retrieval path.
    # Generic terms (BESIII, J/psi, psi(2S)) never trigger this override.
    scan = hep_entities.extract(question)
    if scan.n_specific() >= 2 and decision.mode != "decompose":
        subs = [f"{question}\n\n(focus: {ent})" for ent in scan.specific[:4]]
        decision = RoutingDecision(
            mode="decompose",
            sub_questions=subs,
            rationale=(
                f"HEP entity override: {scan.n_specific()} specific entities "
                f"→ decompose ({', '.join(scan.specific[:4])})"
            ),
        )
        emit(
            "router_hep_override",
            n_specific=scan.n_specific(),
            entities=scan.specific[:4],
        )

    emit(
        "router_decision",
        mode=decision.mode,
        n_sub=len(decision.sub_questions),
        rationale=decision.rationale[:120],
    )
    return decision


__all__ = ["RoutingDecision", "classify", "ROUTER_SYSTEM_PROMPT"]
