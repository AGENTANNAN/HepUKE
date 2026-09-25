"""OpenAI-tool-calling agent loop (Claude-Code independent).

Usage:
    from hepuke.agent import RagAgent, ToolRegistry
    tools = ToolRegistry()
    tools.register_all(build_default_tools(retriever, memo_root))
    agent = RagAgent(llm=llm_client, tools=tools, max_steps=8)
    answer = agent.run("What are the systematic uncertainties in ...?")
"""
from __future__ import annotations

import json
import time
import uuid
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Optional

from hepuke.agent.prompts import SYSTEM_PROMPT
from hepuke.agent.tools import Tool, ToolRegistry
from hepuke.config import ConfidenceConfig
from hepuke.observability.trace import (
    _ctx as _trace_ctx,
    clear_trace_context,
    emit,
    set_trace_context,
)


ASK_USER_TOOL_NAME = "ask_user"

# One-line progress probe fired between the fan-out and the step-review
# pause when `enable_perturn_summary` is on. The parent LLM sees all tool
# results already appended to session.messages; this extra user turn asks
# it to summarise progress in 1–2 sentences WITHOUT emitting a final
# answer. `tools=[]` + `tool_choice="none"` guarantee no tool_call so we
# never fork off another dispatch on the probe.
PROGRESS_PROBE_PROMPT = (
    "请用 1–2 句中文报告当前搜索进展：**已经获得了什么关键信息**，"
    "**还差什么**，**下一步倾向于做什么**。不要给最终答案，也不要复述"
    "原始工具输出，只做进度总结。"
)
ASK_USER_TOOL_SCHEMA: dict = {
    "type": "function",
    "function": {
        "name": ASK_USER_TOOL_NAME,
        "description": (
            "Ask the human user a single clarifying question. Use ONLY when "
            "you cannot make progress without input the user has not given "
            "— e.g. an ambiguous entity name, a missing constraint, or a "
            "choice between valid interpretations. Do NOT use for general "
            "clarification you can infer from context, and do NOT combine "
            "with other tool calls in the same turn (siblings will be "
            "deferred until the user replies). The user's reply is delivered "
            "as this tool's return value."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "question": {
                    "type": "string",
                    "description": (
                        "One concise question, in the user's language. "
                        "Include just enough context that the user can "
                        "answer without re-reading the whole thread."
                    ),
                },
            },
            "required": ["question"],
        },
    },
}


@dataclass
class PendingCall:
    """One tool_call captured mid-turn for deferred dispatch on resume."""
    id: str
    name: str
    arguments: str  # raw JSON string, as the LLM emitted


@dataclass
class AgentSession:
    """Serializable snapshot of an in-flight RagAgent run.

    Persisted between turns when `ask_user` suspends the loop. Everything
    here is JSON-friendly so a future disk-backed session store can drop in
    without touching the loop.
    """
    agent_id: str
    messages: list[dict] = field(default_factory=list)
    trace: list[dict] = field(default_factory=list)
    call_counter: dict[str, int] = field(default_factory=dict)
    tokens_spent: int = 0
    step: int = 0
    limit: int = 8
    routing: Optional[dict] = None
    # Stable id that ties every retrieval / raw-result row across MULTIPLE
    # /agent/chat turns of the same session (ask_user + step-review can
    # split one logical run into N HTTP requests). Defaults to agent_id
    # so in-process run()s still get a non-null id; the fastapi layer may
    # overwrite this with the store's session_id.
    session_id: str = ""
    # Set when the loop is paused waiting for a user reply.
    pending_ask: Optional[dict] = None  # {"tool_call_id": str, "question": str}
    # Non-ask_user tool_calls emitted in the SAME turn as the ask_user.
    # OpenAI protocol demands every tool_call have a matching tool-role
    # message, so we defer these until resume() and dispatch them then.
    pending_deferred: list[dict] = field(default_factory=list)
    # Step-review pause: set when enable_step_review=True and a fan-out
    # round just finished. Resume via resume_from_review(session, reply).
    pending_review: Optional[dict] = None
    # Real-time hit accumulator: every tool_call that returns structured
    # results contributes zero or more {tool, step, doc_id?, row_count?,
    # score?} entries via `_extract_hits`. Populated inside dispatch,
    # BEFORE the raw result is serialised (so we're not bounded by the
    # 400-char `result_preview` in the trace). Consumed by the subagent
    # base.py salvage path when the inner loop hits max_steps.
    observed_hits: list[dict] = field(default_factory=list)
    # Per-round answer posterior p_t (top-B softmax over first content-token
    # logprobs). Populated only when the confidence controller is active
    # AND the LLM client supports the `logprobs=True` kwarg. Length matches
    # the number of parent LLM calls made during this session (one entry
    # per step). Entries may be empty lists on rounds where the server did
    # not return logprobs.
    posteriors: list[list[float]] = field(default_factory=list)


class RagAgent:
    """OpenAI function-calling agent loop.

    Parameters
    ----------
    llm : object with `.chat_with_tools(messages, tools, tool_choice)` returning
          an OpenAI-shape `ChatCompletionMessage` (has `.content` and
          `.tool_calls`).
    tools : `ToolRegistry` populated with agent-callable tools.
    max_steps : hard cap on tool-call rounds (protects against runaway loops).
    system_prompt : optional override of the default SYSTEM_PROMPT.
    parallel_tool_calls : when the assistant emits N tool_calls in one turn,
        dispatch them concurrently via a thread pool instead of one after
        the other. Wall-clock drops to max(latencies) instead of the sum —
        the win for a subagent fan-out (three subagents × ~60s each) is
        real. Set False to serialise (deterministic for debugging).
    max_parallel : upper bound on concurrent tool dispatches. Defaults to
        min(len(tool_calls), 8). Only used when `parallel_tool_calls=True`.
    """

    def __init__(
        self,
        *,
        llm: Any,
        tools: ToolRegistry,
        max_steps: int = 8,
        system_prompt: str | None = None,
        parallel_tool_calls: bool = True,
        max_parallel: int = 8,
        route: bool = False,
        direct_system_prompt: str | None = None,
        decompose_max_parallel: int = 4,
        repeated_tool_call_threshold: int = 3,
        token_budget: int = 0,
        token_budget_enabled: bool = False,
        enable_ask_user: bool = False,
        enable_step_review: bool = False,
        enable_perturn_summary: bool = True,
        memory: Any = None,
        confidence: ConfidenceConfig | None = None,
        round_observer: Any = None,
        collect_logprobs_only: bool = False,
        force_first_tool: Optional[str] = None,
        force_search_after_first: bool = False,
        posterior_llm: Any = None,
    ) -> None:
        self.llm = llm
        # Dedicated posterior client for shadow probes. Defaults to the main
        # llm for backward compatibility so existing tests / callers keep
        # working; production callers plug in a cheap model (e.g. gpt-4.1-mini)
        # so token-logprob harvesting does not require the main agent to run
        # on a logprob-friendly model.
        self.posterior_llm = posterior_llm or llm
        self.tools = tools
        self.max_steps = max_steps
        self.system_prompt = system_prompt or SYSTEM_PROMPT
        self.parallel_tool_calls = parallel_tool_calls
        self.max_parallel = max(1, int(max_parallel))
        self.route = route
        self.direct_system_prompt = direct_system_prompt or (
            "You are a helpful assistant. Answer the user's question "
            "directly and concisely. No knowledge base is available."
        )
        self.decompose_max_parallel = max(1, int(decompose_max_parallel))
        # 0 disables the repeated-call guard; otherwise the Nth+1 identical
        # (tool_name, canonical_args) call in the same run() gets its tool
        # result prefixed with a WARNING banner so the LLM can see it's
        # spinning. This is a hard signal, NOT a prompt hint.
        self.repeated_tool_call_threshold = max(0, int(repeated_tool_call_threshold))
        # Cap on total completion tokens across one run(). When enabled and
        # the accumulated count crosses the threshold, the next step is
        # skipped and the loop returns terminated=token_budget. Hard stop —
        # no prompt hint. Set to 0 (or _enabled=False) to disable.
        self.token_budget = max(0, int(token_budget))
        self.token_budget_enabled = bool(token_budget_enabled) and self.token_budget > 0
        # ask_user: when True, `ask_user` schema is exposed as a tool at
        # depth==0 (never to subagents, never during decompose sub-runs).
        # Emitting an ask_user tool_call suspends the loop; the caller
        # resume()s once the user has replied.
        self.enable_ask_user = bool(enable_ask_user)
        # Step-review: when True, after EVERY tool-call fan-out finishes
        # (all tool-role messages appended), the loop pauses with a
        # `status='review'` summary of what each sub returned, letting the
        # user decide continue / stop / redirect before the next LLM turn.
        # Interactive path only — opt-in per RagAgent instance so the plain
        # /agent/ask endpoint keeps its non-interactive semantics.
        self.enable_step_review = bool(enable_step_review)
        # Per-turn progress summary: when True AND step_review is on, a
        # tool-less LLM probe fires between the fan-out and the review
        # pause, asking the parent LLM to summarise progress in 1–2
        # sentences. The text lands on `pending_review["progress"]` and
        # is surfaced to the caller (see fastapi_app / examples/chat).
        # Off cheaply skips the extra LLM call — one hop per step is not
        # free.
        self.enable_perturn_summary = bool(enable_perturn_summary)
        # Confidence-based retrieval controller (see hepuke/agent/confidence).
        # When config.enabled=True, a per-round controller runs at the end of
        # every _advance iteration and can force an early stop. When False
        # (default), the loop keeps its fixed max_steps behaviour untouched —
        # no controller is constructed and no per-round hook fires.
        self.confidence_config = confidence
        self._confidence_controller = None
        if confidence is not None and confidence.enabled:
            # Local import to avoid the confidence subpackage becoming a hard
            # import-time dependency of loop.py when the feature is off.
            from hepuke.agent.confidence import ConfidenceController

            self._confidence_controller = ConfidenceController(confidence)
        # Round observer: optional callback fired at the end of every
        # _advance iteration (AFTER the confidence controller has run, in
        # the round that just finished dispatching tools). Signature:
        #   ``observer(step_index: int, session: AgentSession) -> None``.
        # Used by the T7 replay collector (hepuke/agent/confidence/
        # collect_predictor_data.py) to snapshot per-round posteriors +
        # observed doc ids without touching the loop's return contract.
        # Exceptions from the observer are logged but never propagate.
        self.round_observer = round_observer
        # Collect-only mode: request logprobs from the LLM and stash
        # posteriors on the session WITHOUT constructing a controller.
        # Used by the T7 replay collector to record unshortened
        # trajectories. Ignored when a real controller is set (in that
        # case posteriors are already collected as part of the pipeline).
        # `top_b` for this mode reuses ConfidenceConfig.top_b via
        # `_effective_top_b`; when no confidence config is passed, defaults
        # to 8 (paper's setting).
        self._collect_logprobs_only = bool(collect_logprobs_only)
        # Force the model's very first tool_choice to this tool name (OpenAI
        # `tool_choice={"type":"function","function":{"name":...}}`). Round 2
        # onward switches back to "auto" and the forced tool is FILTERED OUT
        # of the tool_schemas so the LLM cannot call it again — matches the
        # "rewrite once at round 1, then let the agent plan freely" contract.
        # None disables both behaviours (backwards compatible).
        self.force_first_tool = (force_first_tool or None)
        # force_search_after_first: on every round AFTER the first, force the
        # "search" tool via tool_choice so the LLM cannot self-stop. Used by
        # the "full rounds" ablation arm — the loop always runs to max_steps
        # (the cost ceiling), then the max_steps path synthesises the answer.
        # Off by default; orthogonal to force_first_tool (round 0 rewrite).
        self.force_search_after_first = bool(force_search_after_first)
        # Session-scoped file memory. When set, load_context() is prepended
        # to the system prompt on every OUTER run (depth==0) and the first
        # question of a fresh session is recorded as user_intent. Subagent
        # re-entries use the memory the parent already wired into them via
        # make_subagent_tool(memory=...), not this instance-level attr —
        # keep them separate so decompose sub-runs don't double-write.
        self.memory = memory

    def run(
        self,
        question: str,
        *,
        max_steps: Optional[int] = None,
        session_id: Optional[str] = None,
    ) -> dict:
        """Run the agent loop, return `{'answer': str, 'trace': [...], 'steps': int}`.

        The `trace` field is a Python-side transcript of tool calls (name,
        arguments summary, result summary). Emits are also written to the
        JSONL trace sink if configured.

        Nesting: every run() mints a fresh `agent_id` (an 8-char hex slug)
        and stamps it into thread-local trace context along with the parent
        agent id and a monotonically-increasing depth. A subagent spawned
        via `make_subagent_tool` inherits the parent context first, then
        this method overwrites `agent_id` with its own — so every emitted
        record in the subagent's transcript carries both `parent_agent_id`
        and the child `agent_id`, which is what makes the JSONL/SQLite
        trace disentanglable.

        Routing: when `route=True` and this call is at depth 0 (outermost
        agent, not a subagent-as-tool re-entry), a cheap classifier runs
        first to decide between `direct` (skip tools), `rag` (this loop),
        and `decompose` (split into parallel sub-questions + synthesize).
        Subagent re-entries always skip routing — they've been assigned a
        specific corpus already.
        """
        limit = max_steps or self.max_steps
        parent_ctx = dict(_trace_ctx())
        agent_id = uuid.uuid4().hex[:8]
        depth = int(parent_ctx.get("depth", -1)) + 1
        parent_id = parent_ctx.get("agent_id")

        new_ctx = {
            **parent_ctx,
            "agent_id": agent_id,
            "depth": depth,
        }
        if parent_id:
            new_ctx["parent_agent_id"] = parent_id
        set_trace_context(**new_ctx)

        # Session memory intent-write lives in `_new_session` so both
        # `run()` and the fastapi chat path (which calls `_new_session`
        # directly) behave identically. Nothing to do here.

        try:
            if self.route and depth == 0:
                return self._run_routed(question, limit=limit)
            return self._run_inner(question, limit=limit, session_id=session_id)
        finally:
            # Restore the parent context for the calling thread. clear +
            # re-apply is simpler than tracking removed keys.
            clear_trace_context()
            if parent_ctx:
                set_trace_context(**parent_ctx)

    # ------------------------------------------------------------------
    # routing (only fires at depth == 0)
    # ------------------------------------------------------------------

    def _run_routed(self, question: str, *, limit: int) -> dict:
        """Classify the question, then take the matching path.

        On any router failure we fall through to the standard rag loop —
        the router must never lose the user's question, only miss a
        shortcut.
        """
        # Imported inside the method to avoid a circular import at module
        # load time (router → observability.trace → agent.loop is fine,
        # but keeping the import local documents the dependency direction).
        from hepuke.agent.router import classify

        decision = classify(question, llm=self.llm)

        if decision.mode == "direct":
            return self._run_direct(question, decision=decision)
        if decision.mode == "decompose":
            return self._run_decompose(question, decision=decision, limit=limit)
        # rag — fall through to the standard loop.
        return self._run_inner(question, limit=limit, routing=decision)

    def _run_direct(self, question: str, *, decision) -> dict:
        """Answer without touching the subagent tools.

        Kept intentionally minimal — one chat call with the direct system
        prompt, no tool schemas, no loop. Emits `agent_end` so the trace
        SQL still counts one turn.
        """
        emit("agent_start", question=question, max_steps=0, mode="direct")
        t0 = time.time()
        msg = self.llm.chat_with_tools(
            [
                {"role": "system", "content": self.direct_system_prompt},
                {"role": "user", "content": question},
            ],
            tools=[],
            tool_choice="none",
        )
        dt = round((time.time() - t0) * 1000, 1)
        answer = (getattr(msg, "content", None) or "").strip()
        emit("agent_end", steps=0, latency_ms=dt, answer_len=len(answer), mode="direct")
        return {
            "answer": answer,
            "trace": [],
            "steps": 0,
            "routing": decision.to_dict(),
        }

    def _run_decompose(self, question: str, *, decision, limit: int) -> dict:
        """Run each sub-question through its own RAG loop, then synthesize.

        The child runs happen on a thread pool (wall-clock ≈ max, not sum).
        Trace context is re-hydrated across the thread boundary the same
        way parallel tool dispatch does it, so each child inherits the
        outer agent's `agent_id` as `parent_agent_id`.
        """
        emit(
            "agent_start", question=question, max_steps=limit,
            mode="decompose", n_sub=len(decision.sub_questions),
        )
        # Per-call snapshot — see _dispatch_all() for why this must not
        # live on `self`. When two sub-runs execute concurrently on the
        # same RagAgent instance, an instance-level snapshot would race
        # and workers would inherit the wrong parent_agent_id.
        ctx_snapshot = dict(_trace_ctx())
        subs = list(decision.sub_questions)
        workers = min(len(subs), self.decompose_max_parallel)

        def _run_one(sub_q: str) -> dict:
            if ctx_snapshot:
                set_trace_context(**ctx_snapshot)
            try:
                # A nested run() bumps depth automatically; its own routing
                # is skipped because depth != 0 after the bump.
                return self.run(sub_q, max_steps=limit)
            finally:
                clear_trace_context()

        results: list[dict] = [{}] * len(subs)
        with ThreadPoolExecutor(
            max_workers=workers, thread_name_prefix="decompose",
        ) as pool:
            futures = {pool.submit(_run_one, q): i for i, q in enumerate(subs)}
            for fut in futures:
                results[futures[fut]] = fut.result()

        # Synthesis pass: fold sub-answers into one answer.
        joined = "\n\n".join(
            f"### Sub-question {i+1}: {q}\n{results[i].get('answer','') or '(no answer)'}"
            for i, q in enumerate(subs)
        )
        synth_prompt = (
            "You just received answers to sub-questions decomposed from a "
            "single user question. Combine them into one coherent, "
            "grounded answer. Preserve every doc_id / paper_id / file "
            "reference that appeared. Do not invent facts absent from the "
            "sub-answers."
        )
        t0 = time.time()
        msg = self.llm.chat_with_tools(
            [
                {"role": "system", "content": synth_prompt},
                {
                    "role": "user",
                    "content": (
                        f"Original question: {question}\n\n"
                        f"Sub-answers:\n{joined}"
                    ),
                },
            ],
            tools=[],
            tool_choice="none",
        )
        dt = round((time.time() - t0) * 1000, 1)
        answer = (getattr(msg, "content", None) or "").strip()

        merged_trace: list[dict] = []
        for i, r in enumerate(results):
            for entry in r.get("trace", []) or []:
                merged_trace.append({**entry, "sub_q_idx": i})
        emit(
            "agent_end", steps=sum(r.get("steps", 0) for r in results),
            latency_ms=dt, answer_len=len(answer), mode="decompose",
        )
        return {
            "answer": answer,
            "trace": merged_trace,
            "steps": sum(r.get("steps", 0) for r in results),
            "routing": decision.to_dict(),
            "sub_results": [
                {"question": subs[i], "answer": results[i].get("answer", ""),
                 "steps": results[i].get("steps", 0)}
                for i in range(len(subs))
            ],
        }

    def _run_inner(
        self, question: str, *, limit: int, routing=None,
        session_id: Optional[str] = None,
    ) -> dict:
        """Legacy synchronous entry point.

        Runs the loop to completion (final answer OR max_steps OR token
        budget). If `enable_ask_user=True` and the loop suspends waiting
        for a user reply, returns a `{'status': 'question', 'session': ...}`
        pending result — the caller (typically the /agent/chat endpoint)
        is responsible for round-tripping to the user and calling resume().
        """
        session = self._new_session(
            question, limit=limit, routing=routing, session_id=session_id,
        )
        self.last_session = session
        return self._advance(session)

    # ------------------------------------------------------------------
    # public suspend/resume API
    # ------------------------------------------------------------------

    def resume(self, session: AgentSession, user_reply: str) -> dict:
        """Continue a paused run() with the user's reply.

        Injects the reply as the tool-role response for the pending
        `ask_user` call, dispatches any deferred siblings, then re-enters
        the loop. Raises ValueError if the session is not actually paused.
        """
        if not session.pending_ask:
            raise ValueError("resume() called on a session with no pending ask_user")
        pending_id = session.pending_ask["tool_call_id"]
        # 1. Answer the ask_user tool_call with the user's reply.
        session.messages.append({
            "role": "tool",
            "tool_call_id": pending_id,
            "content": user_reply,
        })
        emit(
            "agent_ask_user_resume",
            step=session.step,
            reply_len=len(user_reply or ""),
        )
        # 2. Dispatch the deferred siblings — their tool_calls are already
        # in the assistant echo, so the protocol requires a tool-role
        # message per id. Reuse the standard dispatch path so trace + repeat
        # detection stay consistent.
        if session.pending_deferred:
            deferred = [_DictToolCall(d) for d in session.pending_deferred]
            counter = Counter(session.call_counter)
            results = self._dispatch_all(
                deferred, step=session.step, trace=session.trace,
                observed=session.observed_hits,
            )
            for tc, result_json in zip(deferred, results):
                content = self._maybe_warn_repeat(
                    tc, result_json, counter, step=session.step,
                )
                session.messages.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "content": content,
                })
            session.call_counter = dict(counter)
        # 3. Clear pending state and continue the loop.
        session.pending_ask = None
        session.pending_deferred = []
        session.step += 1
        return self._advance(session)

    def resume_from_review(self, session: AgentSession, user_reply: str) -> dict:
        """Continue a step-review pause.

        The user's reply is a directive, not an answer to a tool_call:
          * "" / "continue" / "c" / "y" / "继续" / "next" → resume the loop
            unchanged (parent LLM sees the fan-out results and decides).
          * "stop" / "s" / "quit" / "停" / "abort" → terminate immediately,
            surface whatever content the assistant last emitted (typically
            empty in the middle of a fan-out) plus a terminated marker.
          * anything else → treat as a NEW user instruction: append a plain
            user-role message so the parent LLM's next turn can react to it.
        Raises ValueError if the session isn't paused at a review.
        """
        if not session.pending_review:
            raise ValueError(
                "resume_from_review() called on a session with no pending review"
            )
        reply = (user_reply or "").strip()
        emit(
            "agent_step_review_resume",
            step=session.step,
            reply_len=len(reply),
        )
        session.pending_review = None

        stop_tokens = {"stop", "s", "quit", "abort", "停", "终止", "结束"}
        continue_tokens = {"", "continue", "c", "y", "yes", "next", "继续", "ok"}

        if reply.lower() in stop_tokens:
            emit("agent_end", steps=session.step, terminated="user_stop")
            out = {
                "status": "answer",
                "answer": "[agent] 用户在中途终止了搜索。",
                "trace": session.trace,
                "steps": session.step,
                "terminated": "user_stop",
                "observed_hits": list(session.observed_hits),
            }
            if session.routing is not None:
                out["routing"] = session.routing
            return out

        if reply.lower() not in continue_tokens:
            # Redirect: inject the reply as a normal user message so the
            # parent LLM sees it on its next turn.
            session.messages.append({"role": "user", "content": reply})

        return self._advance(session)

    # ------------------------------------------------------------------
    # session construction + main loop
    # ------------------------------------------------------------------

    def _new_session(
        self, question: str, *, limit: int, routing=None,
        session_id: Optional[str] = None,
    ) -> AgentSession:
        # Persist first-question intent AND read back the memory block in
        # the same place so both the run()-path and the fastapi chat-path
        # (which builds a session directly via _new_session/_advance,
        # bypassing run()) get identical memory behaviour. Guard on
        # trace-context depth so decompose sub-runs (self.run(sub_q))
        # don't clobber the outer intent with a sub-question. Best-
        # effort: a broken write must not stop the loop.
        _depth = int(_trace_ctx().get("depth", 0) or 0)
        if self.memory is not None and _depth == 0:
            try:
                if not self.memory.read_user_intent():
                    self.memory.write_user_intent(question)
            except Exception:
                pass
        # Prepend session-memory context (user_intent + dead-ends) to the
        # system prompt when a SessionMemory is wired. The block is a
        # separate ##-titled section so the LLM can read it without the
        # base prompt's rules bleeding into it, and empty when the
        # session has nothing recorded yet — no separator gymnastics
        # needed in the caller.
        system_content = self.system_prompt
        if self.memory is not None:
            try:
                mem_ctx = self.memory.load_context()
            except Exception:
                mem_ctx = ""
            if mem_ctx:
                system_content = mem_ctx + "\n" + self.system_prompt
        session = AgentSession(
            agent_id=uuid.uuid4().hex[:8],
            messages=[
                {"role": "system", "content": system_content},
                {"role": "user", "content": question},
            ],
            limit=limit,
            routing=routing.to_dict() if routing is not None else None,
        )
        # Default session_id to the fresh agent_id. Callers (e.g. fastapi
        # SessionStore) may overwrite with a stable id that spans multiple
        # HTTP turns of the same paused/resumed run.
        session.session_id = session_id or session.agent_id
        emit("agent_start", question=question, max_steps=limit)
        return session

    def _tool_schemas_for_session(self) -> list[dict]:
        schemas = self.tools.as_openai_tools()
        # ask_user is only meaningful at the top level and only when the
        # host (an interactive endpoint) has opted in. Subagents / decompose
        # sub-runs get plain `route=False` re-entries via run(), which set
        # enable_ask_user themselves — we don't need to check depth here.
        if self.enable_ask_user:
            schemas = schemas + [ASK_USER_TOOL_SCHEMA]
        return schemas

    def _advance(self, session: AgentSession) -> dict:
        """Core loop. Runs until final answer, max_steps, token budget,
        or an ask_user pause. All state lives on `session`.
        """
        # Stamp the session id into the trace context so every emit fired
        # from this run (and from worker threads that snapshot ctx) carries
        # it — retrievals / tool_results_raw pull it from the record.
        if session.session_id:
            set_trace_context(session_id=session.session_id)
        tool_schemas = self._tool_schemas_for_session()
        call_counter: Counter[str] = Counter(session.call_counter)

        while session.step < session.limit:
            step = session.step
            # Token-budget backstop. Check BEFORE the LLM call so we don't
            # pay for the step that would push us over — but skip on step 0
            # so a very small budget doesn't kill the run before its first
            # request. Hard runtime cap; the LLM cannot see it.
            if (
                self.token_budget_enabled
                and step > 0
                and session.tokens_spent >= self.token_budget
            ):
                emit(
                    "agent_end", steps=step,
                    terminated="token_budget",
                    tokens_spent=session.tokens_spent,
                    token_budget=self.token_budget,
                )
                session.call_counter = dict(call_counter)
                out = {
                    "status": "answer",
                    "answer": (
                        f"[agent] token budget exceeded "
                        f"({session.tokens_spent}/{self.token_budget}); "
                        f"stopping without final answer."
                    ),
                    "trace": session.trace,
                    "steps": step,
                    "terminated": "token_budget",
                    "tokens_spent": session.tokens_spent,
                    "observed_hits": list(session.observed_hits),
                }
                if session.routing is not None:
                    out["routing"] = session.routing
                return out

            t0 = time.time()
            llm_kwargs: dict = {}
            # force_first_tool contract: round 1 is locked to the named
            # tool via tool_choice; round 2+ runs auto AND removes that
            # tool from the visible schemas so the LLM can't call it again.
            # Applies at depth==0 only (subagents keep the classic path).
            step_tool_schemas = tool_schemas
            step_tool_choice: Any = "auto"
            if self.force_first_tool:
                if step == 0:
                    step_tool_choice = {
                        "type": "function",
                        "function": {"name": self.force_first_tool},
                    }
                else:
                    step_tool_schemas = [
                        s for s in tool_schemas
                        if (s.get("function") or {}).get("name")
                        != self.force_first_tool
                    ]
            # force_search_after_first: lock rounds 2..max_steps to `search`
            # so the LLM always retrieves and cannot end the loop by omitting
            # tool_calls (self-stop). See __init__ docstring.
            if self.force_search_after_first and step > 0:
                step_tool_choice = {
                    "type": "function",
                    "function": {"name": "search"},
                }
            want_posterior = (
                self._confidence_controller is not None
                or self._collect_logprobs_only
            )
            # Main call NEVER requests logprobs — deepseek/HepAI returns 400
            # when the primary tool-choice=auto call has logprobs=True. The
            # shadow probe below runs on self.posterior_llm (cheap model) and
            # is the only place we actually harvest content-token posteriors.
            # `want_posterior` still gates the shadow probe fire below.
            try:
                msg = self.llm.chat_with_tools(
                    session.messages, tools=step_tool_schemas,
                    tool_choice=step_tool_choice,
                    **llm_kwargs,
                )
            except TypeError:
                # Legacy client without logprobs kwargs: retry without.
                msg = self.llm.chat_with_tools(
                    session.messages, tools=step_tool_schemas,
                    tool_choice=step_tool_choice,
                )
            dt = round((time.time() - t0) * 1000, 1)
            step_tokens = getattr(msg, "_hepuke_completion_tokens", None) or 0
            session.tokens_spent += int(step_tokens)

            # Snapshot the round's posterior (may be empty when logprobs
            # weren't returned). One entry per _advance iteration, in order.
            if want_posterior:
                posterior = getattr(msg, "_hepuke_posterior", None) or []
                session.posteriors.append(list(posterior))

            tool_calls = getattr(msg, "tool_calls", None) or []
            if not tool_calls:
                content = getattr(msg, "content", "") or ""
                emit("agent_end", steps=step, latency_ms=dt, answer_len=len(content))
                session.call_counter = dict(call_counter)
                out = {
                    "status": "answer",
                    "answer": content,
                    "trace": session.trace,
                    "steps": step,
                    "observed_hits": list(session.observed_hits),
                }
                if session.routing is not None:
                    out["routing"] = session.routing
                return out

            # Serialise the assistant echo, then push it BEFORE any pause
            # decision — the OpenAI protocol requires every tool_call in
            # the echo to be followed (eventually) by a matching tool-role
            # message, whether we dispatch it now or on resume.
            assistant_msg: dict = {
                "role": "assistant",
                "content": getattr(msg, "content", None) or None,
                "tool_calls": [
                    {
                        "id": tc.id,
                        "type": getattr(tc, "type", "function"),
                        "function": {
                            "name": tc.function.name,
                            "arguments": tc.function.arguments,
                        },
                    }
                    for tc in tool_calls
                ],
            }
            session.messages.append(assistant_msg)

            # ask_user handling — if any tool_call names ask_user, we pause
            # the loop. Only the FIRST ask_user is honored; a second one in
            # the same turn is treated as a normal (non-existent) tool call
            # and errors via dispatch. Non-ask_user siblings become
            # `pending_deferred` so their tool-role messages get created
            # after the user replies (protocol requirement).
            if self.enable_ask_user:
                ask_tc = next(
                    (tc for tc in tool_calls if tc.function.name == ASK_USER_TOOL_NAME),
                    None,
                )
                if ask_tc is not None:
                    question_text = _extract_ask_question(ask_tc)
                    session.pending_ask = {
                        "tool_call_id": ask_tc.id,
                        "question": question_text,
                    }
                    session.pending_deferred = [
                        {
                            "id": tc.id,
                            "type": getattr(tc, "type", "function"),
                            "function": {
                                "name": tc.function.name,
                                "arguments": tc.function.arguments,
                            },
                        }
                        for tc in tool_calls
                        if tc is not ask_tc
                    ]
                    session.call_counter = dict(call_counter)
                    emit(
                        "agent_ask_user_pause", step=step,
                        question_len=len(question_text or ""),
                        deferred=len(session.pending_deferred),
                    )
                    out = {
                        "status": "question",
                        "question": question_text,
                        "steps": step,
                        "trace": session.trace,
                    }
                    if session.routing is not None:
                        out["routing"] = session.routing
                    return out

            dispatched = self._dispatch_all(
                tool_calls, step=step, trace=session.trace,
                observed=session.observed_hits,
            )
            for tc, result_json in zip(tool_calls, dispatched):
                content = self._maybe_warn_repeat(
                    tc, result_json, call_counter, step=step,
                )
                session.messages.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "content": content,
                })

            # Step-review pause: after every fan-out round the caller may
            # want a human check-in. We build a compact per-call summary
            # from what was just dispatched and return; the loop resumes
            # via resume_from_review(). Skipped when disabled or on the
            # very last permitted step (nothing useful to review then).
            if (
                self.enable_step_review
                and session.step + 1 < session.limit
            ):
                summary = _summarise_fanout(tool_calls, dispatched)
                progress_text = self._maybe_probe_progress(session, step=step)
                # Feed the progress note back into the parent LLM's own
                # message history so the NEXT turn can see it. Without
                # this the summary is human-only and the LLM has no
                # memory of what it just said about its own progress —
                # which was the whole point of asking. The note goes on
                # as a plain assistant message so no tool_call protocol
                # invariant is disturbed.
                if progress_text:
                    session.messages.append({
                        "role": "assistant",
                        "content": f"[进度自评 step={step}] {progress_text}",
                    })
                session.pending_review = {
                    "step": step,
                    "calls": summary,
                }
                if progress_text is not None:
                    session.pending_review["progress"] = progress_text
                session.call_counter = dict(call_counter)
                emit(
                    "agent_step_review_pause",
                    step=step,
                    n_calls=len(summary),
                    has_progress=progress_text is not None,
                )
                out: dict = {
                    "status": "review",
                    "step": step,
                    "calls": summary,
                    "trace": session.trace,
                    "steps": step,
                }
                if progress_text is not None:
                    out["progress"] = progress_text
                if session.routing is not None:
                    out["routing"] = session.routing
                # Advance step counter BEFORE returning — when the caller
                # resumes we're already at step+1's LLM turn, not stuck
                # dispatching the same fan-out again.
                session.step += 1
                return out

            session.step += 1

            # Shadow-answer probe: OpenAI only attaches logprobs to CONTENT
            # tokens, so tool_call-only rounds arrive at the observer with an
            # empty posterior. When collect_logprobs_only or the controller is
            # on AND the primary call returned tool_calls, fire ONE extra LLM
            # call with tool_choice="none" — the model is forced to answer,
            # yielding a real content-token posterior — and REPLACE the empty
            # posterior we appended after the primary call. This shadow probe
            # is diagnostic: its tokens are NOT charged to session.tokens_spent
            # and its output is NOT folded back into session.messages.
            want_posterior = (
                self._confidence_controller is not None
                or self._collect_logprobs_only
            )
            if want_posterior and session.posteriors and not session.posteriors[-1]:
                shadow = self._shadow_posterior_probe(session)
                if shadow:
                    session.posteriors[-1] = list(shadow)

            # Round observer: snapshot hook for offline data collection
            # (T7). Fires unconditionally when set — the controller may
            # or may not be active. Kept upstream of the controller stop
            # branch so observers see the FINAL round on early stops.
            if self.round_observer is not None:
                try:
                    self.round_observer(step, session)
                except Exception:
                    # Never let the observer break the loop; it's a
                    # diagnostic side-channel, not a critical path.
                    pass

            # Confidence controller hook: after the fan-out has completed
            # AND we've decided not to pause for step-review, ask the
            # controller whether the loop should stop early. Only fires
            # when self._confidence_controller is set (config.enabled=True).
            if self._confidence_controller is not None:
                decision = self._controller_step(
                    session, tool_calls, dispatched, step_index=step,
                )
                if decision is not None and decision.stop:
                    # Force one final answer on the evidence gathered so
                    # far, so ``terminated="confidence_stop"`` ships a real
                    # answer downstream evaluators can score instead of a
                    # placeholder. Shared with the max_steps path.
                    synth_answer, synth_dt = self._synthesize_final_answer(session)
                    if not synth_answer:
                        # Model refused / errored → keep the caller informed
                        # rather than silently returning empty text.
                        synth_answer = (
                            "[agent] confidence controller stopped the loop; "
                            "the model produced no final answer."
                        )
                    emit(
                        "agent_end",
                        steps=session.step,
                        terminated="confidence_stop",
                        controller_reason=decision.reason,
                        stab=decision.stab,
                        gain=decision.gain,
                        z=decision.z,
                        synth_latency_ms=synth_dt,
                        answer_len=len(synth_answer),
                    )
                    session.call_counter = dict(call_counter)
                    out = {
                        "status": "answer",
                        "answer": synth_answer,
                        "trace": session.trace,
                        "steps": session.step,
                        "terminated": "confidence_stop",
                        "observed_hits": list(session.observed_hits),
                    }
                    if session.routing is not None:
                        out["routing"] = session.routing
                    return out

        # max_steps reached without a final answer. Force one final
        # answer-synthesis call on the evidence already in the message
        # history so any gold docs already retrieved aren't wasted.
        synth_answer, synth_dt = self._synthesize_final_answer(session)
        if not synth_answer:
            synth_answer = "[agent] max_steps 到达，未生成最终答案"
        emit(
            "agent_end", steps=session.limit, terminated="max_steps",
            synth_latency_ms=synth_dt, answer_len=len(synth_answer),
        )
        session.call_counter = dict(call_counter)
        out = {
            "status": "answer",
            "answer": synth_answer,
            "trace": session.trace,
            "steps": session.limit,
            "terminated": "max_steps",
            "observed_hits": list(session.observed_hits),
        }
        if session.routing is not None:
            out["routing"] = session.routing
        return out

    # ------------------------------------------------------------------
    # forced final-answer synthesis (shared by max_steps + confidence_stop)
    # ------------------------------------------------------------------

    def _synthesize_final_answer(self, session: "AgentSession") -> tuple[str, float]:
        """Force one final answer on the evidence already gathered.

        Both early-stop paths (max_steps, confidence_stop) need to turn the
        accumulated tool results into a real answer instead of a placeholder.
        The naive approach — ``chat_with_tools(messages, tools=[],
        tool_choice="none")`` — was observed to 400 on HepAI/deepseek when
        the message history ends on a ``tool`` role and still carries prior
        ``tool_calls`` while ``tools`` is empty (see 5ade9c9c probe: synth
        returned "" → placeholder → short-llm hallucinated "Ben-Hur").

        Fix: append an explicit user turn asking for the final answer, and
        keep passing the real tool schemas with ``tool_choice="none"`` so
        the request stays well-formed. ``tool_choice="none"`` still forces
        content, not another tool call. The synth user turn is transient —
        NOT folded back into ``session.messages`` — so a resume() after this
        wouldn't see a dangling prompt.

        Returns ``(answer, latency_ms)``. ``answer`` is "" on failure so the
        caller can decide on a fallback string.
        """
        synth_prompt = {
            "role": "user",
            "content": (
                "You have gathered enough evidence. Using ONLY the tool "
                "results above, write the final answer now. Cite the "
                "specific doc_id(s) it rests on. Do not call any more tools."
            ),
        }
        probe_messages = session.messages + [synth_prompt]
        t0 = time.time()
        try:
            msg = self.llm.chat_with_tools(
                probe_messages,
                tools=self._tool_schemas_for_session(),
                tool_choice="none",
            )
            answer = (getattr(msg, "content", None) or "").strip()
        except Exception as e:
            emit("final_synth_error", error=f"{type(e).__name__}: {e}")
            answer = ""
        return answer, round((time.time() - t0) * 1000, 1)

    # ------------------------------------------------------------------
    # confidence controller integration
    # ------------------------------------------------------------------

    def _controller_step(
        self,
        session: AgentSession,
        tool_calls: list,
        dispatched: list,
        *,
        step_index: int,
    ):
        """Build a RoundObservation from the just-completed fan-out and step
        the confidence controller. Returns the controller decision, or None
        when logprobs / observations aren't available (safe → no-op)."""

        if self._confidence_controller is None:
            return None

        # Local import keeps loop.py light when the feature is off.
        from hepuke.agent.confidence.fusion import RoundObservation

        # Doc IDs observed this round vs. the running history. observed_hits
        # is the real-time accumulator populated by _dispatch_all → its rows
        # look like {tool, step, doc_id?, ...}; we bucket by session.step.
        current_ids: list[str] = []
        history_ids: list[str] = []
        for hit in session.observed_hits:
            doc = hit.get("doc_id") if isinstance(hit, dict) else None
            if not doc:
                continue
            if hit.get("step") == step_index:
                current_ids.append(str(doc))
            else:
                history_ids.append(str(doc))

        # All-dry: every sub-agent tool_result was empty. Detected by the
        # subagent contract shape {"answer": "", "dry": True, ...}; keep the
        # check tolerant — non-subagent tool calls simply won't match.
        all_dry = False
        if dispatched:
            dry_flags = []
            for res in dispatched:
                if isinstance(res, dict):
                    dry_flags.append(bool(res.get("dry", False)))
                elif isinstance(res, str):
                    # Best-effort: cheap substring probe. Real parse happens
                    # once the controller integrates the salvage layer.
                    dry_flags.append('"dry": true' in res.lower())
                else:
                    dry_flags.append(False)
            all_dry = bool(dry_flags) and all(dry_flags)

        # mean_cost = per-round search cost (search_calls / round_index).
        # Raw cumulative search_calls grows ~linearly with t so a
        # negative α_cost trained on the collector distribution ends up
        # over/under-penalising online runs whenever the per-round
        # retrieval count differs (e.g. the rewrite-first policy sends
        # fewer parallel searches than the collector's plain-agent runs).
        # Normalising by t gives a bounded ratio (0-few) that is stable
        # across the two policies. grid_search_fusion applies the same
        # transform on the training side.
        search_calls = sum(
            1 for e in session.trace
            if isinstance(e, dict) and e.get("tool") == "search"
        )
        round_index = step_index + 1
        mean_cost = float(search_calls) / max(1, round_index)
        obs = RoundObservation(
            round_index=step_index + 1,  # 1-based per paper convention
            posterior=(
                list(session.posteriors[-1]) if session.posteriors else []
            ),
            sub_agent_embeddings=[],
            current_doc_ids=current_ids,
            history_doc_ids=history_ids,
            mean_cost=mean_cost,
            all_dry=all_dry,
        )
        return self._confidence_controller.step(obs)

    # ------------------------------------------------------------------
    # tool-call dispatch (serial or parallel)
    # ------------------------------------------------------------------

    def _dispatch_all(
        self,
        tool_calls: list,
        *,
        step: int,
        trace: list[dict],
        observed: Optional[list[dict]] = None,
    ) -> list[str]:
        """Run every tool_call in one assistant turn and return their result
        JSON strings in the ORIGINAL order.

        Trace and emit are threadsafe (each writes through a shared lock in
        `observability.trace`). Trace list appends are protected by the
        pool's join — every future completes before we mutate `trace`, so
        no lock is needed on the list itself. The append order is the
        original tool_calls order.

        `observed`, when passed, receives the structured hit records
        `_extract_hits` produces from each tool's raw result. Same
        after-the-join invariant as `trace` — appended in original order,
        no per-append lock needed.
        """
        n = len(tool_calls)
        if n == 0:
            return []
        if not self.parallel_tool_calls or n == 1:
            out: list[str] = []
            for tc in tool_calls:
                out.append(self._dispatch_one(
                    tc, step=step, trace=trace, observed=observed,
                ))
            return out

        results: list[Optional[str]] = [None] * n
        traces: list[Optional[dict]] = [None] * n
        hits_per_call: list[list[dict]] = [[] for _ in range(n)]
        workers = min(n, self.max_parallel)
        # Read the current thread's trace context AT SUBMIT TIME. Worker
        # threads start with an empty thread-local; re-hydrating from this
        # snapshot is the only way to keep agent_id / depth on emits fired
        # from the workers. Snapshot per-call (NOT stored on `self`) so
        # concurrent run()s on the same agent don't clobber each other.
        ctx_snapshot = dict(_trace_ctx())

        def _worker(tc, step):
            if ctx_snapshot:
                set_trace_context(**ctx_snapshot)
            try:
                return self._dispatch_one_capture(tc, step)
            finally:
                clear_trace_context()

        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="tool") as pool:
            futures = {
                pool.submit(_worker, tc, step): idx
                for idx, tc in enumerate(tool_calls)
            }
            for fut in futures:
                idx = futures[fut]
                result_json, trace_entry, hits = fut.result()
                results[idx] = result_json
                traces[idx] = trace_entry
                hits_per_call[idx] = hits
        trace.extend(t for t in traces if t is not None)
        if observed is not None:
            for hs in hits_per_call:
                observed.extend(hs)
        return [r or "" for r in results]

    def _dispatch_one(
        self,
        tc,
        *,
        step: int,
        trace: list[dict],
        observed: Optional[list[dict]] = None,
    ) -> str:
        """Serial-path helper: run one tool_call, append its trace entry."""
        result_json, entry, hits = self._dispatch_one_capture(tc, step)
        trace.append(entry)
        if observed is not None and hits:
            observed.extend(hits)
        return result_json

    def _dispatch_one_capture(self, tc, step: int) -> tuple[str, dict, list[dict]]:
        """Run one tool_call, return (result_json, trace_entry, hits).

        `hits` is the `_extract_hits` output built from the RAW result
        (dict / list, un-serialised), so the salvage path can see the
        full hit list — not the 400-char preview stored on `trace_entry`.
        Kept side-effect-free on shared trace state so the parallel path
        can call it from worker threads and merge afterwards in order.

        Also emits two structured events for the SQLite sink's new tables:
          * `tool_result_full` — one row per tool call, full args + raw
            result (12 KB cap) + latency, keyed by `raw_id` (uuid8). Goes
            to `tool_results_raw`.
          * `retrieval_row` — one row per hit / row-count signal, carrying
            rank / score / snippet + a foreign-key `raw_id`. Goes to
            `retrievals`. Zero rows emitted when the tool returned nothing
            structured (e.g. run_cmd) — the raw row is still written.
        """
        name = tc.function.name
        raw_args = tc.function.arguments or "{}"
        try:
            args = json.loads(raw_args) if isinstance(raw_args, str) else raw_args
        except json.JSONDecodeError:
            args = {"_raw": raw_args}
        t1 = time.time()
        try:
            result = self.tools.dispatch(name, args)
            err = None
        except Exception as e:
            result = {"error": f"{type(e).__name__}: {e}"}
            err = str(e)
        dt = round((time.time() - t1) * 1000, 1)
        result_json = _safe_dump(result)
        hits = _extract_hits(name, args, result, step=step) if err is None else []
        raw_id = uuid.uuid4().hex[:16]
        # Stamp raw_id on every hit so the salvage consumer (and any
        # downstream aggregator) can join back to the full result.
        for h in hits:
            h["raw_id"] = raw_id
        entry = {
            "step": step,
            "tool": name,
            "arguments": args,
            "result_preview": result_json[:400],
            "latency_ms": dt,
            "error": err,
            "raw_id": raw_id,
        }
        try:
            args_json = json.dumps(args, ensure_ascii=False, default=str)
        except (TypeError, ValueError):
            args_json = str(args)
        emit(
            "tool_result_full", step=step, tool=name,
            raw_id=raw_id, args_json=args_json,
            result_json=result_json, error=err, latency_ms=dt,
        )
        for h in hits:
            emit(
                "retrieval_row", step=step, tool=name, raw_id=raw_id,
                doc_id=h.get("doc_id"),
                rank=h.get("rank"),
                score=h.get("score"),
                row_count=h.get("row_count"),
                snippet=h.get("snippet"),
            )
        emit(
            "agent_tool_call", step=step, tool=name,
            latency_ms=dt, error=err,
            args_keys=list(args.keys()) if isinstance(args, dict) else None,
            args_summary=_args_summary(args),
        )
        return result_json, entry, hits

    def _shadow_posterior_probe(
        self, session: AgentSession,
    ) -> list[float]:
        """Fire an off-the-record LLM call to harvest a content-token posterior.

        The primary agent turn may have returned tool_calls only, which means
        the server did NOT attach logprobs (OpenAI only logs CONTENT tokens).
        This shadow call forces the model to answer by passing
        `tool_choice="none"` + an empty `tools=[]`, requests `logprobs=True,
        top_logprobs=top_b`, and returns the length-normalised softmax over
        the first content token's top-B.

        Diagnostic-only: tokens are NOT added to session.tokens_spent, and
        the reply is NOT folded back into session.messages. Any failure
        (client rejects logprobs, network error, empty content) returns
        an empty list — the round then keeps its empty posterior.
        """
        top_b = int(self.confidence_config.top_b) if self.confidence_config is not None else 8
        try:
            msg = self.posterior_llm.chat_with_tools(
                session.messages, tools=[], tool_choice="none",
                logprobs=True, top_logprobs=top_b,
                max_tokens=8,  # one content token is enough; small cap saves cost
            )
        except TypeError:
            return []
        except Exception:
            return []
        posterior = getattr(msg, "_hepuke_posterior", None) or []
        return list(posterior)

    def _maybe_probe_progress(
        self, session: AgentSession, *, step: int,
    ) -> str | None:
        """Ask the parent LLM for a 1-2 sentence progress report.

        Returns the trimmed text, or None when disabled / errored. The
        probe sees the FULL conversation (system + user + assistant
        echoes + tool results) but is forced with `tools=[]` +
        `tool_choice="none"` so it can't spawn another dispatch. Any
        exception is swallowed — a broken probe must not break the loop.
        Emits `agent_perturn_summary` for trace visibility.
        """
        if not self.enable_perturn_summary:
            return None
        probe_messages = session.messages + [
            {"role": "user", "content": PROGRESS_PROBE_PROMPT},
        ]
        t0 = time.time()
        try:
            msg = self.llm.chat_with_tools(
                probe_messages, tools=[], tool_choice="none",
            )
        except Exception as e:  # pragma: no cover — best-effort
            emit(
                "agent_perturn_summary", step=step,
                error=f"{type(e).__name__}: {e}",
            )
            return None
        dt = round((time.time() - t0) * 1000, 1)
        text = (getattr(msg, "content", None) or "").strip()
        emit(
            "agent_perturn_summary", step=step,
            latency_ms=dt, text_len=len(text),
        )
        return text or None

    def _maybe_warn_repeat(
        self, tc, result_json: str, counter: Counter, *, step: int,
    ) -> str:
        """Bump the repeat counter for this (tool, args) and, once it crosses
        `repeated_tool_call_threshold`, prepend a WARNING banner to the tool
        result so the LLM sees on-band evidence that it's re-issuing the
        same call. Returns the (possibly banner-prefixed) content string.
        """
        threshold = self.repeated_tool_call_threshold
        if threshold <= 0:
            return result_json
        key = _canonical_call_key(tc)
        counter[key] += 1
        n = counter[key]
        if n < threshold:
            return result_json
        emit("agent_repeat_warning", step=step, tool=tc.function.name, count=n)
        banner = (
            f"[WARNING: you have called `{tc.function.name}` with the same "
            f"arguments {n} times in this run. The result has not changed. "
            f"STOP re-issuing this call; either try a sibling tool or "
            f"return your best answer with what you already have.]\n"
        )
        return banner + result_json


_DOC_ID_KEYS: tuple[str, ...] = (
    "doc_id", "paper_id", "bam_id", "id",
)


def _pluck_doc_id(row: Any) -> Optional[str]:
    """Return the first plausible doc-id string from a row-shaped dict."""
    if not isinstance(row, dict):
        return None
    for k in _DOC_ID_KEYS:
        v = row.get(k)
        if isinstance(v, str):
            s = v.strip()
            # Skip auto-increment integer-y `id` fields — real doc_ids are
            # tokens like "BAM-00713" or "arXiv:2401.12345". A purely
            # numeric string is almost certainly a row PK, not a doc-id.
            if k == "id" and s.isdigit():
                continue
            if s:
                return s
    return None


def _extract_hits(
    tool_name: str, args: Any, result: Any, *, step: int,
) -> list[dict]:
    """Pull structured hit signals from a tool's RAW result.

    Called after dispatch, before serialisation, so we're not bounded by
    the trace's 400-char `result_preview`. Each entry looks like:
        {tool, step, doc_id?, row_count?, score?, source="observed"}

    Every branch is best-effort: a shape we don't recognise returns [].
    Never raises — a broken extractor cannot take down the tool loop.
    """
    try:
        return _extract_hits_impl(tool_name, args, result, step=step)
    except Exception:
        return []


def _extract_hits_impl(
    tool_name: str, args: Any, result: Any, *, step: int,
) -> list[dict]:
    hits: list[dict] = []

    def _push(**fields: Any) -> None:
        entry = {"tool": tool_name, "step": step, "source": "observed"}
        entry.update({k: v for k, v in fields.items() if v is not None})
        hits.append(entry)

    if tool_name == "sql_query" and isinstance(result, dict):
        row_count = result.get("row_count")
        rows = result.get("rows") or []
        seen: set[str] = set()
        for i, row in enumerate(rows):
            did = _pluck_doc_id(row)
            if did and did not in seen:
                seen.add(did)
                _push(doc_id=did, row_count=row_count,
                      rank=i, snippet=_row_snippet(row))
        if not seen and isinstance(row_count, int) and row_count > 0:
            _push(row_count=row_count)
        return hits

    if tool_name == "filter_rows" and isinstance(result, list):
        seen = set()
        for i, row in enumerate(result):
            did = _pluck_doc_id(row)
            if did and did not in seen:
                seen.add(did)
                _push(doc_id=did, rank=i, snippet=_row_snippet(row))
        if not seen and result and isinstance(result[0], dict) \
                and "error" not in result[0]:
            _push(row_count=len(result))
        return hits

    if tool_name == "retrieve_vec" and isinstance(result, list):
        seen = set()
        rank_i = 0
        for hit in result:
            did = _pluck_doc_id(hit)
            if did and did not in seen:
                seen.add(did)
                score = hit.get("score") if isinstance(hit, dict) else None
                _push(doc_id=did, score=score, rank=rank_i,
                      snippet=_row_snippet(hit))
                rank_i += 1
        return hits

    if tool_name == "search" and isinstance(result, list):
        # Main-agent hybrid search tool (registered in agent/defaults.py as
        # `search`). Result is `_compact`-ed hits with score + doc_id + content.
        # Without this branch observed_hits stays empty even when retrieval
        # fires — breaks doc_ids/novelty features for T7 collector runs.
        seen = set()
        rank_i = 0
        for hit in result:
            did = _pluck_doc_id(hit)
            if did and did not in seen:
                seen.add(did)
                score = hit.get("score") if isinstance(hit, dict) else None
                _push(doc_id=did, score=score, rank=rank_i,
                      snippet=_row_snippet(hit))
                rank_i += 1
        return hits

    if tool_name == "search_index" and isinstance(result, dict):
        seen = set()
        rank_i = 0
        for row in result.get("hits") or []:
            did = _pluck_doc_id(row)
            if did and did not in seen:
                seen.add(did)
                _push(doc_id=did, rank=rank_i, snippet=_row_snippet(row))
                rank_i += 1
        return hits

    if tool_name == "grep" and isinstance(result, list):
        # grep returns [{file, line, text, is_match}]. Extract BAM-style
        # doc_ids from the filename when the row is a match.
        import re as _re
        _FILE_ID = _re.compile(r"(BAM-\d{3,5})", _re.IGNORECASE)
        seen = set()
        rank_i = 0
        for row in result:
            if not isinstance(row, dict) or not row.get("is_match"):
                continue
            fname = row.get("file") or ""
            m = _FILE_ID.search(fname)
            if m:
                did = m.group(1).upper()
                if did not in seen:
                    seen.add(did)
                    text = row.get("text") or ""
                    snippet = str(text).strip().replace("\n", " ")[:300]
                    _push(doc_id=did, rank=rank_i, snippet=snippet or None)
                    rank_i += 1
        return hits

    if tool_name in {"read_summary", "read_memo", "preview_memo"}:
        # args carry the doc_id; the result confirms it was found (has
        # text / body). If we can't verify success cheaply just trust
        # the arg — the tool would have raised otherwise.
        did = None
        if isinstance(args, dict):
            v = args.get("doc_id")
            if isinstance(v, str) and v.strip():
                did = v.strip()
        if did:
            ok = True
            if isinstance(result, dict) and result.get("error"):
                ok = False
            if ok:
                snippet = _row_snippet(result) if isinstance(result, dict) else None
                _push(doc_id=did, rank=0, snippet=snippet)
        return hits

    if tool_name == "get_full_doc" and isinstance(result, dict):
        did = result.get("doc_id")
        if isinstance(did, str) and did.strip() and not result.get("error"):
            _push(doc_id=did.strip(), rank=0, snippet=_row_snippet(result))
        return hits

    return hits


def _row_snippet(row: Any, max_len: int = 300) -> str | None:
    """Compact one-line preview for a hit row. Prefers title / abstract_lede
    / text / content over a full JSON dump. Returns None when the row has
    nothing usable, so the caller can drop the column.
    """
    if not isinstance(row, dict):
        s = str(row).strip().replace("\n", " ")
        return (s[:max_len] or None)
    for key in ("title", "abstract_lede", "abstract", "text", "content", "body", "snippet"):
        v = row.get(key)
        if isinstance(v, str) and v.strip():
            s = v.strip().replace("\n", " ")
            return s[:max_len]
    try:
        s = json.dumps(row, ensure_ascii=False, default=str)
    except (TypeError, ValueError):
        s = str(row)
    s = s.replace("\n", " ")
    return s[:max_len] or None


def _canonical_call_key(tc) -> str:
    """Stable hash-friendly key for one tool_call. Args are canonicalised
    via json.dumps(sort_keys=True) so semantically-identical arg dicts in
    different key order collapse to the same key.
    """
    name = tc.function.name
    raw = tc.function.arguments or "{}"
    try:
        args = json.loads(raw) if isinstance(raw, str) else raw
    except json.JSONDecodeError:
        args = {"_raw": raw}
    try:
        canon = json.dumps(args, ensure_ascii=False, sort_keys=True, default=str)
    except (TypeError, ValueError):
        canon = str(args)
    return f"{name}::{canon}"


def _safe_dump(obj: Any, max_len: int = 12_000) -> str:
    try:
        s = json.dumps(obj, ensure_ascii=False, default=str)
    except (TypeError, ValueError):
        s = json.dumps(str(obj), ensure_ascii=False)
    if len(s) > max_len:
        s = s[:max_len] + " ... [truncated]"
    return s


def _args_summary(args: Any, max_len: int = 120) -> str:
    """Compact one-line preview of a tool_call's args for console logs.

    Prefers the most informative single field (question / pattern /
    filter_expr / sql / query / doc_id) when present; falls back to a
    truncated JSON dump. Never raises — sink output must not break the run.
    """
    if not isinstance(args, dict):
        return ""
    for key in ("question", "pattern", "filter_expr", "sql", "query", "doc_id", "expr"):
        v = args.get(key)
        if isinstance(v, str) and v.strip():
            body = v.strip().replace("\n", " ")
            if len(body) > max_len:
                body = body[:max_len] + "…"
            return f"{key}={body}"
    try:
        s = json.dumps(args, ensure_ascii=False, default=str)
    except (TypeError, ValueError):
        s = str(args)
    return s if len(s) <= max_len else s[:max_len] + "…"


def _extract_ask_question(tc) -> str:
    """Pull the `question` string out of an ask_user tool_call. Robust to
    malformed arguments — a garbled JSON call still gets a placeholder
    question rather than raising, because the loop is already committed
    to pausing.
    """
    raw = tc.function.arguments or "{}"
    try:
        args = json.loads(raw) if isinstance(raw, str) else raw
    except json.JSONDecodeError:
        return "(the assistant asked something but its arguments were malformed)"
    if isinstance(args, dict):
        q = args.get("question")
        if isinstance(q, str) and q.strip():
            return q.strip()
    return "(the assistant asked something but did not include a question)"


class _DictToolCall:
    """Adapter that lets a dict-shaped deferred tool_call ride the same
    dispatch paths as the OpenAI SDK's ChatCompletionMessageToolCall.
    """
    def __init__(self, d: dict) -> None:
        self.id = d["id"]
        self.type = d.get("type", "function")
        fn = d.get("function") or {}
        self.function = type("_F", (), {
            "name": fn.get("name", ""),
            "arguments": fn.get("arguments", "{}"),
        })()


def _summarise_fanout(tool_calls: list, dispatched: list[str]) -> list[dict]:
    """Build a per-call summary for a step-review pause.

    For each tool_call in the fan-out, emit:
      {tool, args_summary, dry, evidence_n, answer_preview}
    Pulls `dry` / `evidence` fields out of subagent contract JSON when
    present (that's how the milvus_/memo_/hypernews_db subagents shape
    their return); for plain tools it falls back to a short preview of
    the raw JSON so the caller still sees SOMETHING useful.
    """
    out: list[dict] = []
    for tc, result_json in zip(tool_calls, dispatched):
        raw_args = getattr(tc.function, "arguments", "") or "{}"
        try:
            args = json.loads(raw_args) if isinstance(raw_args, str) else raw_args
        except json.JSONDecodeError:
            args = {}
        entry: dict = {
            "tool": tc.function.name,
            "args_summary": _args_summary(args),
        }
        try:
            parsed = json.loads(result_json) if isinstance(result_json, str) else result_json
        except (json.JSONDecodeError, TypeError):
            parsed = None
        if isinstance(parsed, dict):
            if "dry" in parsed:
                entry["dry"] = bool(parsed.get("dry"))
            evi = parsed.get("evidence")
            if isinstance(evi, list):
                entry["evidence_n"] = len(evi)
            ans = parsed.get("answer")
            if isinstance(ans, str):
                preview = ans.strip().replace("\n", " ")
                entry["answer_preview"] = preview[:200] + ("…" if len(preview) > 200 else "")
        if "answer_preview" not in entry:
            preview = (result_json or "").replace("\n", " ")
            entry["answer_preview"] = preview[:200] + ("…" if len(preview) > 200 else "")
        out.append(entry)
    return out
