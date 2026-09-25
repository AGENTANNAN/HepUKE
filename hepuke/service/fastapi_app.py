"""FastAPI application exposing HepUKE as a REST service.

Run:
    hepuke serve                     # via CLI, most convenient
    uvicorn hepuke.service.fastapi_app:get_app --factory   # explicit

Auth: `Authorization: Bearer <cfg.service.api_key>` when the key is non-empty.
"""
from __future__ import annotations

import json as _json
import logging
import os
import time as _time
import uuid as _uuid
from typing import Any, Callable, List, Optional

from fastapi import Depends, FastAPI, HTTPException, Header, Path
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from hepuke.agent.defaults import build_default_tools
from hepuke.agent.loop import RagAgent
from hepuke.agent.memory import SessionMemory
from hepuke.agent.prompts import build_system_prompt
from hepuke.agent.tools import ToolRegistry
from hepuke.config import Config, load_config
from hepuke.core.documents import dict_to_document
from hepuke.core.index import HepUKE
from hepuke.service import schemas
from hepuke.service.sessions import SessionStore
from hepuke.tools import memo_tools

logger = logging.getLogger(__name__)

# pymilvus writes RPC errors (and full tracebacks) via its own logger even
# when the caller catches the exception. Our tools ALREADY catch and turn
# these into structured error dicts for the LLM, so surfacing them in
# stderr twice is noise. Silence pymilvus at WARNING+; genuine errors are
# still visible in the trace sink.
logging.getLogger("pymilvus").setLevel(logging.ERROR + 10)


class UTF8JSONResponse(JSONResponse):
    """JSONResponse that emits real UTF-8 instead of `\\uXXXX` escapes.

    FastAPI's default `JSONResponse` calls `json.dumps(..., ensure_ascii=True)`,
    which makes every CJK character render as an escape sequence to any
    downstream `json.tool` / logs / curl. Overriding `render()` with
    `ensure_ascii=False` gives the same wire format (still valid JSON) but
    keeps the bytes readable. All routes inherit this by setting
    `default_response_class=UTF8JSONResponse` on the `FastAPI(...)` ctor.
    """

    def render(self, content: Any) -> bytes:
        return _json.dumps(
            content,
            ensure_ascii=False,
            allow_nan=False,
            indent=None,
            separators=(",", ":"),
        ).encode("utf-8")


def _bearer_dependency(cfg: Config) -> Callable:
    async def verify(authorization: str | None = Header(default=None)) -> None:
        expected = cfg.service.api_key
        if not expected:
            return
        if not authorization or not authorization.lower().startswith("bearer "):
            raise HTTPException(status_code=401, detail="Missing bearer token")
        token = authorization.split(None, 1)[1].strip()
        if token != expected:
            raise HTTPException(status_code=403, detail="Invalid bearer token")

    return verify


# ---------------------------------------------------------------------------
# Additional schemas
# ---------------------------------------------------------------------------

class GrepRequest(BaseModel):
    pattern: str
    doc_ids: Optional[List[str]] = None
    context: int = 3
    ignore_case: bool = False
    max_results: int = 200
    path_root: Optional[str] = None


class AskRequest(BaseModel):
    question: str
    max_steps: Optional[int] = None
    collection: Optional[str] = None


class ChatRequest(BaseModel):
    """One turn on the /agent/chat conversation.

    Modes:
      * fresh turn: `question` set, `session_id` / `user_reply` unset.
      * resume: `session_id` set, `user_reply` set to the human's answer
        to the previous ask_user question. `question` is ignored.
    """
    question: Optional[str] = None
    session_id: Optional[str] = None
    user_reply: Optional[str] = None
    max_steps: Optional[int] = None
    collection: Optional[str] = None


class UpdateMemoRequest(BaseModel):
    doc_id: str
    text: str
    metadata: dict = {}
    split_chunk: bool = True


def _shape_chat_response(result: dict, session_id: str | None, store: "SessionStore") -> dict:
    """Turn a RagAgent run/resume result into the /agent/chat wire shape.

    Client contract:
      status='answer'   → {session_id: null, status, answer, steps, ...}
      status='question' → {session_id: <sid>, status, question, steps, ...}
      status='review'   → {session_id: <sid>, status, calls, step, ...}
    Any 'routing' / 'trace' / 'terminated' from the loop pass through.
    """
    status = result.get("status", "answer")
    out: dict = {"status": status, "session_id": session_id}
    if status == "question":
        out["question"] = result.get("question", "")
    elif status == "review":
        # Fan-out summary — one entry per tool_call in the just-finished
        # round. Client should render it and prompt the user for
        # continue/stop/redirect.
        out["calls"] = result.get("calls", [])
        out["step"] = result.get("step")
        # Optional 1-2 sentence progress summary from the parent LLM,
        # produced right before this pause when `enable_perturn_summary`
        # is on. Absent when the feature is off or the probe errored.
        if result.get("progress") is not None:
            out["progress"] = result["progress"]
    else:
        out["answer"] = result.get("answer", "")
    for k in ("trace", "steps", "routing", "terminated", "tokens_spent"):
        if k in result:
            out[k] = result[k]
    # Clean up finished sessions so they don't sit until TTL reaps them.
    if status == "answer" and session_id is not None:
        store.drop(session_id)
    return out


def create_app(config_path: str | None = None) -> FastAPI:
    cfg = load_config(config_path or os.environ.get("HEPUKE_CONFIG"))
    rag = HepUKE(cfg)
    prefix = cfg.service.route_prefix.rstrip("/")
    app = FastAPI(
        title="HepUKE",
        version="0.3.0",
        default_response_class=UTF8JSONResponse,
    )
    auth = Depends(_bearer_dependency(cfg))

    memo_root = cfg.wiki.memo_md_root

    # Interactive-chat session store (paused RagAgents awaiting ask_user
    # replies). In-memory, process-local; disappears on restart. Sized for
    # a handful of concurrent conversations — not a general session cache.
    session_store = SessionStore(max_sessions=256, ttl_seconds=3600.0)

    def _resolve_token_budget() -> tuple[int, bool]:
        """Env → config precedence for the token-budget backstop.
        Shared by /agent/ask and /agent/chat.
        """
        env_enabled = os.environ.get("HEPUKE_TOKEN_BUDGET_ENABLED")
        if env_enabled is not None:
            enabled = env_enabled.strip().lower() not in {
                "", "0", "false", "no", "off",
            }
        else:
            enabled = cfg.agent.token_budget_enabled
        env_budget = os.environ.get("HEPUKE_TOKEN_BUDGET")
        if env_budget:
            try:
                budget = max(0, int(env_budget))
            except ValueError:
                budget = cfg.agent.token_budget
        else:
            budget = cfg.agent.token_budget
        return budget, enabled

    def _build_registry(collection: str, memory: Optional[SessionMemory] = None) -> ToolRegistry:
        registry = ToolRegistry()
        registry.register_all(
            build_default_tools(
                rag.retriever,
                memo_root=memo_root,
                default_collection=collection,
                schema_profile=rag.schema_profile,
                use_subagents=cfg.agent.use_subagents,
                llm=rag.llm,
                hypernews_db_path=cfg.agent.hypernews_db_path,
                memo_index_path=cfg.wiki.memo_index_path,
                memo_summaries_root=cfg.wiki.memo_summaries_root,
                rewriter=getattr(rag, "rewriter", None),
                memory=memory,
                default_with_rerank=cfg.retrieval.with_reranker,
            )
        )
        return registry

    def _make_memory() -> Optional[SessionMemory]:
        """Fresh per-request `SessionMemory` under `HEPUKE_AGENT_MEMORY_ROOT`,
        or `None` when the env var is unset (memory disabled, matches the
        pre-wiring behaviour). One uuid-named subdirectory per /agent/chat
        fresh turn; a paused session keeps the same directory across
        resumes because the `RagAgent` instance is stored intact in
        `SessionStore`. Direct/decompose branches also get a directory —
        it becomes garbage on process exit, which is cheap enough not to
        justify a delayed-creation refactor of the router path."""
        root = os.environ.get("HEPUKE_AGENT_MEMORY_ROOT")
        if not root:
            return None
        session_dir = os.path.join(root, _uuid.uuid4().hex[:12])
        return SessionMemory(session_dir)

    # ---- startup warm-up ------------------------------------------------
    # deepseek/HepAI prefix caching operates at 2 KiB block granularity —
    # a cold prefix costs ~2× latency + full prompt tokens. Fire one tiny
    # request with the exact router SYSTEM prompt at boot so the cache
    # entry is populated before the first real user request. Failure is
    # non-fatal — the first request will just pay the miss.
    if cfg.agent.enabled and cfg.agent.route:
        try:
            from hepuke.agent.router import ROUTER_SYSTEM_PROMPT
            t0 = _time.time()
            rag.llm.chat_with_tools(
                [
                    {"role": "system", "content": ROUTER_SYSTEM_PROMPT},
                    {"role": "user", "content": "warmup"},
                ],
                tools=[], tool_choice="none",
                temperature=0, max_tokens=8,
            )
            logger.info(
                "router prompt warmed in %.0fms (cache ready for first request)",
                (_time.time() - t0) * 1000,
            )
        except Exception as e:  # pragma: no cover — depends on live LLM
            logger.warning("router warm-up skipped: %s", e)

    def _ensure_collection(name: str) -> None:
        if rag.collection_name != name:
            rag.connect_collection(name)

    @app.get(f"{prefix}/health", tags=["meta"])
    def health() -> dict:
        return {"ok": True, "collection": rag.collection_name}

    @app.get(f"{prefix}/agent/capabilities", tags=["agent"])
    def agent_capabilities() -> dict:
        """List the knowledge bases actually mounted right now.

        Reflects the real `ToolRegistry` (not a hand-rolled description),
        so if a subagent gets added/dropped this endpoint tracks it.
        Also surfaces the on-disk paths / db files so the user can tell
        which corpus each subagent owns without reading source.
        """
        registry = _build_registry(cfg.milvus.default_collection)
        tools_info: list[dict] = []
        for tool_name in registry.names():
            tool = registry.get(tool_name)
            tools_info.append({
                "name": tool.name,
                "description": tool.description,
            })
        return {
            "collection": rag.collection_name,
            "use_subagents": cfg.agent.use_subagents,
            "tools": tools_info,
            "mounts": {
                "milvus_collection": rag.collection_name,
                "memo_md_root": cfg.wiki.memo_md_root,
                "memo_index_path": cfg.wiki.memo_index_path,
                "memo_summaries_root": cfg.wiki.memo_summaries_root,
                "hypernews_db_path": cfg.agent.hypernews_db_path,
            },
        }

    # ---- collections ---------------------------------------------------
    @app.get(f"{prefix}/collections", dependencies=[auth], tags=["collections"])
    def list_collections() -> dict:
        return {"collections": rag.get_collections()}

    @app.post(f"{prefix}/collections/{{name}}", dependencies=[auth], tags=["collections"])
    def create_collection(name: str = Path(...)) -> dict:
        return {"collections": rag.create_collection(name)}

    @app.delete(f"{prefix}/collections/{{name}}", dependencies=[auth], tags=["collections"])
    def drop_collection(name: str = Path(...)) -> dict:
        return {"collections": rag.drop_collection(name)}

    # ---- documents -----------------------------------------------------
    @app.post(
        f"{prefix}/{{collection}}/insert",
        dependencies=[auth],
        response_model=schemas.InsertResponse,
        tags=["docs"],
    )
    def insert(collection: str, req: schemas.InsertRequest) -> schemas.InsertResponse:
        _ensure_collection(collection)
        try:
            document = dict_to_document(req.doc.model_dump(exclude_none=True))
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
        msg, doc_id = rag.insert(document, split_chunk=req.split_chunk)
        return schemas.InsertResponse(msg=msg, doc_id=doc_id)

    @app.delete(
        f"{prefix}/{{collection}}/docs/{{doc_id}}",
        dependencies=[auth],
        tags=["docs"],
    )
    def delete(collection: str, doc_id: str) -> schemas.SimpleResponse:
        _ensure_collection(collection)
        result = rag.delete(doc_id)
        return schemas.SimpleResponse(ok=result is True, detail=result)

    @app.post(
        f"{prefix}/{{collection}}/docs",
        dependencies=[auth],
        tags=["docs"],
    )
    def lookup_docs(collection: str, req: schemas.DocLookupRequest) -> dict:
        _ensure_collection(collection)
        return {"docs": rag.get_doc(req.key, req.value)}

    # ---- retrieval / query --------------------------------------------
    @app.post(
        f"{prefix}/{{collection}}/retrieve",
        dependencies=[auth],
        response_model=schemas.RetrieveResponse,
        tags=["query"],
    )
    def retrieve(collection: str, req: schemas.RetrieveRequest) -> schemas.RetrieveResponse:
        _ensure_collection(collection)
        hits = rag.retrieve(
            req.content,
            similarity_top_k=req.similarity_top_k,
            filters=req.filters,
            condition=req.condition,
            with_reranker=req.with_reranker,
        )
        return schemas.RetrieveResponse(
            nodes=[
                schemas.RetrievedNode(
                    score=h.get("score"),
                    text=h.get("content", ""),
                    metadata={k: h.get(k) for k in
                              ("doc_id", "section_type", "section_title", "title")
                              if h.get(k) is not None},
                )
                for h in hits
            ]
        )

    @app.post(
        f"{prefix}/{{collection}}/query",
        dependencies=[auth],
        response_model=schemas.QueryResponse,
        tags=["query"],
    )
    def query(collection: str, req: schemas.QueryRequest) -> schemas.QueryResponse:
        _ensure_collection(collection)
        out = rag.query(
            req.content,
            similarity_top_k=req.similarity_top_k,
            filters=req.filters,
            condition=req.condition,
            with_reranker=req.with_reranker,
        )
        return schemas.QueryResponse(
            answer=out["answer"],
            sources=[schemas.RetrievedNode(**s) for s in out["sources"]],
        )

    # ---- memo tools ---------------------------------------------------
    @app.post(f"{prefix}/memo/grep", dependencies=[auth], tags=["memo"])
    def memo_grep(req: GrepRequest) -> dict:
        try:
            hits = memo_tools.grep(
                pattern=req.pattern,
                path_root=req.path_root or memo_root,
                doc_ids=req.doc_ids,
                context=req.context,
                ignore_case=req.ignore_case,
                max_results=req.max_results,
            )
        except memo_tools.MemoToolError as e:
            raise HTTPException(status_code=400, detail=str(e))
        return {"hits": hits}

    @app.get(f"{prefix}/memo/list", dependencies=[auth], tags=["memo"])
    def memo_list(pattern: str = "*.md", path_root: Optional[str] = None) -> dict:
        try:
            docs = memo_tools.list_memos(path_root=path_root or memo_root, pattern=pattern)
        except memo_tools.MemoToolError as e:
            raise HTTPException(status_code=400, detail=str(e))
        return {"docs": docs}

    @app.get(f"{prefix}/memo/{{doc_id}}", dependencies=[auth], tags=["memo"])
    def memo_read(doc_id: str, collection: Optional[str] = None) -> dict:
        if collection:
            _ensure_collection(collection)
        return memo_tools.read_memo(
            doc_id=doc_id,
            retriever=rag.retriever if collection else None,
            collection=collection,
            path_root=memo_root,
        )

    @app.get(
        f"{prefix}/memo/{{doc_id}}/{{section}}",
        dependencies=[auth],
        tags=["memo"],
    )
    def memo_section(doc_id: str, section: str, collection: str) -> dict:
        _ensure_collection(collection)
        return memo_tools.get_section(
            doc_id=doc_id, section=section,
            retriever=rag.retriever, collection=collection,
        )

    @app.post(
        f"{prefix}/{{collection}}/memo/update",
        dependencies=[auth],
        tags=["memo"],
    )
    def memo_update(collection: str, req: UpdateMemoRequest) -> dict:
        _ensure_collection(collection)
        from hepuke.document import Document
        d = Document(text=req.text, metadata=req.metadata, doc_id=req.doc_id)
        # Delete old copy first for idempotency.
        rag.delete(req.doc_id)
        msg, doc_id = rag.insert(d, split_chunk=req.split_chunk)
        return {"msg": msg, "doc_id": doc_id}

    @app.post(f"{prefix}/{{collection}}/memo/persist", dependencies=[auth], tags=["memo"])
    def memo_persist(collection: str) -> dict:
        _ensure_collection(collection)
        rag.persist()
        return {"ok": True}

    # ---- agent --------------------------------------------------------
    @app.post(f"{prefix}/agent/ask", dependencies=[auth], tags=["agent"])
    def agent_ask(req: AskRequest) -> dict:
        collection = req.collection or cfg.milvus.default_collection
        _ensure_collection(collection)
        memory = _make_memory()
        if memory is not None:
            try:
                memory.write_user_intent(req.question)
            except Exception:
                pass
        registry = _build_registry(collection, memory=memory)
        token_budget, token_budget_enabled = _resolve_token_budget()
        agent = RagAgent(
            llm=rag.llm,
            tools=registry,
            max_steps=req.max_steps or cfg.agent.max_steps,
            system_prompt=build_system_prompt(
                rag.schema_profile, use_subagents=cfg.agent.use_subagents,
            ),
            parallel_tool_calls=cfg.agent.parallel_tool_calls,
            max_parallel=cfg.agent.max_parallel,
            route=cfg.agent.route,
            decompose_max_parallel=cfg.agent.decompose_max_parallel,
            repeated_tool_call_threshold=cfg.agent.repeated_tool_call_threshold,
            token_budget=token_budget,
            token_budget_enabled=token_budget_enabled,
            memory=memory,
        )
        return agent.run(req.question)

    @app.post(f"{prefix}/agent/chat", dependencies=[auth], tags=["agent"])
    def agent_chat(req: ChatRequest) -> dict:
        """Interactive chat with ask_user round-trips.

        Two turn shapes:
          * FRESH — client sends `question`, gets either a final answer
            (status='answer') or a clarifying question (status='question' +
            session_id).
          * RESUME — client sends `session_id` + `user_reply` to answer a
            prior ask_user; response follows the same shape.

        Session state (paused RagAgent + conversation) is process-local; a
        missing session_id returns 404 and the client must restart.
        """
        # Resume path takes precedence when session_id is provided.
        if req.session_id:
            entry = session_store.get(req.session_id)
            if entry is None:
                raise HTTPException(status_code=404, detail="unknown or expired session")
            if req.user_reply is None:
                raise HTTPException(
                    status_code=400,
                    detail="resume requires user_reply",
                )
            # Pick the right resume path from session state — the same
            # (session_id, user_reply) shape covers both ask_user replies
            # and step-review continue/stop/redirect. This keeps the wire
            # protocol tiny; the state on the server disambiguates.
            try:
                if entry.session.pending_ask:
                    result = entry.agent.resume(entry.session, req.user_reply)
                elif entry.session.pending_review:
                    result = entry.agent.resume_from_review(
                        entry.session, req.user_reply,
                    )
                else:
                    raise HTTPException(
                        status_code=400,
                        detail="session has no pending pause to resume",
                    )
            except ValueError as e:
                raise HTTPException(status_code=400, detail=str(e))
            # If the resume yielded ANOTHER pause (e.g. the next fan-out
            # round also triggers step-review), keep the same session_id
            # so the client can keep round-tripping without new state.
            if result.get("status") in ("question", "review"):
                return _shape_chat_response(result, req.session_id, session_store)
            return _shape_chat_response(result, req.session_id, session_store)

        # Fresh turn.
        if not req.question:
            raise HTTPException(status_code=400, detail="question required for a fresh chat turn")
        collection = req.collection or cfg.milvus.default_collection
        _ensure_collection(collection)
        memory = _make_memory()
        # Write intent at the HTTP entry point so all three router branches
        # (direct / decompose / rag) get intent persisted. `_new_session`
        # inside RagAgent covers the rag branch only; direct/decompose
        # bypass it and used to leave `memory/` empty.
        if memory is not None:
            try:
                memory.write_user_intent(req.question)
            except Exception:
                pass
        registry = _build_registry(collection, memory=memory)
        token_budget, token_budget_enabled = _resolve_token_budget()
        agent = RagAgent(
            llm=rag.llm,
            tools=registry,
            max_steps=req.max_steps or cfg.agent.max_steps,
            system_prompt=build_system_prompt(
                rag.schema_profile, use_subagents=cfg.agent.use_subagents,
            ),
            parallel_tool_calls=cfg.agent.parallel_tool_calls,
            max_parallel=cfg.agent.max_parallel,
            route=cfg.agent.route,
            decompose_max_parallel=cfg.agent.decompose_max_parallel,
            repeated_tool_call_threshold=cfg.agent.repeated_tool_call_threshold,
            token_budget=token_budget,
            token_budget_enabled=token_budget_enabled,
            enable_ask_user=True,
            enable_step_review=True,
            enable_perturn_summary=cfg.agent.enable_perturn_summary,
            memory=memory,
        )

        # Run router first so decompose and HEP-entity override fire on the
        # interactive path too. Direct / decompose complete synchronously
        # (no interactive pause), so we surface the answer as-is. Only the
        # single-topic `rag` branch owns a session that may pause via
        # ask_user or step_review — for THAT branch we manage the loop
        # ourselves so we can hand a live AgentSession to the store.
        if agent.route:
            from hepuke.agent.router import classify as _classify
            try:
                decision = _classify(req.question, llm=agent.llm)
            except Exception:  # pragma: no cover — router must never break
                decision = None
        else:
            decision = None

        if decision is not None and decision.mode == "direct":
            result = agent._run_direct(req.question, decision=decision)
            result.setdefault("status", "answer")
            return _shape_chat_response(result, None, session_store)
        if decision is not None and decision.mode == "decompose":
            result = agent._run_decompose(
                req.question, decision=decision, limit=agent.max_steps,
            )
            result.setdefault("status", "answer")
            return _shape_chat_response(result, None, session_store)

        # rag branch (or no router). Own the session so pause can persist.
        session = agent._new_session(
            req.question, limit=agent.max_steps, routing=decision,
        )
        # Mint the store session_id BEFORE _advance runs so every emit /
        # retrieval_row / tool_result_full row is stamped with the same id
        # that the store will later use as its key. Without this the first
        # turn's rows would carry the AgentSession's fallback (agent_id)
        # and later turns would carry the store sid — impossible to join.
        import uuid as _uuid
        session.session_id = _uuid.uuid4().hex
        result = agent._advance(session)
        if result.get("status") in ("question", "review"):
            sid = session_store.create(agent, session)
            return _shape_chat_response(result, sid, session_store)
        return _shape_chat_response(result, None, session_store)

    return app


def get_app() -> FastAPI:
    """Factory entrypoint for `uvicorn --factory`."""
    return create_app()


if __name__ == "__main__":  # pragma: no cover
    import uvicorn
    cfg = load_config()
    uvicorn.run(get_app(), host=cfg.service.host, port=cfg.service.port, factory=False)
