"""Default tool bindings: wraps HybridRetriever + memo_tools into the registry."""
from __future__ import annotations

from typing import Any, Optional

from hepuke.agent.tools import Tool, ToolRegistry
from hepuke.tools import memo_tools
from hepuke.tools import shell_tools


_VALID_COLLECTION_RE = __import__("re").compile(r"^[A-Za-z0-9_]+$")
# Forbid vector-field selection in filter_rows output (huge, useless to the LLM).
_VECTOR_FIELD_SUFFIXES = ("_vector", "_vec", "embedding", "embeddings")


def _search_tool(
    retriever,
    default_collection: Optional[str],
    *,
    default_with_rerank: bool = True,
    return_top: Optional[int] = None,
) -> Tool:
    def _call(
        question: str,
        collection: Optional[str] = None,
        top_k: int = 5,
        section_types: Optional[list] = None,
        doc_id: Optional[str] = None,
        with_rerank: Optional[bool] = None,
    ) -> list:
        # LLMs sometimes stuff a document id (e.g. "BAM-00003_...") into
        # `collection`. Milvus only accepts `[A-Za-z0-9_]+`. If the caller
        # passed something that can't be a collection name, treat it as a
        # doc_id filter and fall back to the default collection.
        if collection and not _VALID_COLLECTION_RE.match(collection):
            if not doc_id:
                doc_id = collection
            collection = None
        coll = collection or default_collection
        if not coll:
            return [{"error": "collection is required and no default configured"}]
        # If the LLM didn't specify with_rerank explicitly, use the default
        # baked in at tool-build time from `config.retrieval.with_reranker`.
        use_rerank = default_with_rerank if with_rerank is None else with_rerank

        import logging as _logging
        _logging.getLogger(__name__).info(
            "[search-tool-debug] q=%r top_k=%r return_top=%r use_rerank=%r "
            "doc_id=%r collection=%r section_types=%r",
            (question[:80] + "…") if isinstance(question, str) and len(question) > 80 else question,
            top_k, return_top, use_rerank, doc_id, coll, section_types,
        )

        # Pure single-query search: no rewrite, no fan-out. Query rewriting
        # is now an explicit `rewrite` tool the agent calls once at round 1
        # so the parent LLM sees N candidate queries as separate messages
        # and picks which to actually search. Baking a hidden fan-out into
        # every `search` collapsed those N angles into a single RRF-merged
        # top-k that the agent could not reason about — see 5ade9c9c probe.
        hits = retriever.search(
            question=question,
            collection=coll,
            top_k=top_k,
            section_types=section_types,
            doc_id=doc_id,
            with_rerank=use_rerank,
        )
        _logging.getLogger(__name__).info(
            "[search-tool-debug/post] q=%r n_hits=%d top3_doc_ids=%s",
            (question[:60] + "…") if isinstance(question, str) and len(question) > 60 else question,
            len(hits),
            [h.get("doc_id") for h in hits[:3]],
        )
        # `return_top` narrows what the AGENT sees without touching recall:
        # the retriever still fetched `top_k*4` and reranked the full pool,
        # so `hits[0]` is the highest cross-encoder match. Returning only
        # the top result per search makes each call "narrow and deep" — the
        # agent reads one strongly-relevant doc's full content instead of
        # skimming N shallow hits. Round-2 fan-out still gives breadth (N
        # lanes = N top-1 docs). None → return the full top_k (default).
        if return_top is not None:
            hits = hits[:return_top]
        return [_compact(h) for h in hits]

    parameters = {
        "type": "object",
        "properties": {
            "question": {"type": "string", "description": "the query text"},
            "collection": {
                "type": "string",
                "description": (
                    "Milvus collection name — usually leave unset to use the "
                    "default. Do NOT put a document id here; use `doc_id` for "
                    "that. Only alphanumerics and underscores allowed."
                ),
            },
            "top_k": {"type": "integer", "default": 5, "minimum": 1, "maximum": 20},
            "section_types": {"type": "array", "items": {"type": "string"},
                              "description": "whitelist of section_type values"},
            "doc_id": {"type": "string",
                       "description": "restrict search to one document id (e.g. 'BAM-00003_...')"},
            "with_rerank": {"type": "boolean", "default": False},
        },
        "required": ["question"],
    }
    return Tool(
        name="search",
        description=(
            "Hybrid BGE-M3 dense+sparse semantic search over the collection's "
            "content_fields. Best for paraphrase / concept queries. "
            "Prefer over grep/filter_rows when the user asks about ideas, "
            "not exact tokens or metadata. Returns hits with score + content."
        ),
        parameters=parameters,
        call=_call,
    )


def _rewrite_tool(rewriter) -> Tool:
    """Expose the query rewriter as an explicit tool the agent calls once.

    Contract: return N candidate retrieval queries (label + query text +
    optional expansions). Does NOT perform search. The agent picks which
    candidates to actually run through the `search` tool in the next round.

    Why a tool and not a subagent: the parent LLM must see each candidate
    as an independent choice so it can (a) skip candidates that obviously
    drift, (b) fan out `search` calls in parallel with different queries
    on round 2, and (c) reason across the returned top-k of each lane.
    See 5ade9c9c probe: the previous in-search fan-out RRF-merged three
    lanes into one top-k, which hid entity-B lane failures from the LLM
    and let it lock onto a wrong assumption.
    """
    def _call(question: str) -> list[dict]:
        try:
            cands = rewriter.rewrite_multi(question)
        except Exception as e:
            return [{"error": f"rewriter failed: {type(e).__name__}: {e}"}]
        out = []
        for c in cands:
            out.append({
                "label": getattr(c, "label", "") or "rewrite",
                "query": getattr(c, "dense_text", "") or "",
                "expansions": list(getattr(c, "expansions", []) or []),
            })
        return out

    return Tool(
        name="rewrite",
        description=(
            "Rewrite the user question into 3-5 complementary English "
            "retrieval query candidates (faithful phrasing, entity-anchored, "
            "concept-broadened). Does NOT search — returns the candidate "
            "list as [{label, query, expansions}, ...]. Call this ONCE at "
            "the start of a new question, then in the next round call "
            "`search` (possibly in parallel) using each candidate's `query` "
            "field verbatim as the `search.question` argument."
        ),
        parameters={
            "type": "object",
            "properties": {
                "question": {
                    "type": "string",
                    "description": "The user's original question, verbatim.",
                },
            },
            "required": ["question"],
        },
        call=_call,
    )


def _grep_tool(memo_root: str) -> Tool:
    def _call(
        pattern: str,
        doc_ids: Optional[list] = None,
        context: int = 3,
        ignore_case: bool = False,
        max_results: int = 50,
        path_root: Optional[str] = None,
    ) -> list:
        return memo_tools.grep(
            pattern=pattern,
            path_root=path_root or memo_root,
            doc_ids=doc_ids,
            context=context,
            ignore_case=ignore_case,
            max_results=max_results,
        )

    parameters = {
        "type": "object",
        "properties": {
            "pattern": {"type": "string",
                        "description": "POSIX ERE regex; passed as argv, no shell interpretation"},
            "doc_ids": {"type": "array", "items": {"type": "string"},
                        "description": "if given, only search {path_root}/{doc_id}.md"},
            "context": {"type": "integer", "default": 3, "minimum": 0, "maximum": 20},
            "ignore_case": {"type": "boolean", "default": False},
            "max_results": {"type": "integer", "default": 50, "minimum": 1, "maximum": 200},
        },
        "required": ["pattern"],
    }
    return Tool(
        name="grep",
        description=(
            "POSIX ERE regex over the on-disk markdown corpus. Use for exact "
            "tokens / LaTeX symbols / numbers that `search` would blur. "
            "Prefer over `run_cmd grep` — this returns structured "
            "{file,line,text} hits."
        ),
        parameters=parameters,
        call=_call,
    )


def _list_memos_tool(memo_root: str) -> Tool:
    def _call(pattern: str = "*.md", path_root: Optional[str] = None) -> dict:
        files = memo_tools.list_memos(
            path_root=path_root or memo_root, pattern=pattern,
        )
        return {"count": len(files), "pattern": pattern, "files": files}

    parameters = {
        "type": "object",
        "properties": {
            "pattern": {"type": "string", "default": "*.md"},
        },
    }
    return Tool(
        name="list_memos",
        description=(
            "List document filenames on disk (fnmatch glob). Use to see what "
            "the disk corpus contains before grep/read_memo. Cheap; call "
            "once per session."
        ),
        parameters=parameters,
        call=_call,
    )


def _read_memo_tool(retriever, memo_root: str, default_collection: Optional[str]) -> Tool:
    def _call(
        doc_id: str,
        collection: Optional[str] = None,
    ) -> dict:
        if collection and not _VALID_COLLECTION_RE.match(collection):
            collection = None
        # Normalize: LLMs often pass the on-disk filename (e.g.
        # "BAM-00713_Memo_SP.md") when the Milvus row's id is a short prefix
        # ("BAM-00713"). Try the raw value first, then strip a `.md` suffix,
        # then keep only a `BAM-\d+` head as a last resort so we still hit the
        # metadata row.
        candidates = [doc_id]
        if doc_id.endswith(".md"):
            candidates.append(doc_id[:-3])
        import re as _re
        m = _re.match(r"^(BAM-\d+)", doc_id)
        if m and m.group(1) not in candidates:
            candidates.append(m.group(1))

        last: dict = {}
        for cand in candidates:
            got = memo_tools.read_memo(
                doc_id=cand,
                retriever=retriever,
                collection=collection or default_collection,
                path_root=memo_root,
            )
            if got.get("text") or got.get("source") in ("milvus", "disk"):
                return got
            last = got
        return last or {"source": "none", "doc_id": doc_id, "text": "", "error": "doc not found"}

    parameters = {
        "type": "object",
        "properties": {
            "doc_id": {"type": "string"},
            "collection": {"type": "string"},
        },
        "required": ["doc_id"],
    }
    return Tool(
        name="read_memo",
        description=(
            "Fetch the full text of one document. Use after search/filter/grep "
            "narrows to a specific doc_id. Milvus parent-chunk join when the "
            "collection is chunked; disk fallback otherwise."
        ),
        parameters=parameters,
        call=_call,
    )


def _get_section_tool(
    retriever,
    default_collection: Optional[str],
    *,
    schema_profile: Any = None,
) -> Tool:
    def _call(doc_id: str, section: str, collection: Optional[str] = None) -> dict:
        if collection and not _VALID_COLLECTION_RE.match(collection):
            collection = None
        coll = collection or default_collection
        # Guard: this collection has no chunk_type/section_type — the underlying
        # query builds a `chunk_type == "parent"` filter that Milvus will reject.
        if schema_profile is not None and coll == getattr(schema_profile, "collection", None):
            if not getattr(schema_profile, "has_chunking", False) or not getattr(schema_profile, "has_section_type", False):
                return {
                    "error": (
                        f"get_section is unavailable on `{coll}`: "
                        "collection has no chunk_type/section_type columns. "
                        "Use `read_memo` or `filter_rows` instead."
                    ),
                    "doc_id": doc_id,
                    "section": section,
                }
        return memo_tools.get_section(
            doc_id=doc_id,
            section=section,
            retriever=retriever,
            collection=coll,
        )

    parameters = {
        "type": "object",
        "properties": {
            "doc_id": {"type": "string"},
            "section": {"type": "string",
                        "description": "e.g. 'systematic_uncertainty', 'method', 'results'"},
            "collection": {"type": "string"},
        },
        "required": ["doc_id", "section"],
    }
    return Tool(
        name="get_section",
        description=(
            "Fetch one section (`section_type`) of a chunked document. "
            "Only usable when the schema line says `has_chunking=true` and "
            "`has_section_type=true`; otherwise returns an error dict and you "
            "should fall back to `read_memo`."
        ),
        parameters=parameters,
        call=_call,
    )


def _filter_rows_tool(
    retriever,
    default_collection: Optional[str],
    *,
    schema_profile: Any = None,
    max_limit: int = 200,
) -> Tool:
    """Milvus scalar-filter passthrough.

    Exposes `MilvusClient.query(filter=..., output_fields=..., limit=...)` to
    the agent, with three guards:
      * vector fields silently stripped from output_fields (they are huge and
        useless to the LLM);
      * `limit` capped at `max_limit`;
      * collection name still restricted to `[A-Za-z0-9_]+`.
    """

    def _call(
        filter_expr: str,
        output_fields: Optional[list] = None,
        limit: int = 20,
        collection: Optional[str] = None,
    ) -> list[dict]:
        if collection and not _VALID_COLLECTION_RE.match(collection):
            collection = None
        coll = collection or default_collection
        if not coll:
            return [{"error": "collection is required and no default configured"}]
        limit = max(1, min(int(limit or 20), max_limit))

        if not output_fields:
            if schema_profile is not None and coll == getattr(schema_profile, "collection", None):
                output_fields = list(getattr(schema_profile, "output_fields", []) or [])
            else:
                output_fields = None
        if output_fields:
            output_fields = [
                f for f in output_fields
                if not any(f.endswith(sfx) or f == sfx for sfx in _VECTOR_FIELD_SUFFIXES)
            ]

        try:
            retriever.store.load_collection(coll)
            rows = retriever.store.client.query(
                collection_name=coll,
                filter=filter_expr,
                output_fields=output_fields or None,
                limit=limit,
            )
        except Exception as e:
            return [{"error": f"{type(e).__name__}: {e}", "filter": filter_expr}]
        # Row dicts may still carry vectors on some servers; drop them.
        cleaned: list[dict] = []
        for r in rows:
            cleaned.append({
                k: v for k, v in r.items()
                if not any(k.endswith(sfx) or k == sfx for sfx in _VECTOR_FIELD_SUFFIXES)
            })
        return cleaned

    parameters = {
        "type": "object",
        "properties": {
            "filter_expr": {
                "type": "string",
                "description": (
                    "Milvus boolean expression over scalar fields. Examples: "
                    "`author like \"%Jane Q. Author%\"`, "
                    "`status == \"PUBLISHED\" and journal like \"%PRL%\"`, "
                    "`id in [\"BAM-00713\", \"BAM-00831\"]`. Only reference "
                    "columns listed in the current-collection schema preamble."
                ),
            },
            "output_fields": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "Scalar fields to return. Defaults to all non-vector "
                    "fields of the current collection. Never request vector "
                    "columns (they are stripped anyway)."
                ),
            },
            "limit": {
                "type": "integer",
                "default": 20,
                "minimum": 1,
                "maximum": max_limit,
            },
            "collection": {
                "type": "string",
                "description": (
                    "Usually leave unset. Only alphanumerics and underscores "
                    "allowed."
                ),
            },
        },
        "required": ["filter_expr"],
    }
    return Tool(
        name="filter_rows",
        description=(
            "Exact Milvus scalar filter over the current collection (no "
            "vectors). Use when the user names a person / referee / id / "
            "status / journal — semantic `search` cannot match scalar payload "
            "fields. Syntax: `col == \"x\"`, `col like \"%x%\"`, "
            "`col in [\"a\",\"b\"]`, combined with `and`/`or`."
        ),
        parameters=parameters,
        call=_call,
    )


def _run_cmd_tool(memo_root: str, *, timeout: int = 20) -> Tool:
    def _call(argv: list, timeout_s: Optional[int] = None) -> dict:
        try:
            return shell_tools.run_cmd(
                argv,
                path_root=memo_root,
                timeout=int(timeout_s) if timeout_s else timeout,
            )
        except shell_tools.MemoToolError as e:
            return {"error": str(e), "argv": argv}

    parameters = {
        "type": "object",
        "properties": {
            "argv": {
                "type": "array",
                "items": {"type": "string"},
                "minItems": 1,
                "description": (
                    "Read-only argv. argv[0] must be one of "
                    "`grep`, `find`, `wc`, `head`, `tail`, `cat`, `ls`, "
                    "`sort`, `uniq`, `cut`, `tr`. Path arguments are resolved "
                    "under the memo corpus root and rejected if they escape. "
                    "Destructive flags (`find -delete`, `find -exec`, "
                    "`sort -o`, `grep --output-file`, `tail -f`, ...) are "
                    "rejected up-front. No shell metacharacters — no pipes, "
                    "no redirects, no `sh -c`."
                ),
            },
            "timeout_s": {
                "type": "integer",
                "default": timeout,
                "minimum": 1,
                "maximum": 60,
            },
        },
        "required": ["argv"],
    }
    return Tool(
        name="run_cmd",
        description=(
            "Read-only shell surface: `argv[0]` in "
            "{grep,find,wc,head,tail,cat,ls,sort,uniq,cut,tr}. Use for file "
            "listings, counts, name globs, quick peeks — anything the "
            "specialised tools don't cover. No pipes/redirects: chain by "
            "issuing multiple calls; do the join step in your reasoning "
            "using the returned stdout. Destructive flags are refused."
        ),
        parameters=parameters,
        call=_call,
    )


def build_default_tools(
    retriever,
    memo_root: str,
    *,
    default_collection: Optional[str] = None,
    schema_profile: Any = None,
    use_subagents: bool = False,
    llm: Any = None,
    hypernews_db_path: Optional[str] = None,
    memo_index_path: Optional[str] = None,
    memo_summaries_root: Optional[str] = None,
    rewriter: Any = None,
    memory: Any = None,
    default_with_rerank: bool = False,
    search_return_top: Optional[int] = None,
) -> list[Tool]:
    """Assemble the parent-agent toolset.

    When `use_subagents=True` and `llm` is provided:
      * `milvus_search` — one Milvus collection, hybrid recall + scalar
        filter + full-doc pull. Schema is compiled into its own prompt.
      * `memo_search`   — on-disk markdown corpus. Grep / list / read /
        run_cmd, all disk-only.
      * `hypernews_db`  — HyperNews SQLite mirror. Only registered when
        `hypernews_db_path` points at an existing file.

    The default (`use_subagents=False`) keeps the flat 7-tool surface for
    backwards compatibility.
    """
    tools: list[Tool] = []

    if use_subagents and llm is not None and default_collection:
        from hepuke.agent.subagents import (
            build_hypernews_db_subagent,
            build_memo_search_subagent,
            build_milvus_search_subagent,
        )
        tools.append(build_milvus_search_subagent(
            llm=llm,
            retriever=retriever,
            collection=default_collection,
            schema_profile=schema_profile,
            rewriter=rewriter,
            memory=memory,
        ))
        tools.append(build_memo_search_subagent(
            llm=llm, memo_root=memo_root,
            memo_index_path=memo_index_path,
            memo_summaries_root=memo_summaries_root,
            rewriter=rewriter,
            memory=memory,
        ))
        if hypernews_db_path:
            import os as _os
            if _os.path.isfile(hypernews_db_path):
                tools.append(build_hypernews_db_subagent(
                    llm=llm, db_path=hypernews_db_path,
                    rewriter=rewriter,
                    memory=memory,
                ))
    else:
        base = [
            _search_tool(retriever, default_collection, default_with_rerank=default_with_rerank,
                         return_top=search_return_top),
            _filter_rows_tool(retriever, default_collection, schema_profile=schema_profile),
            _get_section_tool(retriever, default_collection, schema_profile=schema_profile),
            _grep_tool(memo_root),
            _run_cmd_tool(memo_root),
            _list_memos_tool(memo_root),
            _read_memo_tool(retriever, memo_root, default_collection),
        ]
        # Expose rewrite as a first-class tool when a rewriter is wired in.
        # The parent agent is expected to call rewrite once at round 1
        # (typically forced via tool_choice) and then fan out search calls.
        if rewriter is not None:
            base.insert(0, _rewrite_tool(rewriter))
        tools.extend(base)
    return tools


def _compact(hit: dict) -> dict:
    """Drop vectors + heavy fields from a hit before sending to the model."""
    keep = (
        "score", "doc_id", "bam_id", "paper_id", "chunk_type", "parent_id",
        "section_type", "section_number", "section_title", "title",
        "authors", "date", "source_url", "content", "child_content",
    )
    return {k: hit.get(k) for k in keep if k in hit or hit.get(k) is not None}
