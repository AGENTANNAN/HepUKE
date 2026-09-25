"""`milvus_search` subagent — specialist over one Milvus collection.

Owns the currently-connected collection (typically `hypernews_hybrid` or a
chunked P0 layout). Three tools; the subagent knows the schema shape and
picks between them:

  * `retrieve_vec`  → BGE-M3 hybrid dense+sparse over `content_fields`.
                      Right for paraphrase / concept queries.
  * `filter_rows`   → exact scalar filter over `output_fields`. Right when
                      the user names a person / status / journal / id —
                      semantic search cannot match scalar payload.
  * `get_full_doc`  → parent-chunk join (or single-row fetch on flat
                      metadata collections) for a known doc_id.

The main agent talks to this subagent through a single Tool. The
per-collection schema (which columns are vectorised vs scalar, whether the
layout is chunked) is compiled into the subagent's system prompt at build
time, so the main agent doesn't need to remember any of it.
"""
from __future__ import annotations

import json
import re
from typing import Any, Optional

from hepuke.agent.subagents.base import make_subagent_tool
from hepuke.agent.tools import Tool


_VALID_COLLECTION_RE = re.compile(r"^[A-Za-z0-9_]+$")
_VECTOR_FIELD_SUFFIXES = ("_vector", "_vec", "embedding", "embeddings")


# ---------------------------------------------------------------------------
# tools (each closure captures `retriever` + `collection` — the subagent is
# pinned to ONE collection; cross-collection routing is the parent's job)
# ---------------------------------------------------------------------------

def _retrieve_vec_tool(retriever, collection: str, *, rewriter=None) -> Tool:
    def _call(
        question: str,
        top_k: int = 5,
        section_types: Optional[list] = None,
        doc_id: Optional[str] = None,
        with_rerank: bool = False,
    ) -> list:
        # Multi-query rewrite: when a rewriter is wired, let it produce
        # several English candidate queries (faithful / concept / entity-
        # expansion) so an EN corpus is searched under the term family a
        # paper would actually use — not just the string the user typed.
        # `retriever.search` fans out per pair and RRF-merges. If the
        # rewriter fails, its own fallback yields ONE candidate, so this
        # path degrades gracefully to legacy single-query behaviour.
        query_pairs = None
        if rewriter is not None:
            try:
                cands = rewriter.rewrite_multi(question)
                query_pairs = [c.as_pair() for c in cands]
            except Exception:
                query_pairs = None
        hits = retriever.search(
            question=question,
            collection=collection,
            top_k=top_k,
            section_types=section_types,
            doc_id=doc_id,
            with_rerank=with_rerank,
            query_pairs=query_pairs,
        )
        return [_compact(h) for h in hits]

    parameters = {
        "type": "object",
        "properties": {
            "question": {"type": "string"},
            "top_k": {"type": "integer", "default": 5, "minimum": 1, "maximum": 20},
            "section_types": {
                "type": "array", "items": {"type": "string"},
                "description": "whitelist of section_type values (chunked layouts only)",
            },
            "doc_id": {
                "type": "string",
                "description": "restrict search to one document id",
            },
            "with_rerank": {"type": "boolean", "default": False},
        },
        "required": ["question"],
    }
    return Tool(
        name="retrieve_vec",
        description=(
            "BGE-M3 hybrid dense+sparse over the collection's vectorised "
            "content fields. Use for concept / paraphrase queries. Cannot "
            "match scalar payload — use `filter_rows` for those."
        ),
        parameters=parameters,
        call=_call,
    )


def _filter_rows_tool(
    retriever,
    collection: str,
    *,
    schema_profile: Any,
    max_limit: int = 200,
) -> Tool:
    def _call(
        filter_expr: str,
        output_fields: Optional[list] = None,
        limit: int = 20,
    ) -> list[dict]:
        n = max(1, min(int(limit or 20), max_limit))
        if not output_fields:
            output_fields = list(getattr(schema_profile, "output_fields", []) or [])
        if output_fields:
            output_fields = [
                f for f in output_fields
                if not any(f.endswith(sfx) or f == sfx for sfx in _VECTOR_FIELD_SUFFIXES)
            ]
        try:
            retriever.store.load_collection(collection)
            rows = retriever.store.client.query(
                collection_name=collection,
                filter=filter_expr,
                output_fields=output_fields or None,
                limit=n,
            )
        except Exception as e:
            # Repeat the doc_id_field + allowed column list in the error so
            # the LLM's next attempt is guided by the SAME schema info the
            # system prompt tried to embed — LLMs frequently ignore prompt
            # rules but always read tool-error text. This is what turns
            # "pk == \"BAM-00831\"" (wrong column) into a self-heal on the
            # next round.
            doc_id_field = getattr(schema_profile, "doc_id_field", "doc_id") \
                if schema_profile is not None else "doc_id"
            allowed = list(getattr(schema_profile, "output_fields", []) or []) \
                if schema_profile is not None else []
            hint = (
                f" — reminder: the document identifier column in this "
                f"collection is `{doc_id_field}` (VARCHAR, always quoted). "
                f"Allowed columns: {allowed}. Only reference these; "
                f"`pk` / `doc_id` / `paper_id` / `bam_id` may not exist here."
            )
            return [{
                "error": f"{type(e).__name__}: {e}{hint}",
                "filter": filter_expr,
                "doc_id_field": doc_id_field,
                "allowed_columns": allowed,
            }]
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
                    "`id in [\"BAM-00713\", \"BAM-00831\"]`."
                ),
            },
            "output_fields": {
                "type": "array", "items": {"type": "string"},
                "description": (
                    "Scalar fields to return. Defaults to all non-vector "
                    "fields for this collection."
                ),
            },
            "limit": {
                "type": "integer", "default": 20,
                "minimum": 1, "maximum": max_limit,
            },
        },
        "required": ["filter_expr"],
    }
    return Tool(
        name="filter_rows",
        description=(
            "Exact Milvus scalar filter. Use when the user names a person / "
            "id / status / journal — semantic search cannot match scalar "
            "payload. Syntax: `col == \"x\"`, `col like \"%x%\"`, "
            "`col in [...]`, combined with `and`/`or`."
        ),
        parameters=parameters,
        call=_call,
    )


def _get_full_doc_tool(retriever, collection: str) -> Tool:
    def _call(doc_id: str, section_order: Optional[list] = None) -> dict:
        try:
            got = retriever.get_full_doc(
                doc_id, collection=collection, section_order=section_order,
            )
        except Exception as e:
            return {"error": f"{type(e).__name__}: {e}", "doc_id": doc_id}
        if not got:
            return {"doc_id": doc_id, "sections": [], "error": "doc not found"}
        # Reduce serialised payload for the LLM: keep metadata + concatenated
        # section text.
        body = "\n\n".join(
            f"## {s.get('section_title') or s.get('section_type', '')}\n"
            f"{s.get('content', '')}"
            for s in got.get("sections", [])
        )
        return {
            "doc_id": got.get("doc_id", doc_id),
            "title": got.get("title"),
            "authors": got.get("authors"),
            "date": got.get("date"),
            "n_sections": len(got.get("sections", [])),
            "text": body[:60_000],
        }

    parameters = {
        "type": "object",
        "properties": {
            "doc_id": {"type": "string"},
            "section_order": {
                "type": "array", "items": {"type": "string"},
                "description": "optional order over section_type values",
            },
        },
        "required": ["doc_id"],
    }
    return Tool(
        name="get_full_doc",
        description=(
            "Fetch one document's parent chunks joined into a single text "
            "blob. On flat metadata collections returns the matching row's "
            "content. Use after retrieve_vec/filter_rows narrows to a doc_id."
        ),
        parameters=parameters,
        call=_call,
    )


def build_milvus_search_tools(
    retriever, collection: str, *, schema_profile: Any, rewriter=None,
) -> list[Tool]:
    """Three-tool bundle for the `milvus_search` subagent."""
    return [
        _retrieve_vec_tool(retriever, collection, rewriter=rewriter),
        _filter_rows_tool(retriever, collection, schema_profile=schema_profile),
        _get_full_doc_tool(retriever, collection),
    ]


# ---------------------------------------------------------------------------
# system prompt — schema is compiled in at build time so the parent doesn't
# have to carry it.
# ---------------------------------------------------------------------------

_SYSTEM_PROMPT_TEMPLATE = """\
You are the Milvus-search specialist over collection `{collection}`. You
have three tools; nothing else exists.

Schema (frozen at connect time — DO NOT reference other columns):
{schema_block}

Routing rules:
  * `retrieve_vec` — semantic / paraphrase / concept queries. Only matches
                     the vectorised content fields above.
  * `filter_rows`  — exact scalar filter (author, status, journal, id …).
                     Semantic search CANNOT find scalar payload; if the
                     user names a person or metadata value, start here.
                     Syntax: `col == "x"`, `col like "%x%"`, `col in [...]`,
                     combined with `and` / `or`. Never reference a column
                     that isn't listed in Schema above.
                     **The document identifier column is `{doc_id_field}` —
                     it is NOT called `doc_id`, `pk`, `paper_id`, or
                     `bam_id` unless the Schema above says so. Referencing
                     the wrong column name is the #1 way this tool fails.**
                     **`{doc_id_field}` is a VARCHAR — always quote its
                     values: `{doc_id_field} == "BAM-00831"`, never
                     `{doc_id_field} == BAM-00831`. `like` only works on
                     VARCHAR columns.**
  * `get_full_doc` — after you have a concrete doc_id, pull the full body
                     (parent-chunk join on chunked layouts; single row
                     otherwise).

Cost / budget rules:
  * Do not repeat a tool call with identical arguments.
  * If `retrieve_vec` and `filter_rows` both return empty, stop — that is
    the signal to report `dry=true`. The parent will route to a sibling
    corpus.

Contract — your FINAL message MUST be a JSON object, no prose, no fences:
  {{
    "answer":   "<one-paragraph natural-language answer, or null>",
    "evidence": [ {{"source": "milvus", "collection": "{collection}",
                    "doc_id": "...", "snippet": "..."}}, ... ],
    "dry":      true|false
  }}
Every claim in `answer` must trace to a `doc_id` in `evidence`. If you
cannot ground it, drop the claim. Never invent doc_ids.
"""


def _build_schema_block(schema_profile: Any) -> str:
    if schema_profile is None:
        return "  (schema unknown — collection was not profiled at connect time)"
    content = list(getattr(schema_profile, "content_fields", []) or [])
    output = list(getattr(schema_profile, "output_fields", []) or [])
    doc_id = getattr(schema_profile, "doc_id_field", "doc_id")
    scalars = [f for f in output if f not in content and f != doc_id]
    return (
        f"  doc_id_field    : {doc_id}\n"
        f"  content_fields  : {content}   # vectorised → use retrieve_vec\n"
        f"  scalar_fields   : {scalars}   # payload → use filter_rows\n"
        f"  has_chunking    : {bool(getattr(schema_profile, 'has_chunking', False))}\n"
        f"  has_section_type: {bool(getattr(schema_profile, 'has_section_type', False))}"
    )


def build_milvus_search_system_prompt(collection: str, schema_profile: Any) -> str:
    doc_id_field = (
        getattr(schema_profile, "doc_id_field", "doc_id")
        if schema_profile is not None
        else "doc_id"
    )
    return _SYSTEM_PROMPT_TEMPLATE.format(
        collection=collection,
        schema_block=_build_schema_block(schema_profile),
        doc_id_field=doc_id_field,
    )


MILVUS_SEARCH_DESCRIPTION = (
    "Search a Milvus collection: hybrid vector recall over content fields, "
    "exact scalar filter over payload fields (author / status / journal / "
    "id), or full-document fetch by doc_id. Use for author / metadata "
    "lookups and semantic recall. NOT for reading memo body markdown "
    "(use `memo_search` for that)."
)


def build_milvus_search_subagent(
    *,
    llm: Any,
    retriever,
    collection: str,
    schema_profile: Any,
    max_steps: int = 6,
    rewriter=None,
    memory=None,
) -> Tool:
    """Factory: return the milvus-search subagent as one Tool the parent calls.

    `rewriter`, when set, is wired at the SUBAGENT ENTRY (make_subagent_tool)
    — one rewrite fans out into N parallel inner runs, each with its own
    English candidate as the user prompt. The inner `retrieve_vec` tool
    is intentionally built WITHOUT a rewriter to avoid double fan-out
    (N × M Milvus queries instead of N + M). `memory` is passed to the
    wrapper for dry-run dead_end recording.
    """
    return make_subagent_tool(
        name="milvus_search",
        description=MILVUS_SEARCH_DESCRIPTION,
        llm=llm,
        tools=build_milvus_search_tools(
            retriever, collection,
            schema_profile=schema_profile,
            rewriter=None,
        ),
        system_prompt=build_milvus_search_system_prompt(collection, schema_profile),
        max_steps=max_steps,
        memory=memory,
        rewriter=rewriter,
    )


def _compact(hit: dict) -> dict:
    keep = (
        "score", "doc_id", "id", "bam_id",
        "chunk_type", "parent_id", "section_type",
        "section_number", "section_title",
        "title", "authors", "author", "date", "journal", "status",
        "source_url", "content", "child_content",
    )
    return {k: hit.get(k) for k in keep if k in hit or hit.get(k) is not None}


__all__ = [
    "build_milvus_search_subagent",
    "build_milvus_search_tools",
    "build_milvus_search_system_prompt",
    "MILVUS_SEARCH_DESCRIPTION",
]
