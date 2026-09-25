"""`hypernews_db` subagent — specialist over the HyperNews SQLite mirror.

Owns the HyperNews SQLite mirror named by `agent.hypernews_db_path`
(not distributed with this release). Four
tables (joins on `paper_id`):

  * `paper`       (~1.2k) : per-paper metadata + `referee_names_json`,
                            `work_description`, `latest_memo_local_path`.
  * `email`       (~110k) : threaded HyperNews replies. `raw_content`,
                            `sender_role`, `depth`, `time_str`, extracted
                            `core_issue` / `details`.
  * `pdfdocument` (~1.5k) : PDFs attached to threads. Parsed markdown in
                            `markdown_content`.
  * `reviewqa`    (~4.5k) : extracted Q&A pairs (`question_text`,
                            `answer_text`, `question_author`,
                            `answer_author`, `phase`, `round_number`).

Two tools:
  * `sql_query`       : read-only SELECT / WITH. Injects LIMIT, byte-caps.
  * `describe_schema` : bootstrap table/column listing so the subagent
                        doesn't waste turns on `sqlite_master`.
"""
from __future__ import annotations

from typing import Any, Optional

from hepuke.agent.subagents.base import make_subagent_tool
from hepuke.agent.tools import Tool
from hepuke.tools import sql_tools


HYPERNEWS_DB_SYSTEM_PROMPT = """\
You are the HyperNews SQLite specialist. The database holds the full
history of referee-author discussion around each BESIII paper. Two tools;
nothing else exists.

Schema (SQLite; join on `paper_id`; author names live in author columns):

  paper (id, paper_id "BAM-xxxxx", url, is_fully_parsed, last_updated,
         latest_memo_local_path, referee_names_json TEXT '[...]',
         work_description TEXT)
    -- `paper_id` is the human id ("BAM-00713"); `id` is the integer FK.

  email (id, paper_id FK, url, author, time_str, depth, raw_content,
         pdf_links_json, is_extracted, core_issue, details,
         sender_role /* 'referee'|'author'|'unknown' */, is_indexed)
    -- threaded HyperNews replies. Depth 0 = original post.

  pdfdocument (id, email_id FK, paper_id FK, name, url, pdf_type,
               version_hint, local_path, markdown_content,
               is_downloaded, is_parsed, is_indexed)
    -- PDFs attached to threads; `markdown_content` is the parsed body.

  reviewqa (id, paper_id FK, question_text, answer_text, phase,
            round_number, question_author, answer_author, source_url,
            memo_excerpt, memo_topic, rationale_for_question,
            answer_source /* 'email' | 'memo' | ... */, extracted_at)
    -- extracted Q&A pairs. Cheaper first hop than trawling `email`.

Routing rules:
  * Author lookups ("what has X worked on"): join through `paper` +
    `reviewqa` or `email` on `paper_id`, filter on the author column with
    `LIKE '%X%'`.
  * Discussion topics per paper: query `reviewqa` first — it's the
    distilled form. Only fall back to `email.raw_content` when you need
    the exact wording.
  * Referee names for a paper: `paper.referee_names_json` is a JSON
    array; use `json_each` (`SELECT value FROM paper, json_each(paper.
    referee_names_json) WHERE paper.paper_id='BAM-00713'`).
  * Always project a small set of columns — never `SELECT *` on `email`
    (rows are large: `raw_content`, `details`, `core_issue`).
  * Always add `LIMIT` (the sandbox will inject 500 if you forget).
  * Prefer `describe_schema` once, THEN sql_query. Do not probe
    `sqlite_master` yourself.

Cost / stop conditions — read carefully:
  * BUDGET: at most 4 sql_query calls per question. If you have not
    answered after 4 queries, STOP and return `dry=true`.
  * ZERO-RESULT RULE: if `SELECT ... FROM paper WHERE paper_id = 'BAM-XXX'`
    returns 0 rows, that paper is NOT in the database. Return `dry=true`
    IMMEDIATELY. Do NOT try synonyms, variants, LIKE patterns, or other
    tables — the paper simply isn't here.
  * DO NOT probe multiple tables looking for a paper that doesn't exist in
    `paper`. The four tables are all keyed off `paper.id`, so if `paper`
    doesn't have it, nothing else will either.
  * If a query fails with an error, fix it and retry ONCE. Two consecutive
    errors ⇒ return `dry=true`.

Sandbox limits (hard):
  * Only `SELECT` or `WITH` (CTE). Multi-statement, `ATTACH`, `PRAGMA`,
    `INSERT`, `UPDATE`, `DELETE`, `CREATE` are refused. The file is opened
    read-only regardless.
  * Wall-clock ≤ 8 s. Row cap 500. Byte cap 128 KiB.

Contract — your FINAL message MUST be a JSON object, no prose, no fences:
  {
    "answer":   "<one-paragraph natural-language answer, or null>",
    "evidence": [ {"source": "hypernews_db", "paper_id": "BAM-...",
                   "table": "reviewqa", "snippet": "..."}, ... ],
    "dry":      true|false
  }
Every claim in `answer` must trace to a `paper_id` in `evidence`. Never
invent paper_ids. Set `dry=true` when every query returned zero rows so
the parent can route to a sibling corpus.
"""


def _sql_query_tool(db_path: str, *, row_cap: int, max_bytes: int, timeout_s: float) -> Tool:
    def _call(sql: str, row_cap_override: Optional[int] = None) -> dict:
        try:
            return sql_tools.sql_query(
                sql,
                db_path=db_path,
                row_cap=int(row_cap_override) if row_cap_override else row_cap,
                max_bytes=max_bytes,
                timeout_s=timeout_s,
            )
        except sql_tools.SqlToolError as e:
            return {"error": str(e), "sql": sql}

    parameters = {
        "type": "object",
        "properties": {
            "sql": {
                "type": "string",
                "description": (
                    "One SELECT or WITH statement. Multi-statement input is "
                    "rejected. `LIMIT` will be injected if you omit it. "
                    "The file is opened read-only — DDL/DML would fail even "
                    "if the guard missed."
                ),
            },
            "row_cap_override": {
                "type": "integer",
                "description": (
                    "Lower the row cap for exploratory queries (default 500)."
                ),
                "minimum": 1,
                "maximum": row_cap,
            },
        },
        "required": ["sql"],
    }
    return Tool(
        name="sql_query",
        description=(
            "Run one read-only SELECT / WITH against the HyperNews SQLite "
            "database. Returns {columns, rows, row_count, truncated, "
            "elapsed_ms}. Byte-capped at 128 KiB; row-capped at 500."
        ),
        parameters=parameters,
        call=_call,
    )


def _describe_schema_tool(db_path: str) -> Tool:
    def _call(include_counts: bool = False) -> dict:
        try:
            return sql_tools.describe_schema(
                db_path=db_path, include_counts=include_counts,
            )
        except sql_tools.SqlToolError as e:
            return {"error": str(e)}

    parameters = {
        "type": "object",
        "properties": {
            "include_counts": {
                "type": "boolean", "default": False,
                "description": "add per-table `row_count`; costs one SELECT COUNT(*) per table",
            },
        },
    }
    return Tool(
        name="describe_schema",
        description=(
            "Return the HyperNews SQLite schema — table names, DDL, column "
            "types. Call once at the start; don't re-probe."
        ),
        parameters=parameters,
        call=_call,
    )


def build_hypernews_db_tools(
    db_path: str,
    *,
    row_cap: int = 500,
    max_bytes: int = 128 * 1024,
    timeout_s: float = 8.0,
) -> list[Tool]:
    return [
        _describe_schema_tool(db_path),
        _sql_query_tool(
            db_path, row_cap=row_cap, max_bytes=max_bytes, timeout_s=timeout_s,
        ),
    ]


HYPERNEWS_DB_DESCRIPTION = (
    "Query the HyperNews SQLite database: per-paper metadata, referee/"
    "author threaded discussion, attached PDFs, and extracted Q&A pairs. "
    "Use for author activity across papers, referee identity, discussion "
    "history, or Q&A around a specific paper_id. Read-only SELECT/WITH "
    "only. NOT a source of memo body text (use `memo_search`) and NOT a "
    "source of vector similarity (use `milvus_search`)."
)


def build_hypernews_db_subagent(
    *,
    llm: Any,
    db_path: str,
    max_steps: int = 5,
    row_cap: int = 500,
    max_bytes: int = 128 * 1024,
    timeout_s: float = 8.0,
    rewriter: Any = None,
    memory: Any = None,
) -> Tool:
    """Factory: return the hypernews_db subagent as one Tool the parent calls.

    `rewriter` — see `make_subagent_tool`. Note hypernews_db mostly runs
    SQL against structured columns, so the rewriter's expansions may or
    may not help; keeping the same wiring uniform means the caller
    doesn't have to think about it.
    """
    return make_subagent_tool(
        name="hypernews_db",
        description=HYPERNEWS_DB_DESCRIPTION,
        llm=llm,
        tools=build_hypernews_db_tools(
            db_path, row_cap=row_cap, max_bytes=max_bytes, timeout_s=timeout_s,
        ),
        system_prompt=HYPERNEWS_DB_SYSTEM_PROMPT,
        max_steps=max_steps,
        rewriter=rewriter,
        memory=memory,
    )


__all__ = [
    "build_hypernews_db_subagent",
    "build_hypernews_db_tools",
    "HYPERNEWS_DB_SYSTEM_PROMPT",
    "HYPERNEWS_DB_DESCRIPTION",
]
