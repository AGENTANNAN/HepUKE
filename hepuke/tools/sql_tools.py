"""Read-only SQL sandbox exposed to the agent.

`sql_query(sql, db_path, ...)` lets a subagent (typically `hypernews_db`)
issue *read-only* SQLite queries with layered guards. Same
"application-layer sandbox" tier as `shell_tools.run_cmd`: no OS-level
containment, but enough guards that a hostile prompt cannot mutate or
exfiltrate.

Guards (all MUST hold; any one failing rejects the call):

  1. **Connection is read-only.** `sqlite3.connect("file:PATH?mode=ro",
     uri=True)` — SQLite opens the file O_RDONLY, so DDL/DML raise
     `OperationalError` before touching any bytes. Even a bypass of the
     statement guard cannot mutate the db.
  2. **Single statement.** We reject anything with a `;` that is not
     trailing whitespace / comment. This is what stops
     `SELECT ...; ATTACH DATABASE ...`.
  3. **Leading keyword whitelist.** After comment/whitespace strip, the
     first token must be `SELECT` or `WITH`. This rejects `PRAGMA`,
     `ATTACH`, `DETACH`, `INSERT`, `UPDATE`, etc. before the parser sees
     them.
  4. **LIMIT injection.** If the statement has no `LIMIT` clause we
     append one (`row_cap`). Prevents 110k-row scans from serialising into
     the LLM context.
  5. **Wall-clock timeout.** `set_progress_handler` aborts long-running
     queries; `sqlite3.connect(timeout=...)` covers lock waits.
  6. **Byte cap on serialised results.** Enforced after row fetch — the
     wrapper truncates the returned row list and flags `truncated=True`
     rather than returning half a JSON blob.

**What this does NOT protect against, by design:**
  * Reading arbitrary tables/columns in the pointed-at db. If the file
    holds secrets, don't point the subagent at it.
  * `ATTACH DATABASE '...' AS x` — blocked by the leading-keyword
    whitelist AND the read-only file handle; both would need to fall for
    an attach to land.
  * Fork bombs / cpu abuse via recursive CTEs — mitigated by the progress
    handler and row cap only, not by a resource cgroup.
"""
from __future__ import annotations

import os
import re
import sqlite3
import time
from typing import Any

from hepuke.observability.trace import emit


class SqlToolError(RuntimeError):
    pass


_COMMENT_LINE = re.compile(r"--[^\n]*")
_COMMENT_BLOCK = re.compile(r"/\*.*?\*/", re.DOTALL)
_LEADING_KEYWORD = re.compile(r"^\s*([A-Za-z]+)")
_HAS_LIMIT = re.compile(r"\blimit\b", re.IGNORECASE)

_ALLOWED_LEAD = frozenset({"SELECT", "WITH"})


def _strip_comments(sql: str) -> str:
    """Remove `-- ...` and `/* ... */` comments.

    Not SQL-perfect (won't handle comment markers inside string literals),
    but the goal isn't to reformat — the goal is to prevent
    `SELECT 1 -- ; DROP` from smuggling a second statement past the
    single-statement guard. After stripping, the guard sees the trailing
    `;` and either accepts (only whitespace after it) or rejects.
    """
    sql = _COMMENT_BLOCK.sub(" ", sql)
    sql = _COMMENT_LINE.sub(" ", sql)
    return sql


def _validate(sql: str) -> str:
    """Return a normalised single SELECT/WITH statement or raise."""
    if not isinstance(sql, str):
        raise SqlToolError(f"sql must be a string, got {type(sql)!r}")
    stripped = _strip_comments(sql).strip()
    if not stripped:
        raise SqlToolError("empty sql")

    # Multi-statement check: split on `;`; every piece past the first must
    # be pure whitespace once comments are gone.
    parts = stripped.split(";")
    core = parts[0].strip()
    for extra in parts[1:]:
        if extra.strip():
            raise SqlToolError(
                "only one statement per call is allowed; multi-statement SQL "
                "is rejected"
            )
    if not core:
        raise SqlToolError("empty sql after comment strip")

    m = _LEADING_KEYWORD.match(core)
    if not m:
        raise SqlToolError("could not identify leading keyword")
    lead = m.group(1).upper()
    if lead not in _ALLOWED_LEAD:
        raise SqlToolError(
            f"leading keyword `{lead}` is not allowed. Only "
            f"{sorted(_ALLOWED_LEAD)} accepted (this is a read-only sandbox)."
        )
    return core


def _inject_limit(sql: str, row_cap: int) -> str:
    """Append `LIMIT row_cap` when the caller didn't specify one."""
    if _HAS_LIMIT.search(sql):
        return sql
    return f"{sql}\nLIMIT {int(row_cap)}"


def _byte_cap_rows(rows: list[dict], max_bytes: int) -> tuple[list[dict], bool]:
    """Truncate the row list once its JSON serialisation would exceed max_bytes."""
    import json as _json
    kept: list[dict] = []
    running = 0
    for r in rows:
        s = _json.dumps(r, ensure_ascii=False, default=str)
        if running + len(s) > max_bytes:
            return kept, True
        kept.append(r)
        running += len(s) + 1  # comma
    return kept, False


def sql_query(
    sql: str,
    *,
    db_path: str | os.PathLike,
    row_cap: int = 500,
    max_bytes: int = 128 * 1024,
    timeout_s: float = 8.0,
) -> dict:
    """Run one read-only SELECT / WITH statement against a SQLite file.

    Parameters
    ----------
    sql : SQL text. Multi-statement input is rejected; the leading keyword
        must be `SELECT` or `WITH`.
    db_path : path to a SQLite database. Opened read-only via `mode=ro`.
    row_cap : if the statement has no `LIMIT`, we inject `LIMIT row_cap`.
        Also the hard cap on rows returned even when a `LIMIT` is present.
    max_bytes : serialised-JSON byte cap on the returned rows.
    timeout_s : wall-clock cap; `sqlite3.connect(timeout=...)` handles lock
        waits, and `set_progress_handler` aborts long-running queries.

    Returns
    -------
    dict
        `{sql, columns, rows, row_count, truncated, elapsed_ms}`.
        `truncated=True` when either the row cap OR the byte cap fired.
        Validation errors raise `SqlToolError` synchronously (never reach
        the DB).
    """
    validated = _validate(sql)
    limited = _inject_limit(validated, row_cap)

    db_path = str(db_path)
    if not os.path.isfile(db_path):
        raise SqlToolError(f"db_path does not exist: {db_path}")

    uri = f"file:{db_path}?mode=ro"

    t0 = time.time()
    deadline = t0 + max(0.5, float(timeout_s))
    conn: sqlite3.Connection | None = None
    try:
        conn = sqlite3.connect(uri, uri=True, timeout=timeout_s)
        conn.row_factory = sqlite3.Row

        def _progress():
            # Called every ~1000 VM instructions. Non-zero return aborts.
            return 1 if time.time() > deadline else 0

        conn.set_progress_handler(_progress, 1000)

        try:
            cursor = conn.execute(limited)
        except sqlite3.OperationalError as e:
            raise SqlToolError(f"OperationalError: {e}") from e
        except sqlite3.DatabaseError as e:
            raise SqlToolError(f"DatabaseError: {e}") from e

        try:
            fetched = cursor.fetchmany(int(row_cap) + 1)
        except sqlite3.OperationalError as e:
            # Progress-handler abort surfaces here.
            raise SqlToolError(f"query aborted (timeout after {timeout_s}s): {e}") from e

        columns = [d[0] for d in (cursor.description or [])]
        row_capped = len(fetched) > row_cap
        rows = [dict(r) for r in fetched[:row_cap]]
        byte_capped_rows, byte_capped = _byte_cap_rows(rows, max_bytes)

        elapsed_ms = round((time.time() - t0) * 1000, 1)
        result = {
            "sql": limited,
            "columns": columns,
            "row_count": len(byte_capped_rows),
            "truncated": row_capped or byte_capped,
            "elapsed_ms": elapsed_ms,
            "rows": byte_capped_rows,
        }
        emit(
            "tool_call", name="sql_query",
            row_count=len(byte_capped_rows),
            truncated=result["truncated"],
            elapsed_ms=elapsed_ms,
        )
        return result
    finally:
        if conn is not None:
            conn.close()


def describe_schema(
    *,
    db_path: str | os.PathLike,
    include_counts: bool = False,
    timeout_s: float = 4.0,
) -> dict:
    """Return `{tables: [{name, sql, columns, row_count?}, ...]}`.

    Small helper that lets a subagent bootstrap its schema knowledge without
    burning `sql_query` calls on `sqlite_master`.
    """
    db_path = str(db_path)
    if not os.path.isfile(db_path):
        raise SqlToolError(f"db_path does not exist: {db_path}")
    uri = f"file:{db_path}?mode=ro"
    conn = sqlite3.connect(uri, uri=True, timeout=timeout_s)
    try:
        rows = conn.execute(
            "SELECT name, sql FROM sqlite_master "
            "WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
        ).fetchall()
        tables: list[dict] = []
        for name, ddl in rows:
            cols = conn.execute(f'PRAGMA table_info("{name}")').fetchall()
            entry: dict[str, Any] = {
                "name": name,
                "sql": ddl,
                "columns": [
                    {"name": c[1], "type": c[2], "notnull": bool(c[3]),
                     "pk": bool(c[5])}
                    for c in cols
                ],
            }
            if include_counts:
                try:
                    entry["row_count"] = conn.execute(
                        f'SELECT COUNT(*) FROM "{name}"'
                    ).fetchone()[0]
                except sqlite3.OperationalError:
                    entry["row_count"] = None
            tables.append(entry)
        return {"tables": tables}
    finally:
        conn.close()


__all__ = ["sql_query", "describe_schema", "SqlToolError"]
