"""Trace sinks — pluggable writers for trace records.

Two sinks ship:
  * ``JsonlSink`` — legacy append-JSONL (one line per record).
  * ``SqliteSink`` — persistent SQLite database with indexed columns for the
    common filter dimensions (``run_id`` / ``kind`` / ``case_id`` /
    ``suite`` / ``config``) and the full record kept as JSON in ``payload``.

Both are thread-safe.

The SQLite database is designed to *accumulate* across runs. Each
``benchmarks.runner.run_suite`` invocation stamps a ``run_id`` into the trace
context, so a single DB file can hold history for many suites/configs and
the runner's eval pass queries with ``WHERE run_id = ?``.
"""
from __future__ import annotations

import json
import sqlite3
import sys
import threading
import time
from pathlib import Path
from typing import Any, Iterable, Iterator, Protocol


class TraceSink(Protocol):
    def write(self, record: dict[str, Any]) -> None: ...
    def close(self) -> None: ...


# ---------------------------------------------------------------------------
# JSONL
# ---------------------------------------------------------------------------

class JsonlSink:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def write(self, record: dict[str, Any]) -> None:
        line = json.dumps(record, ensure_ascii=False, default=str)
        with self._lock:
            with self.path.open("a", encoding="utf-8") as f:
                f.write(line + "\n")

    def close(self) -> None:  # pragma: no cover - nothing to release
        pass


# ---------------------------------------------------------------------------
# Console (human-readable one-liner per event, depth-indented)
# ---------------------------------------------------------------------------

# ANSI colors — kept minimal. Disabled automatically when stderr is not a tty.
_COLOR = {
    "reset": "\x1b[0m",
    "dim":   "\x1b[2m",
    "bold":  "\x1b[1m",
    "cyan":  "\x1b[36m",
    "green": "\x1b[32m",
    "yellow":"\x1b[33m",
    "red":   "\x1b[31m",
    "magenta":"\x1b[35m",
    "blue":  "\x1b[34m",
}


# Fixed prefix tag + colour per subagent NAME. When 3 subagents fan out
# concurrently their log lines interleave; a stable name → tag/colour
# makes it possible to read at a glance. Names outside this table fall
# back to a hashed slot below.
_SUB_TAGS: dict[str, tuple[str, str]] = {
    "milvus_search":  ("MILV", "cyan"),
    "memo_search":    ("MEMO", "yellow"),
    "hypernews_db":   ("HN",   "magenta"),
}

_FALLBACK_COLORS = ("green", "blue", "red", "magenta", "cyan", "yellow")


def _sub_tag(name: str | None) -> tuple[str, str] | None:
    """Return (tag, colour) for a subagent name. None when we're at the
    outer agent (no subagent context stamped)."""
    if not name:
        return None
    if name in _SUB_TAGS:
        return _SUB_TAGS[name]
    # Deterministic fallback for unknown names: 4-char upper prefix +
    # hashed colour slot. Stable across ticks so the same name always
    # picks the same colour in a session.
    tag = name.upper().replace("_", "")[:4] or "SUB"
    colour = _FALLBACK_COLORS[hash(name) % len(_FALLBACK_COLORS)]
    return tag, colour


def _c(name: str, s: str, *, enabled: bool) -> str:
    return f"{_COLOR[name]}{s}{_COLOR['reset']}" if enabled else s


class ConsoleSink:
    """Compact console renderer for live agent inspection.

    Each event → one line, `depth`-indented, colored by kind. Verbose kinds
    (retrieve, embed_call, low-level tool_call variants) are suppressed by
    default so the stream stays readable during a real query.

    Only meant for interactive `hepuke serve` — the JSONL / SQLite sinks
    keep the full record for later analysis.
    """

    _VISIBLE = {
        "router_decision",
        "router_fallback",
        "router_error",
        "agent_start",
        "agent_end",
        "subagent_start",
        "subagent_end",
        "agent_tool_call",
        "llm_call",
    }

    def __init__(self, *, stream=None, color: bool | None = None) -> None:
        self.stream = stream or sys.stderr
        if color is None:
            color = hasattr(self.stream, "isatty") and self.stream.isatty()
        self.color = bool(color)
        self._lock = threading.Lock()

    def close(self) -> None:  # pragma: no cover
        pass

    def write(self, record: dict[str, Any]) -> None:
        kind = record.get("kind", "")
        if kind not in self._VISIBLE:
            return
        depth = int(record.get("depth", 0) or 0)
        indent = "  " * depth
        line = self._format(record, kind, indent)
        if line is None:
            return
        # Subagent tag: prepended to EVERY line emitted while the trace
        # context has `subagent_name` set. Three concurrent subagents'
        # lines interleave in the console; without a stable prefix
        # they're unreadable. The tag replaces the leading indent so
        # depth is still visible via horizontal position.
        #
        # `subagent_start` / `subagent_end` are the boundary events —
        # they belong to the PARENT agent (they announce which sub is
        # being called), so they don't get the sub tag even if the
        # emit happens to fire under a stamped context.
        if kind not in {"subagent_start", "subagent_end"}:
            sub_name = record.get("subagent_name")
            tag_pair = _sub_tag(sub_name) if sub_name else None
            if tag_pair is not None:
                tag, colour = tag_pair
                tag_str = _c(colour, f"[{tag}]", enabled=self.color)
                line = f"{tag_str} {line.lstrip()}"
        with self._lock:
            try:
                self.stream.write(line + "\n")
                self.stream.flush()
            except Exception:  # pragma: no cover
                pass

    # -- formatters --------------------------------------------------------

    def _short(self, aid: str | None) -> str:
        return (aid or "--------")[:6]

    def _format(self, r: dict, kind: str, indent: str) -> str | None:
        aid = self._short(r.get("agent_id"))
        c = self.color

        if kind == "router_decision":
            mode = r.get("mode", "?")
            n_sub = r.get("n_sub", 0)
            rat = (r.get("rationale") or "")[:70]
            head = _c("magenta", f"[router] mode={mode}", enabled=c)
            # `n_sub` is only meaningful for decompose; suppress otherwise.
            sub_part = f"n_sub={n_sub} · " if mode == "decompose" else ""
            tail = _c("dim", f"{sub_part}{rat}", enabled=c)
            return f"{indent}{head} {tail}"

        if kind == "router_fallback":
            return (
                f"{indent}"
                + _c("yellow", "[router] fallback → rag", enabled=c)
                + " "
                + _c("dim", f"reason={r.get('reason','?')}", enabled=c)
            )

        if kind == "router_error":
            return (
                f"{indent}"
                + _c("red", "[router] error", enabled=c)
                + " "
                + _c("dim", str(r.get("error", ""))[:120], enabled=c)
            )

        if kind == "agent_start":
            q = (r.get("question") or "").replace("\n", " ")[:200]
            mode = r.get("mode")
            depth = int(r.get("depth", 0) or 0)
            role = "AGENT" if depth == 0 else "sub"
            head = _c("cyan", f"[{role} {aid}] 目标", enabled=c)
            extra = f" mode={mode}" if mode else ""
            return f"{indent}{head}{extra} " + _c("dim", f'"{q}"', enabled=c)

        if kind == "agent_end":
            steps = r.get("steps", "?")
            alen = r.get("answer_len", 0)
            dt = r.get("latency_ms")
            depth = int(r.get("depth", 0) or 0)
            role = "AGENT" if depth == 0 else "sub"
            terminated = r.get("terminated")
            marker = _c("red", " [MAX_STEPS]", enabled=c) if terminated else ""
            tail = f"steps={steps} answer={alen}b"
            if dt:
                tail += f" · {dt}ms"
            return (
                f"{indent}"
                + _c("cyan", f"[{role} {aid}] end", enabled=c)
                + f"{marker} "
                + _c("dim", tail, enabled=c)
            )

        if kind == "subagent_start":
            name = r.get("name", "?")
            q = (r.get("question") or "").replace("\n", " ").strip()
            hint = r.get("doc_id_hint")
            body = q[:160] if q else ""
            hint_part = f" +hint={hint}" if hint else ""
            head = _c("blue", f"  → 派 {name}{hint_part}", enabled=c)
            if body:
                return f"{indent}{head} " + _c("dim", f'"{body}"', enabled=c)
            return f"{indent}{head}"

        if kind == "subagent_end":
            name = r.get("name", "?")
            dry = r.get("dry", False)
            steps = r.get("steps", 0)
            evi = r.get("evidence_n", 0)
            color = "yellow" if dry else "green"
            status = "DRY" if dry else "ok"
            return (
                f"{indent}"
                + _c(color, f"  ← {name} {status}", enabled=c)
                + " "
                + _c("dim", f"steps={steps} evidence={evi}", enabled=c)
            )

        if kind == "agent_tool_call":
            tool = r.get("tool", "?")
            dt = r.get("latency_ms", "?")
            err = r.get("error")
            summary = r.get("args_summary") or ""
            head_c = "red" if err else "green"
            head = _c(head_c, f"  · tool {tool}", enabled=c)
            tail = f"{dt}ms" + (f" ERR={err[:40]}" if err else "")
            if summary:
                tail += f" | {summary}"
            return f"{indent}{head} " + _c("dim", tail, enabled=c)

        if kind == "llm_call":
            pt = r.get("prompt_tokens", 0) or 0
            ct = r.get("cached_tokens", 0) or 0
            comp = r.get("completion_tokens", 0) or 0
            hit = f" ({int(ct/pt*100)}% cached)" if pt and ct else ""
            return (
                f"{indent}"
                + _c("dim", f"  · llm prompt={pt} completion={comp}{hit}", enabled=c)
            )
        return None


# ---------------------------------------------------------------------------
# SQLite
# ---------------------------------------------------------------------------

_SCHEMA = """
CREATE TABLE IF NOT EXISTS events (
    id       INTEGER PRIMARY KEY AUTOINCREMENT,
    ts       REAL NOT NULL,
    kind     TEXT NOT NULL,
    run_id   TEXT,
    case_id  TEXT,
    suite    TEXT,
    config   TEXT,
    payload  TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS ix_events_run   ON events(run_id);
CREATE INDEX IF NOT EXISTS ix_events_kind  ON events(kind, ts);
CREATE INDEX IF NOT EXISTS ix_events_case  ON events(case_id);
CREATE INDEX IF NOT EXISTS ix_events_suite ON events(suite, config);

-- Per-hit / per-signal rows extracted from tool results.
-- One row per doc_id (or row-count signal) returned by one tool call.
-- Compact so aggregations (recall@K / entropy / rank stability) stay cheap.
CREATE TABLE IF NOT EXISTS retrievals (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    ts             REAL NOT NULL,
    session_id     TEXT,
    run_id         TEXT,
    agent_id       TEXT,
    subagent_name  TEXT,
    step           INTEGER,
    tool           TEXT,
    doc_id         TEXT,
    rank           INTEGER,
    score          REAL,
    row_count      INTEGER,
    snippet        TEXT,
    raw_id         TEXT
);
CREATE INDEX IF NOT EXISTS ix_retr_session ON retrievals(session_id);
CREATE INDEX IF NOT EXISTS ix_retr_run     ON retrievals(run_id);
CREATE INDEX IF NOT EXISTS ix_retr_agent   ON retrievals(agent_id);
CREATE INDEX IF NOT EXISTS ix_retr_doc     ON retrievals(doc_id);
CREATE INDEX IF NOT EXISTS ix_retr_raw     ON retrievals(raw_id);

-- Full raw tool result, for audit / replay. Keyed by raw_id; retrievals.raw_id
-- references it. Kept 12KB-capped by the loop (_safe_dump).
CREATE TABLE IF NOT EXISTS tool_results_raw (
    raw_id       TEXT PRIMARY KEY,
    ts           REAL NOT NULL,
    session_id   TEXT,
    run_id       TEXT,
    agent_id     TEXT,
    step         INTEGER,
    tool         TEXT,
    args_json    TEXT,
    result_json  TEXT,
    error        TEXT,
    latency_ms   REAL
);
CREATE INDEX IF NOT EXISTS ix_raw_session ON tool_results_raw(session_id);
CREATE INDEX IF NOT EXISTS ix_raw_run     ON tool_results_raw(run_id);
CREATE INDEX IF NOT EXISTS ix_raw_agent   ON tool_results_raw(agent_id);
"""


class SqliteSink:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        # ``check_same_thread=False`` + our own lock: safe for the low volume
        # we produce. ``isolation_level=None`` = autocommit; combined with WAL
        # this gives durable writes without per-call BEGIN/COMMIT bookkeeping.
        self._conn = sqlite3.connect(
            str(self.path), check_same_thread=False, isolation_level=None,
        )
        self._conn.execute("PRAGMA journal_mode=WAL;")
        self._conn.execute("PRAGMA synchronous=NORMAL;")
        self._conn.executescript(_SCHEMA)

    def write(self, record: dict[str, Any]) -> None:
        kind = str(record.get("kind") or "")
        # New structured tables: retrieval_row / tool_result_full events are
        # written to `retrievals` / `tool_results_raw` respectively, NOT the
        # generic `events` payload. They carry benchmark-critical columns
        # (rank / score / snippet / raw_id fk) that would otherwise be
        # trapped inside a JSON blob.
        if kind == "retrieval_row":
            self._write_retrieval_row(record)
            return
        if kind == "tool_result_full":
            self._write_raw_result(record)
            return
        payload = json.dumps(record, ensure_ascii=False, default=str)
        row = (
            float(record.get("ts") or 0.0),
            kind,
            _opt_str(record.get("run_id")),
            _opt_str(record.get("case_id")),
            _opt_str(record.get("suite")),
            _opt_str(record.get("config")),
            payload,
        )
        with self._lock:
            self._conn.execute(
                "INSERT INTO events (ts, kind, run_id, case_id, suite, config, payload) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                row,
            )

    def _write_retrieval_row(self, record: dict[str, Any]) -> None:
        row = (
            float(record.get("ts") or 0.0),
            _opt_str(record.get("session_id")),
            _opt_str(record.get("run_id")),
            _opt_str(record.get("agent_id")),
            _opt_str(record.get("subagent_name")),
            _opt_int(record.get("step")),
            _opt_str(record.get("tool")),
            _opt_str(record.get("doc_id")),
            _opt_int(record.get("rank")),
            _opt_float(record.get("score")),
            _opt_int(record.get("row_count")),
            _opt_str(record.get("snippet")),
            _opt_str(record.get("raw_id")),
        )
        with self._lock:
            self._conn.execute(
                "INSERT INTO retrievals ("
                "ts, session_id, run_id, agent_id, subagent_name, step, tool, "
                "doc_id, rank, score, row_count, snippet, raw_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                row,
            )

    def _write_raw_result(self, record: dict[str, Any]) -> None:
        row = (
            _opt_str(record.get("raw_id")),
            float(record.get("ts") or 0.0),
            _opt_str(record.get("session_id")),
            _opt_str(record.get("run_id")),
            _opt_str(record.get("agent_id")),
            _opt_int(record.get("step")),
            _opt_str(record.get("tool")),
            _opt_str(record.get("args_json")),
            _opt_str(record.get("result_json")),
            _opt_str(record.get("error")),
            _opt_float(record.get("latency_ms")),
        )
        with self._lock:
            # `raw_id` collisions shouldn't happen (uuid8) but a duplicate
            # emit for the same tool call must not crash the loop; INSERT
            # OR REPLACE keeps the most recent snapshot.
            self._conn.execute(
                "INSERT OR REPLACE INTO tool_results_raw ("
                "raw_id, ts, session_id, run_id, agent_id, step, tool, "
                "args_json, result_json, error, latency_ms) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                row,
            )

    def close(self) -> None:
        with self._lock:
            try:
                self._conn.close()
            except Exception:  # pragma: no cover
                pass


def _opt_str(v: Any) -> str | None:
    if v is None:
        return None
    return str(v)


def _opt_int(v: Any) -> int | None:
    if v is None:
        return None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _opt_float(v: Any) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Read side
# ---------------------------------------------------------------------------

def is_sqlite_path(path: str | Path) -> bool:
    """Best-effort suffix sniff: ``.db``, ``.sqlite``, ``.sqlite3``."""
    return Path(path).suffix.lower() in (".db", ".sqlite", ".sqlite3")


def iter_sqlite_events(
    path: str | Path,
    *,
    run_id: str | None = None,
    kind: str | None = None,
    case_id: str | None = None,
    suite: str | None = None,
    config: str | None = None,
) -> Iterator[dict[str, Any]]:
    """Yield records (as dicts) from a SQLite trace DB, ordered by ``ts``.

    Filters compose with AND. Every yielded dict is the full JSON payload
    that was originally passed to :func:`hepuke.observability.trace.emit`.
    """
    filters = {
        "run_id": run_id, "kind": kind, "case_id": case_id,
        "suite": suite, "config": config,
    }
    clauses: list[str] = []
    params: list[Any] = []
    for col, val in filters.items():
        if val is not None:
            clauses.append(f"{col} = ?")
            params.append(val)
    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""

    conn = sqlite3.connect(str(path))
    try:
        cur = conn.execute(
            f"SELECT payload FROM events {where} ORDER BY ts",
            params,
        )
        for (payload,) in cur:
            try:
                yield json.loads(payload)
            except json.JSONDecodeError:
                continue
    finally:
        conn.close()


def known_run_ids(path: str | Path) -> list[str]:
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.execute(
            "SELECT run_id, COUNT(*) FROM events "
            "WHERE run_id IS NOT NULL GROUP BY run_id ORDER BY MAX(ts) DESC"
        )
        return [row[0] for row in cur]
    finally:
        conn.close()


# retrievals / tool_results_raw readers ---------------------------------------

_RETRIEVAL_COLS = (
    "id", "ts", "session_id", "run_id", "agent_id", "subagent_name",
    "step", "tool", "doc_id", "rank", "score", "row_count", "snippet",
    "raw_id",
)

_RAW_COLS = (
    "raw_id", "ts", "session_id", "run_id", "agent_id", "step", "tool",
    "args_json", "result_json", "error", "latency_ms",
)


def iter_retrieval_rows(
    path: str | Path,
    *,
    session_id: str | None = None,
    run_id: str | None = None,
    agent_id: str | None = None,
    doc_id: str | None = None,
    tool: str | None = None,
) -> Iterator[dict[str, Any]]:
    """Yield structured retrieval rows from a SQLite trace DB, ordered by ``ts``.

    Every column becomes a dict key; NULLs come through as ``None``. All
    filters compose with AND.
    """
    filters = {
        "session_id": session_id, "run_id": run_id, "agent_id": agent_id,
        "doc_id": doc_id, "tool": tool,
    }
    clauses: list[str] = []
    params: list[Any] = []
    for col, val in filters.items():
        if val is not None:
            clauses.append(f"{col} = ?")
            params.append(val)
    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
    cols = ", ".join(_RETRIEVAL_COLS)
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.execute(
            f"SELECT {cols} FROM retrievals {where} ORDER BY ts, id",
            params,
        )
        for row in cur:
            yield dict(zip(_RETRIEVAL_COLS, row))
    finally:
        conn.close()


def iter_tool_results_raw(
    path: str | Path,
    *,
    session_id: str | None = None,
    run_id: str | None = None,
    agent_id: str | None = None,
    raw_id: str | None = None,
    tool: str | None = None,
) -> Iterator[dict[str, Any]]:
    """Yield raw tool-result rows, ordered by ``ts``."""
    filters = {
        "session_id": session_id, "run_id": run_id, "agent_id": agent_id,
        "raw_id": raw_id, "tool": tool,
    }
    clauses: list[str] = []
    params: list[Any] = []
    for col, val in filters.items():
        if val is not None:
            clauses.append(f"{col} = ?")
            params.append(val)
    where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
    cols = ", ".join(_RAW_COLS)
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.execute(
            f"SELECT {cols} FROM tool_results_raw {where} ORDER BY ts, raw_id",
            params,
        )
        for row in cur:
            yield dict(zip(_RAW_COLS, row))
    finally:
        conn.close()
