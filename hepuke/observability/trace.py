"""Trace context + pluggable-sink emitter.

Interface (all frozen):
    configure_trace(*, enabled=None, file=None, include_context=None,
                    backend=None, db_path=None)
    set_trace_context(**kw)
    clear_trace_context()
    emit(kind, **fields)

Configuration comes from three sources, in precedence order (later wins):
  1. ``TraceConfig`` (config.yaml) → wired by ``HepUKE.__init__``
  2. HEPUKE__TRACE__* environment variables (see ``_bootstrap_from_env``)
  3. explicit ``configure_trace(...)`` calls (e.g. ``benchmarks/runner.py``)

Backends
--------
``backend='jsonl'``  → append-JSONL at ``file`` (legacy).
``backend='sqlite'`` → SQLite DB at ``db_path`` (default).
``backend='both'``   → both sinks fire in order.

Backwards compat: callers that only pass ``file=`` still work; the module
defaults to ``sqlite`` at ``TraceConfig.db_path`` unless ``backend`` is set.
"""
from __future__ import annotations

import json  # noqa: F401  (kept for downstream import compat)
import os
import threading
import time
from pathlib import Path
from typing import Any

from hepuke.observability.sinks import ConsoleSink, JsonlSink, SqliteSink, TraceSink


_STATE = threading.local()
_WRITE_LOCK = threading.Lock()

_CONFIG: dict[str, Any] = {
    "enabled": False,
    "backend": "sqlite",
    "file": "./trace/rag_trace.jsonl",
    "db_path": "./storage/traces.db",
    "include_context": True,
    "console": False,  # attach ConsoleSink alongside the primary backend
}
_SINKS: list[TraceSink] = []
_SINK_KEY: tuple | None = None  # (backend, file, db_path, console) for change detection


# ---------------------------------------------------------------------------
# configure
# ---------------------------------------------------------------------------

def configure_trace(
    *,
    enabled: bool | None = None,
    file: str | None = None,
    include_context: bool | None = None,
    backend: str | None = None,
    db_path: str | None = None,
    console: bool | None = None,
) -> None:
    """Update config; rebuild sinks lazily on next :func:`emit`."""
    global _SINK_KEY
    if enabled is not None:
        _CONFIG["enabled"] = bool(enabled)
    if file is not None:
        _CONFIG["file"] = str(file)
    if include_context is not None:
        _CONFIG["include_context"] = bool(include_context)
    if backend is not None:
        if backend not in ("jsonl", "sqlite", "both"):
            raise ValueError(f"unknown trace backend: {backend!r}")
        _CONFIG["backend"] = backend
    if db_path is not None:
        _CONFIG["db_path"] = str(db_path)
    if console is not None:
        _CONFIG["console"] = bool(console)
    # Invalidate current sinks; ``_ensure_sinks`` will rebuild on next emit.
    _SINK_KEY = None
    _close_sinks()


def _ensure_sinks() -> list[TraceSink]:
    global _SINK_KEY, _SINKS
    key = (
        _CONFIG["backend"], _CONFIG["file"], _CONFIG["db_path"],
        _CONFIG["console"],
    )
    if _SINK_KEY == key and _SINKS:
        return _SINKS
    _close_sinks()
    sinks: list[TraceSink] = []
    if _CONFIG["backend"] in ("jsonl", "both"):
        sinks.append(JsonlSink(_CONFIG["file"]))
    if _CONFIG["backend"] in ("sqlite", "both"):
        sinks.append(SqliteSink(_CONFIG["db_path"]))
    if _CONFIG["console"]:
        sinks.append(ConsoleSink())
    _SINKS = sinks
    _SINK_KEY = key
    return _SINKS


def _close_sinks() -> None:
    global _SINKS
    for s in _SINKS:
        try:
            s.close()
        except Exception:  # pragma: no cover
            pass
    _SINKS = []


# ---------------------------------------------------------------------------
# context
# ---------------------------------------------------------------------------

def _ctx() -> dict[str, Any]:
    return getattr(_STATE, "ctx", {})


def set_trace_context(**kwargs: Any) -> None:
    ctx = dict(_ctx())
    ctx.update(kwargs)
    _STATE.ctx = ctx


def clear_trace_context() -> None:
    _STATE.ctx = {}


# ---------------------------------------------------------------------------
# emit
# ---------------------------------------------------------------------------

def emit(kind: str, **fields: Any) -> None:
    """Append a single trace record to all configured sinks. No-op when off."""
    if not _CONFIG["enabled"]:
        return
    record: dict[str, Any] = {"ts": time.time(), "kind": kind}
    if _CONFIG["include_context"]:
        record.update(_ctx())
    record.update(fields)
    with _WRITE_LOCK:
        sinks = _ensure_sinks()
        for s in sinks:
            try:
                s.write(record)
            except Exception:  # pragma: no cover
                # Never let a trace failure kill the caller.
                continue


# ---------------------------------------------------------------------------
# introspection (small helpers for CLI / tests)
# ---------------------------------------------------------------------------

def current_config() -> dict[str, Any]:
    """Return a snapshot of the active trace config."""
    return dict(_CONFIG)


# ---------------------------------------------------------------------------
# env bootstrap
# ---------------------------------------------------------------------------

def _bootstrap_from_env() -> None:
    """Pick up TraceConfig-style env at import time.

    Env keys (all optional):
      HEPUKE__TRACE__ENABLED   ("1"/"true"/"yes" → on)
      HEPUKE__TRACE__BACKEND   ("jsonl"/"sqlite"/"both")
      HEPUKE__TRACE__FILE      (path for JsonlSink)
      HEPUKE__TRACE__DB_PATH   (path for SqliteSink)
      HEPUKE__TRACE__CONSOLE   ("1"/"true"/"yes" → also print to stderr)
    """
    enabled = os.environ.get("HEPUKE__TRACE__ENABLED", "").lower() in ("1", "true", "yes")
    console = os.environ.get("HEPUKE__TRACE__CONSOLE", "").lower() in ("1", "true", "yes")
    if not (enabled or console):
        return
    kwargs: dict[str, Any] = {"enabled": True}
    backend = os.environ.get("HEPUKE__TRACE__BACKEND")
    if backend:
        kwargs["backend"] = backend.lower()
    f = os.environ.get("HEPUKE__TRACE__FILE")
    if f:
        kwargs["file"] = f
    db = os.environ.get("HEPUKE__TRACE__DB_PATH")
    if db:
        kwargs["db_path"] = db
    if console:
        kwargs["console"] = True
    configure_trace(**kwargs)


_bootstrap_from_env()
