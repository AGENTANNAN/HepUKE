"""In-process session store for the /agent/chat endpoint.

Holds paused RagAgent runs between HTTP turns. Deliberately in-memory —
this is Phase 1 of the interactive path; a disk/redis backend can drop in
later behind the same interface. Threadsafe (single big lock is fine, the
critical section is a dict op).

Bounded by:
  * `max_sessions` — LRU-evict the least recently touched entry.
  * `ttl_seconds` — reap on read.

The store never persists across process restarts; the front-end must be
prepared to see a session_id go missing (return 404 / start over).
"""
from __future__ import annotations

import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class _Entry:
    session: Any  # AgentSession (avoid importing to keep this module light)
    agent: Any   # RagAgent (same reason)
    created_at: float
    last_used: float
    turn_count: int = 0


class SessionStore:
    def __init__(self, *, max_sessions: int = 256, ttl_seconds: float = 3600.0) -> None:
        self._data: dict[str, _Entry] = {}
        self._lock = threading.Lock()
        self.max_sessions = int(max_sessions)
        self.ttl_seconds = float(ttl_seconds)

    def create(self, agent: Any, session: Any) -> str:
        """Register a paused session, return the new id.

        If the session already carries a non-empty ``session_id`` (set by
        the caller before running the loop, so ALL trace rows can be
        stamped with the same id across turns), reuse it; otherwise mint a
        fresh one. Duplicate ids overwrite silently — the previous entry
        for the same id is dropped.
        """
        sid = getattr(session, "session_id", None) or uuid.uuid4().hex
        now = time.time()
        with self._lock:
            self._reap_locked(now)
            self._evict_if_full_locked()
            self._data[sid] = _Entry(
                session=session, agent=agent,
                created_at=now, last_used=now, turn_count=1,
            )
        return sid

    def get(self, sid: str) -> Optional[_Entry]:
        """Return the entry and refresh its last-used stamp."""
        now = time.time()
        with self._lock:
            entry = self._data.get(sid)
            if entry is None:
                return None
            # TTL reap check — treat expired as gone.
            if now - entry.last_used > self.ttl_seconds:
                self._data.pop(sid, None)
                return None
            entry.last_used = now
            entry.turn_count += 1
            return entry

    def drop(self, sid: str) -> None:
        with self._lock:
            self._data.pop(sid, None)

    def size(self) -> int:
        with self._lock:
            return len(self._data)

    # --- internals ------------------------------------------------------

    def _reap_locked(self, now: float) -> None:
        expired = [
            sid for sid, e in self._data.items()
            if now - e.last_used > self.ttl_seconds
        ]
        for sid in expired:
            self._data.pop(sid, None)

    def _evict_if_full_locked(self) -> None:
        if len(self._data) < self.max_sessions:
            return
        # LRU by last_used.
        oldest = min(self._data.items(), key=lambda kv: kv[1].last_used)
        self._data.pop(oldest[0], None)
