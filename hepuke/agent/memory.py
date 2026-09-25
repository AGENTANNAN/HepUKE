"""Per-session file-backed memory for the RAG agent.

Minimal starter (see [[design_rag_memory]]): two file types only —
`user_intent.md` and `dead_end_<hash>.md`. Everything is written eagerly
under `<session_root>/memory/`, so a subsequent turn (or a subsequent
process, once SessionStore gets a disk backend) can pick up context by
re-reading the directory.

Why files not a dict:
  * Human-inspectable — the user can `cat` a session's memory to see
    exactly what the agent is carrying between turns.
  * Trivially persistable — no schema migration when session storage
    goes to disk / redis / anything else.
  * Cheap enough at this scale (a session accumulates dozens of files,
    not thousands).

Deliberately NOT included in this starter:
  * `finding_*` — needs evidence-distillation logic. Wait until the
    trace-fallback evidence contract has real-machine mileage.
  * Semantic retrieval — a small session doesn't need it; a full-corpus
    embed pass just to filter a handful of files is overkill.
  * Eviction / TTL / global promote — one session lives inside its
    SessionStore TTL already; that's enough for the minimal starter.
"""
from __future__ import annotations

import hashlib
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


_DEAD_END_FILE_RE = re.compile(r"^dead_end_[0-9a-f]{8}\.md$")


def _hash_question(subagent: str, question: str) -> str:
    """Stable short hash for a (subagent, normalised question) pair.

    Whitespace + case normalised so trivial re-wordings collide and we
    don't record the same dead-end twice. 8 hex chars is plenty at the
    per-session cardinality we expect (dozens, not millions).
    """
    norm = re.sub(r"\s+", " ", (question or "").strip().lower())
    key = f"{subagent}||{norm}".encode("utf-8")
    return hashlib.sha1(key).hexdigest()[:8]


@dataclass
class DeadEndEntry:
    subagent: str
    question: str
    reason: str
    steps: int
    created_at: float
    path: Path


class SessionMemory:
    """File-backed memory for one agent session.

    Not threadsafe against concurrent writes to the SAME session — the
    fan-out inside a single run() writes to different files (one per
    subagent × question), so collisions are rare, and the write path
    is `open('w')` which is atomic on POSIX for small files. Cross-
    session isolation is by directory.
    """

    def __init__(self, session_root: str | os.PathLike[str]) -> None:
        self.root = Path(session_root)
        self.memory_dir = self.root / "memory"
        self.memory_dir.mkdir(parents=True, exist_ok=True)
        self._user_intent_path = self.memory_dir / "user_intent.md"

    # ------------------------------------------------------------------
    # user_intent
    # ------------------------------------------------------------------

    def write_user_intent(self, text: str) -> None:
        """Overwrite `user_intent.md`. First non-empty question of a
        session is the intent; subsequent overwrites are only expected
        when the user reformulates via ask_user resume."""
        payload = (text or "").strip()
        if not payload:
            return
        body = (
            "---\n"
            f"created_at: {_iso_now()}\n"
            "type: user_intent\n"
            "---\n\n"
            f"{payload}\n"
        )
        self._user_intent_path.write_text(body, encoding="utf-8")

    def read_user_intent(self) -> Optional[str]:
        if not self._user_intent_path.exists():
            return None
        return self._user_intent_path.read_text(encoding="utf-8")

    # ------------------------------------------------------------------
    # dead_end
    # ------------------------------------------------------------------

    def record_dead_end(
        self,
        subagent: str,
        question: str,
        reason: str,
        *,
        steps: int = 0,
    ) -> Path:
        """Idempotent by (subagent, normalised question).

        Re-recording the same dead-end just refreshes the file — a
        subagent that tried this exact question again and failed
        again is still one dead-end, not two. The `reason` field is
        updated to the latest attempt's reason so we always keep the
        most recent explanation.
        """
        h = _hash_question(subagent, question)
        path = self.memory_dir / f"dead_end_{h}.md"
        body = (
            "---\n"
            f"created_at: {_iso_now()}\n"
            "type: dead_end\n"
            f"subagent: {subagent}\n"
            f"steps: {int(steps)}\n"
            "---\n\n"
            f"**Question**: {question}\n\n"
            f"**Reason**: {reason}\n"
        )
        path.write_text(body, encoding="utf-8")
        return path

    def has_dead_end(self, subagent: str, question: str) -> bool:
        """True if the same (subagent, question) already recorded dead.
        Callers can use this to short-circuit a re-dispatch."""
        h = _hash_question(subagent, question)
        return (self.memory_dir / f"dead_end_{h}.md").exists()

    def list_dead_ends(self) -> list[DeadEndEntry]:
        """Enumerate all dead-end files, freshest-first (by mtime)."""
        entries: list[DeadEndEntry] = []
        for p in sorted(
            self.memory_dir.glob("dead_end_*.md"),
            key=lambda x: x.stat().st_mtime,
            reverse=True,
        ):
            if not _DEAD_END_FILE_RE.match(p.name):
                continue
            entries.append(_parse_dead_end(p))
        return entries

    # ------------------------------------------------------------------
    # context injection
    # ------------------------------------------------------------------

    def load_context(self, *, max_dead_ends: int = 8) -> str:
        """Return a system-prompt-ready block summarising session memory.

        Empty string when there is nothing to inject — callers should
        prepend directly without a separator check.
        """
        parts: list[str] = []
        intent = self._read_user_intent_body()
        if intent:
            parts.append("## Session intent\n" + intent)

        dead_ends = self.list_dead_ends()[:max_dead_ends]
        if dead_ends:
            lines = ["## Known dead-ends (do NOT re-dispatch these)"]
            for e in dead_ends:
                lines.append(
                    f"- `{e.subagent}` on: {e.question!r} — {e.reason}"
                )
            parts.append("\n".join(lines))

        if not parts:
            return ""
        return "\n\n".join(parts) + "\n"

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------

    def _read_user_intent_body(self) -> Optional[str]:
        raw = self.read_user_intent()
        if not raw:
            return None
        # Strip frontmatter block if present.
        if raw.startswith("---"):
            end = raw.find("\n---", 3)
            if end != -1:
                raw = raw[end + 4:]
        return raw.strip() or None


def _iso_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _parse_dead_end(path: Path) -> DeadEndEntry:
    text = path.read_text(encoding="utf-8")
    subagent = _grep_frontmatter(text, "subagent") or "?"
    steps_raw = _grep_frontmatter(text, "steps") or "0"
    try:
        steps = int(steps_raw)
    except ValueError:
        steps = 0
    created_raw = _grep_frontmatter(text, "created_at") or ""
    # We don't need to parse the timestamp back; the mtime is the sort key.
    created_at = path.stat().st_mtime
    del created_raw
    question = _grep_body_field(text, "Question") or ""
    reason = _grep_body_field(text, "Reason") or ""
    return DeadEndEntry(
        subagent=subagent,
        question=question,
        reason=reason,
        steps=steps,
        created_at=created_at,
        path=path,
    )


def _grep_frontmatter(text: str, key: str) -> Optional[str]:
    m = re.search(rf"^{re.escape(key)}:\s*(.+)$", text, re.MULTILINE)
    return m.group(1).strip() if m else None


def _grep_body_field(text: str, label: str) -> Optional[str]:
    m = re.search(rf"\*\*{re.escape(label)}\*\*:\s*(.+)$", text, re.MULTILINE)
    return m.group(1).strip() if m else None


__all__ = ["SessionMemory", "DeadEndEntry"]
