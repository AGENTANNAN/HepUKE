"""Memo / document CLI-style tools exposed to the agent loop.

These wrap plain-old shell operations (grep, ls, cat) into JSON-returning
functions the agent can call. Design tenets:
  * `grep` uses `subprocess.run([...], shell=False)` — no shell metachars
    ever reach a shell interpreter. The `pattern` is passed as an argv element
    directly to `grep -E`.
  * `path_root` is resolved and confined to a caller-configured base; escape
    attempts (`..`) are rejected before any I/O.
  * `read_memo` prefers Milvus (`HybridRetriever.get_full_doc`) so section
    ordering is honored; falls back to `{path_root}/{doc_id}.md` on disk.

SANDBOX POLICY — do NOT broaden this without a real sandbox.
  The `grep` tool below is INTENTIONALLY the only "shell" surface the agent
  gets. `argv[0]` is hard-coded to "grep", `shell=False`, and the corpus root
  is confined by `_safe_root`/`_confine`. So the LLM CANNOT run `rm`, `cat`,
  `ls`, pipes (`|`), redirects (`>`), `;`, `&&`, `$(...)`; those characters are
  bytes passed to grep's regex engine, never expanded by `/bin/sh`.
  If a future task needs a wider shell (e.g. `find` + `wc`), do NOT flip
  `shell=True`. Add a separate `run_cmd` tool with:
    * `shell=False`
    * argv[0] whitelist (e.g. `{"grep","find","wc","head","tail"}`)
    * argv-only quoting, per-arg validation
    * timeout + result-size caps (already present here)
    * *and* a real OS-level sandbox (Firejail / Bubblewrap / a docker run with
      `--network=none --read-only`) if the corpus root is not already
      immutable. Anything that opens `rm`/mutation MUST be sandboxed.
"""
from __future__ import annotations

import fnmatch
import os
import subprocess
from pathlib import Path
from typing import Any, Iterable, Optional

from hepuke.observability.trace import emit


class MemoToolError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# path safety
# ---------------------------------------------------------------------------

def safe_root(path_root: str | os.PathLike) -> Path:
    root = Path(path_root).expanduser().resolve()
    if not root.exists():
        raise MemoToolError(f"path_root does not exist: {root}")
    if not root.is_dir():
        raise MemoToolError(f"path_root is not a directory: {root}")
    return root


def confine(root: Path, candidate: Path) -> Path:
    resolved = candidate.expanduser().resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        raise MemoToolError(
            f"path {resolved} escapes path_root {root}"
        )
    return resolved


# Legacy aliases (used across the package before the rename).
_safe_root = safe_root
_confine = confine


# ---------------------------------------------------------------------------
# grep
# ---------------------------------------------------------------------------

def grep(
    pattern: str,
    *,
    path_root: str,
    doc_ids: Optional[Iterable[str]] = None,
    include: str = "*.md",
    context: int = 3,
    max_results: int = 200,
    ignore_case: bool = False,
    timeout: int = 30,
) -> list[dict]:
    """Run `grep -rnE --include={include} -C{context}` under `path_root`.

    Parameters
    ----------
    pattern : ERE regular expression (passed as argv, never a shell string).
    doc_ids : if given, only search `{path_root}/{doc_id}.md` for each id.
    ignore_case : add `-i`.

    Returns
    -------
    List of `{'file': relpath, 'line': int, 'text': str, 'is_match': bool}`.
    """
    root = _safe_root(path_root)

    argv: list[str] = ["grep", "-HrnE", f"--include={include}", f"-C{context}"]
    if ignore_case:
        argv.append("-i")
    argv.append("--")
    argv.append(pattern)

    if doc_ids:
        targets: list[str] = []
        for did in doc_ids:
            # Try direct file first; fall back to a filename glob under root.
            direct = _confine(root, root / f"{did}.md")
            if direct.is_file():
                targets.append(str(direct))
                continue
            hits = list(root.rglob(f"{did}*.md"))
            for h in hits:
                _confine(root, h)
                targets.append(str(h))
        if not targets:
            emit("tool_call", name="grep", pattern=pattern, n_results=0, note="no doc_id match")
            return []
        argv.extend(targets)
    else:
        argv.append(str(root))

    try:
        proc = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        raise MemoToolError(f"grep timed out after {timeout}s")

    if proc.returncode not in (0, 1):
        # 0=match, 1=no match, >=2=error
        raise MemoToolError(
            f"grep failed (rc={proc.returncode}): {proc.stderr.strip()[:200]}"
        )

    out: list[dict] = []
    for raw in proc.stdout.splitlines():
        if not raw or raw == "--":
            continue
        # `file:line:text` for match lines, `file-line-text` for context lines.
        m_match = _split_grep_line(raw, sep=":")
        m_ctx = _split_grep_line(raw, sep="-")
        if m_match:
            file, line, text = m_match
            is_match = True
        elif m_ctx:
            file, line, text = m_ctx
            is_match = False
        else:
            continue
        try:
            rel = str(Path(file).resolve().relative_to(root))
        except (ValueError, OSError):
            rel = file
        out.append({"file": rel, "line": line, "text": text, "is_match": is_match})
        if len(out) >= max_results:
            break

    emit(
        "tool_call",
        name="grep",
        pattern=pattern,
        n_results=len(out),
        n_matches=sum(1 for h in out if h["is_match"]),
    )
    return out


def _split_grep_line(line: str, *, sep: str) -> Optional[tuple[str, int, str]]:
    # Expect `path{sep}lineno{sep}text` with lineno an int. The path itself may
    # contain colons; iterate from the left and require the second field to be
    # a bare integer.
    first = line.find(sep)
    if first == -1:
        return None
    second = line.find(sep, first + 1)
    if second == -1:
        return None
    file = line[:first]
    lineno_str = line[first + 1:second]
    text = line[second + 1:]
    try:
        lineno = int(lineno_str)
    except ValueError:
        return None
    return file, lineno, text


# ---------------------------------------------------------------------------
# list_memos
# ---------------------------------------------------------------------------

def list_memos(
    *,
    path_root: str,
    pattern: str = "*.md",
    limit: int = 5000,
) -> list[str]:
    """Return sorted relative paths of files matching `pattern` under `path_root`.

    `pattern` is a shell fnmatch glob (`BAM-*.md`, `*.md`, `arxiv/*.md`).
    """
    root = _safe_root(path_root)
    hits: list[str] = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root)
        if fnmatch.fnmatch(str(rel), pattern) or fnmatch.fnmatch(p.name, pattern):
            hits.append(str(rel))
            if len(hits) >= limit:
                break
    emit("tool_call", name="list_memos", pattern=pattern, n_results=len(hits))
    return hits


# ---------------------------------------------------------------------------
# read_memo / get_section
# ---------------------------------------------------------------------------

def read_memo(
    doc_id: str,
    *,
    retriever: Any = None,
    collection: Optional[str] = None,
    path_root: Optional[str] = None,
    section_order: Optional[list[str]] = None,
    max_chars: int = 200_000,
) -> dict:
    """Fetch a full document.

    Precedence:
      1. If `retriever` + `collection` given, call `retriever.get_full_doc()`
         and return `{'source': 'milvus', 'doc_id', 'title', 'authors',
         'date', 'sections': [...], 'text': concatenated}`.
      2. Else read `{path_root}/{doc_id}.md` from disk and return
         `{'source': 'disk', 'doc_id', 'path', 'text'}`.
    """
    if retriever is not None and collection:
        try:
            got = retriever.get_full_doc(
                doc_id, collection=collection, section_order=section_order,
            )
        except Exception as e:
            emit("tool_call", name="read_memo", doc_id=doc_id, error=str(e))
            got = None
        if got:
            body = "\n\n".join(
                f"## {s.get('section_title') or s.get('section_type', '')}"
                f"\n{s.get('content', '')}"
                for s in got.get("sections", [])
            )
            emit(
                "tool_call", name="read_memo", doc_id=doc_id,
                source="milvus", n_sections=len(got.get("sections", [])),
            )
            return {"source": "milvus", **got, "text": body[:max_chars]}

    if path_root:
        root = _safe_root(path_root)
        candidate = _confine(root, root / f"{doc_id}.md")
        if candidate.is_file():
            text = candidate.read_text(encoding="utf-8")[:max_chars]
            emit(
                "tool_call", name="read_memo", doc_id=doc_id, source="disk",
                bytes=len(text),
            )
            return {
                "source": "disk", "doc_id": doc_id,
                "path": str(candidate.relative_to(root)),
                "text": text,
            }
        # Try a loose glob fallback.
        matches = list(root.rglob(f"{doc_id}*.md"))
        if matches:
            _confine(root, matches[0])
            text = matches[0].read_text(encoding="utf-8")[:max_chars]
            emit(
                "tool_call", name="read_memo", doc_id=doc_id,
                source="disk_glob", bytes=len(text),
            )
            return {
                "source": "disk", "doc_id": doc_id,
                "path": str(matches[0].relative_to(root)),
                "text": text,
            }

    emit("tool_call", name="read_memo", doc_id=doc_id, error="not_found")
    return {"source": "none", "doc_id": doc_id, "text": "", "error": "doc not found"}


def get_section(
    doc_id: str,
    section: str,
    *,
    retriever: Any,
    collection: str,
    max_chars: int = 60_000,
) -> dict:
    """Fetch parent chunks matching `(doc_id, section_type)` and return joined text.

    Requires a chunked schema (`chunk_type` + `section_type` columns). Uses
    `retriever.doc_id_field` so it works with `bam_id` / `paper_id` / `id`
    primary keys — not just the P0 `doc_id`.
    """
    schema_field_names = _schema_field_names(retriever.store, collection)
    if "chunk_type" not in schema_field_names or "section_type" not in schema_field_names:
        return {
            "doc_id": doc_id,
            "section": section,
            "n_rows": 0,
            "text": "",
            "error": (
                f"collection `{collection}` has no chunk_type/section_type "
                "columns; get_section is not applicable here."
            ),
        }

    retriever.store.load_collection(collection)
    doc_id_field = getattr(retriever, "doc_id_field", "doc_id")
    rows = retriever.store.client.query(
        collection_name=collection,
        filter=(
            f'chunk_type == "parent" and '
            f'{doc_id_field} == "{doc_id}" and '
            f'section_type == "{section}"'
        ),
        output_fields=["section_number", "section_title", "content"],
        limit=50,
    )
    text = "\n\n".join(r.get("content", "") for r in rows)[:max_chars]
    emit(
        "tool_call", name="get_section",
        doc_id=doc_id, section=section, n_rows=len(rows),
    )
    return {
        "doc_id": doc_id, "section": section,
        "n_rows": len(rows), "text": text,
    }


def _schema_field_names(store, collection: str) -> set[str]:
    """Cheap schema probe used by get_section to avoid crashing on plain
    metadata collections (`hypernews_hybrid` and friends)."""
    try:
        desc = store.describe_collection(collection)
    except Exception:
        return set()
    return {f["name"] for f in desc.get("fields", [])}
