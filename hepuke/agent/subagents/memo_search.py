"""`memo_search` subagent — specialist over the on-disk markdown corpus.

Owns the memo tree named by `wiki.memo_md_root` plus two navigation
aids one level up: `INDEX.md` (single-table of every doc's id/title/abstract
lede) and `summaries/<id>.md` (~1KB abstract + heading map). The subagent
prefers those over body-text grep for any discovery-type question.

Deliberately NOT wired to Milvus or SQLite — cross-corpus joins happen at
the parent level. If the parent needs metadata about a doc, it calls the
`milvus_search` subagent; if it needs full body text, it calls this one.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from hepuke.agent.subagents.base import make_subagent_tool
from hepuke.agent.tools import Tool
from hepuke.tools import memo_tools
from hepuke.tools import shell_tools


MEMO_SEARCH_SYSTEM_PROMPT = """\
You are the memo-search specialist over an on-disk markdown corpus of BESIII
memos (filenames like `BAM-00713_Memo_SP.md`). You have SEVEN tools; nothing
else exists.

CORPUS LAYOUT — three tiers, cheapest first:
  1. `INDEX.md`             one table row per doc (1150 rows), with
                             id, title, abstract lede, bytes. THIS IS THE
                             TABLE OF CONTENTS OF THE WHOLE CORPUS. Access
                             it via `search_index`.
  2. `summaries/<id>.md`    ~1 KB per doc: abstract + heading map.
                             Access it via `read_summary(doc_id)`.
  3. `markdowns/<id>*.md`   full body ~50 KB per doc. Access via
                             `preview_memo` (cheap) / `read_memo` (expensive).

Routing rules — ALWAYS start at the cheapest tier that could answer:
  * DISCOVERY question ("哪些论文关于 X" / "list papers about X" /
    "how many BAMs mention X")
        → **MUST** call `search_index` FIRST. INDEX.md contains every doc's
          title and abstract lede in one file — if a topic exists in the
          corpus, it shows up here. Do NOT grep body files for discovery.
  * KNOWN doc_id, metadata-only ("who wrote BAM-00713" / "abstract of X")
        → `read_summary(doc_id)` — 1 KB, usually enough.
  * KNOWN doc_id, needs body content (methodology, systematic errors,
    specific numbers) → `preview_memo(doc_id)` first to confirm relevance,
          then `read_memo` or `grep` only if the preview looks promising.
  * NEEDS an exact token / LaTeX / number inside body text and you already
    have a candidate doc_id set → `grep(pattern, doc_ids=[...])`. Do NOT
    grep the whole corpus without doc_ids — it is slow and dilutes hits.
  * COUNTS / filename globs → `list_memos` / `run_cmd`.

Tool inventory:
  * `search_index`  — grep INDEX.md by keyword. Returns [{doc_id,title,
                       abstract_lede,bytes}]. **First hop for discovery.**
  * `read_summary`  — read summaries/<doc_id>.md (~1 KB abstract + heads).
  * `preview_memo`  — first ~500 chars + headings of a body file (~1 KB).
  * `read_memo`     — full body pull (~50 KB). EXPENSIVE.
  * `grep`          — POSIX ERE regex over body files. Structured
                       {file,line,text}. Narrow with `doc_ids` when you can.
  * `list_memos`    — filename glob (`BAM-*.md`), for counting/discovery.
  * `run_cmd`       — wc / head / tail / find / sort / uniq / cut / tr /
                       cat / ls / grep. Read-only, no pipes.

Cost / stop conditions:
  * `read_memo` is a last resort. Never call it speculatively.
  * Budget: at most 5 tool calls per question. If not done after 5,
    return `dry=true`.

Contract — your FINAL message MUST be a JSON object, no prose, no fences:
  {
    "answer":   "<one-paragraph natural-language answer, or null>",
    "evidence": [ {"source": "memo", "doc_id": "...", "path": "...",
                   "snippet": "..."}, ... ],
    "dry":      true|false
  }
`dry=true` means the corpus had nothing relevant — set this when every
tier you tried returned zero rows. The parent decides whether to escalate.
Never invent doc_ids or file paths.

Grounding: every claim in `answer` must trace to a `doc_id` in `evidence`.
If you cannot ground it, drop the claim.
"""


def _grep_tool(memo_root: str) -> Tool:
    def _call(
        pattern: str,
        doc_ids: Optional[list] = None,
        context: int = 3,
        ignore_case: bool = False,
        max_results: int = 50,
    ) -> list:
        return memo_tools.grep(
            pattern=pattern,
            path_root=memo_root,
            doc_ids=doc_ids,
            context=context,
            ignore_case=ignore_case,
            max_results=max_results,
        )

    parameters = {
        "type": "object",
        "properties": {
            "pattern": {"type": "string",
                        "description": "POSIX ERE regex; passed as argv"},
            "doc_ids": {"type": "array", "items": {"type": "string"},
                        "description": "narrow to `{memo_root}/{doc_id}*.md`"},
            "context": {"type": "integer", "default": 3, "minimum": 0, "maximum": 20},
            "ignore_case": {"type": "boolean", "default": False},
            "max_results": {"type": "integer", "default": 50, "minimum": 1, "maximum": 200},
        },
        "required": ["pattern"],
    }
    return Tool(
        name="grep",
        description=(
            "POSIX ERE regex over the memo markdown corpus. Structured "
            "{file,line,text} hits. Prefer over `run_cmd grep`."
        ),
        parameters=parameters,
        call=_call,
    )


def _list_memos_tool(memo_root: str) -> Tool:
    def _call(pattern: str = "*.md") -> dict:
        files = memo_tools.list_memos(path_root=memo_root, pattern=pattern)
        return {"count": len(files), "pattern": pattern, "files": files}

    parameters = {
        "type": "object",
        "properties": {
            "pattern": {"type": "string", "default": "*.md",
                        "description": "fnmatch glob, e.g. `BAM-007*.md`"},
        },
    }
    return Tool(
        name="list_memos",
        description=(
            "List markdown filenames under the memo root (fnmatch glob). "
            "Cheap — use for counting / discovery."
        ),
        parameters=parameters,
        call=_call,
    )


def _preview_memo_tool(memo_root: str, *, head_chars: int = 500) -> Tool:
    """Cheap preview: title + first N chars + first-level headings.

    Meant as the FIRST hop when the subagent has a `doc_id` and wants to
    decide whether it needs the full body. Returns ~1 KB instead of the
    ~50 KB `read_memo` costs. Only if the preview says the doc is
    relevant should the subagent follow up with `read_memo`.
    """
    import re as _re

    _HEADING_RE = _re.compile(r"^(#{1,3})\s+(.+?)\s*$", _re.MULTILINE)

    def _call(doc_id: str) -> dict:
        candidates = [doc_id]
        if doc_id.endswith(".md"):
            candidates.append(doc_id[:-3])
        m = _re.match(r"^(BAM-\d+)", doc_id)
        if m and m.group(1) not in candidates:
            candidates.append(m.group(1))

        for cand in candidates:
            got = memo_tools.read_memo(
                doc_id=cand, retriever=None, collection=None,
                path_root=memo_root,
            )
            if got.get("source") == "disk" and got.get("text"):
                text = got["text"]
                title = ""
                for line in text.splitlines():
                    stripped = line.strip()
                    if stripped.startswith("#"):
                        title = stripped.lstrip("# ").strip()
                        break
                headings = [
                    f"{'#' * len(m.group(1))} {m.group(2).strip()}"
                    for m in _HEADING_RE.finditer(text)
                ][:30]
                return {
                    "source": "disk",
                    "doc_id": got.get("doc_id", cand),
                    "path": got.get("path"),
                    "title": title,
                    "head": text[:head_chars],
                    "headings": headings,
                    "total_chars": len(text),
                    "hint": (
                        "Call `read_memo` on this doc_id only if the head "
                        "or a heading looks directly relevant to the "
                        "question. Otherwise stop and answer from the "
                        "preview."
                    ),
                }
        return {
            "source": "none", "doc_id": doc_id,
            "error": "doc not found on disk",
        }

    parameters = {
        "type": "object",
        "properties": {
            "doc_id": {
                "type": "string",
                "description": "e.g. 'BAM-00713' — filename suffix tolerated",
            },
        },
        "required": ["doc_id"],
    }
    return Tool(
        name="preview_memo",
        description=(
            "Cheap preview of one memo: title + first ~500 chars + all "
            "headings (~1 KB total). ALWAYS call this before `read_memo` "
            "so you can decide whether the full body is worth pulling. "
            "For metadata-only questions (author, referees, publication "
            "status) the preview is usually enough."
        ),
        parameters=parameters,
        call=_call,
    )


def _read_memo_tool(memo_root: str) -> Tool:
    def _call(doc_id: str) -> dict:
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
                retriever=None,
                collection=None,
                path_root=memo_root,
            )
            if got.get("text") or got.get("source") == "disk":
                return got
            last = got
        return last or {
            "source": "none", "doc_id": doc_id, "text": "",
            "error": "doc not found on disk",
        }

    parameters = {
        "type": "object",
        "properties": {
            "doc_id": {"type": "string",
                       "description": "e.g. 'BAM-00713' — filename suffix tolerated"},
        },
        "required": ["doc_id"],
    }
    return Tool(
        name="read_memo",
        description=(
            "Read the FULL markdown body from disk (~50 KB). Expensive — "
            "only call after `preview_memo` shows the doc really is "
            "relevant AND the answer needs body text (methodology, "
            "systematic uncertainties, table numbers). For author / "
            "referee / status questions, the preview is enough."
        ),
        parameters=parameters,
        call=_call,
    )


def _search_index_tool(memo_index_path: str) -> Tool:
    """grep INDEX.md by keyword — the cheapest first hop for discovery.

    INDEX.md is a single markdown table (ID | Title | Abstract lede | Bytes)
    with one row per doc. Regex is case-insensitive. Returns parsed rows so
    the agent doesn't have to re-parse `|`-separated markdown.
    """
    import re as _re
    _ROW_RE = _re.compile(r"^\|\s*(BAM-\d+)\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|\s*([\d,]+)\s*\|\s*$")

    def _call(pattern: str, max_results: int = 50) -> dict:
        idx = Path(memo_index_path)
        if not idx.is_file():
            return {"error": f"INDEX.md not found at {memo_index_path}", "hits": []}
        try:
            regex = _re.compile(pattern, _re.IGNORECASE)
        except _re.error as e:
            return {"error": f"invalid regex: {e}", "hits": []}
        hits: list[dict] = []
        try:
            text = idx.read_text(encoding="utf-8", errors="replace")
        except OSError as e:
            return {"error": str(e), "hits": []}
        for line in text.splitlines():
            if not line.startswith("| BAM-"):
                continue
            if not regex.search(line):
                continue
            m = _ROW_RE.match(line)
            if not m:
                continue
            hits.append({
                "doc_id": m.group(1),
                "title": m.group(2),
                "abstract_lede": m.group(3),
                "bytes": m.group(4),
            })
            if len(hits) >= int(max_results):
                break
        return {"pattern": pattern, "n_results": len(hits), "hits": hits}

    parameters = {
        "type": "object",
        "properties": {
            "pattern": {
                "type": "string",
                "description": (
                    "Case-insensitive Python regex (searched against each "
                    "INDEX.md row). Prefer specific tokens over broad OR "
                    "chains — 'Zc\\\\(3900\\\\)' is better than "
                    "'BESIII|Zc|charmonium'."
                ),
            },
            "max_results": {"type": "integer", "default": 50, "minimum": 1, "maximum": 200},
        },
        "required": ["pattern"],
    }
    return Tool(
        name="search_index",
        description=(
            "Search INDEX.md (the corpus TOC — one row per doc with id, "
            "title, abstract lede, bytes) by regex. **First hop for any "
            "'which papers about X' / 'list docs on Y' question.** Returns "
            "structured rows; no body-file reads."
        ),
        parameters=parameters,
        call=_call,
    )


def _read_summary_tool(memo_summaries_root: str) -> Tool:
    """Read summaries/<doc_id>.md — ~1 KB abstract + heading map per doc."""
    def _call(doc_id: str) -> dict:
        root = Path(memo_summaries_root)
        if not root.is_dir():
            return {"error": f"summaries root not found at {memo_summaries_root}"}
        # Accept `BAM-00097`, `BAM-00097.md`, or a bare id.
        base = doc_id[:-3] if doc_id.endswith(".md") else doc_id
        candidate = root / f"{base}.md"
        if not candidate.is_file():
            return {"error": "summary not found", "doc_id": base}
        try:
            text = candidate.read_text(encoding="utf-8", errors="replace")
        except OSError as e:
            return {"error": str(e), "doc_id": base}
        return {
            "doc_id": base,
            "path": str(candidate.relative_to(root)),
            "text": text[:8000],
            "bytes": len(text),
        }

    parameters = {
        "type": "object",
        "properties": {
            "doc_id": {
                "type": "string",
                "description": "e.g. 'BAM-00097'",
            },
        },
        "required": ["doc_id"],
    }
    return Tool(
        name="read_summary",
        description=(
            "Read summaries/<doc_id>.md (~1 KB abstract + section heading "
            "map). Perfect for 'what is this paper about' or previewing "
            "candidates found via `search_index`. Cheaper than "
            "`preview_memo` and much cheaper than `read_memo`."
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
                    "`argv[0]` in {grep,find,wc,head,tail,cat,ls,sort,"
                    "uniq,cut,tr}. No pipes / redirects / `-exec`."
                ),
            },
            "timeout_s": {
                "type": "integer", "default": timeout,
                "minimum": 1, "maximum": 60,
            },
        },
        "required": ["argv"],
    }
    return Tool(
        name="run_cmd",
        description=(
            "Read-only shell over the memo root. Use for counts, name globs, "
            "quick peeks not covered by the specialised tools."
        ),
        parameters=parameters,
        call=_call,
    )


def build_memo_search_tools(
    memo_root: str,
    *,
    memo_index_path: Optional[str] = None,
    memo_summaries_root: Optional[str] = None,
) -> list[Tool]:
    """The seven-tool bundle the `memo_search` subagent is allowed to call.

    `memo_index_path` and `memo_summaries_root` fall back to canonical
    locations relative to `memo_root` (parent dir / INDEX.md and
    parent dir / summaries) when not given. Missing files are handled at
    call time by returning an error dict rather than crashing at build.
    """
    parent = Path(memo_root).parent
    idx_path = memo_index_path or str(parent / "INDEX.md")
    summ_root = memo_summaries_root or str(parent / "summaries")
    return [
        _search_index_tool(idx_path),
        _read_summary_tool(summ_root),
        _grep_tool(memo_root),
        _list_memos_tool(memo_root),
        _preview_memo_tool(memo_root),
        _read_memo_tool(memo_root),
        _run_cmd_tool(memo_root),
    ]


MEMO_SEARCH_DESCRIPTION = (
    "Search the on-disk memo markdown corpus. Use for full-body content, "
    "exact tokens / LaTeX / numbers inside memos, file counts, or when you "
    "have a `doc_id` and need the whole document. NOT for author / status / "
    "journal metadata (use `milvus_search`) and NOT for HyperNews history "
    "(use `hypernews_db` when it lands)."
)


def build_memo_search_subagent(
    *,
    llm: Any,
    memo_root: str,
    max_steps: int = 6,
    memo_index_path: Optional[str] = None,
    memo_summaries_root: Optional[str] = None,
    rewriter: Any = None,
    memory: Any = None,
) -> Tool:
    """Factory: return the memo-search subagent as one Tool the parent calls.

    `rewriter`, when set, fans the parent's question out into N English
    candidates at the subagent entry — each candidate runs its own inner
    loop concurrently. See `make_subagent_tool` for the full contract.
    """
    return make_subagent_tool(
        name="memo_search",
        description=MEMO_SEARCH_DESCRIPTION,
        llm=llm,
        tools=build_memo_search_tools(
            memo_root,
            memo_index_path=memo_index_path,
            memo_summaries_root=memo_summaries_root,
        ),
        system_prompt=MEMO_SEARCH_SYSTEM_PROMPT,
        max_steps=max_steps,
        rewriter=rewriter,
        memory=memory,
    )


__all__ = [
    "build_memo_search_subagent",
    "build_memo_search_tools",
    "MEMO_SEARCH_SYSTEM_PROMPT",
    "MEMO_SEARCH_DESCRIPTION",
]
