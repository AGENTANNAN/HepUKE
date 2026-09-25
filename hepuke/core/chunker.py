"""Generic Markdown chunker + secondary section-aware splitter.

Ported from rag-besiii `rag/memo_chunker.py`, kept **domain-agnostic**:
  * `chunk_pdf_markdown()` — split by markdown H1-H3, sub-split large sections,
    merge tiny tails. Verbatim from rag-besiii.
  * `clean_content_for_store()` / `clean_content_for_embed()` — strip HTML/MD
    tables, images, bare URLs, display/inline LaTeX. Verbatim.
  * `make_embed_text(chunk, prefix)` — build the text fed to BGE-M3. The BESIII
    hard-coded `[BESIII Memo | {sec}]` prefix is now a parameter.
  * `SecondaryChunker` — section-aware size-bounded splitter with sidecar
    asset extraction ([TABLE:n] / [FIG:n]). Verbatim except the per-section
    `_SECONDARY_MAX_CHARS` map is a class attr that callers can override.
  * `SectionParser` Protocol + `PaperSectionChunker` — pluggable section
    identification. `BES3MemoParser` is **not** ported (it's BESIII-specific);
    keep it in rag-besiii or subclass here.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterable, List, Protocol, Tuple

# ---------------------------------------------------------------------------
# Primary chunker: split by heading + paragraph merge
# ---------------------------------------------------------------------------

MAX_CHUNK_CHARS = 1500
MIN_CHUNK_CHARS = 80


def chunk_pdf_markdown(markdown: str) -> List[str]:
    """Split raw markdown into vector-ready chunks.

    Strategy: split by H1-H3 headings; for oversize sections, sub-split by
    paragraphs; then merge tail chunks below `MIN_CHUNK_CHARS` back in.
    """
    if not markdown or not markdown.strip():
        return []

    sections = re.split(r"(?=^#{1,3} )", markdown, flags=re.MULTILINE)
    sections = [s.strip() for s in sections if s.strip()]

    chunks: list[str] = []
    for section in sections:
        if len(section) <= MAX_CHUNK_CHARS:
            chunks.append(section)
            continue
        paragraphs = re.split(r"\n{2,}", section)
        current = ""
        for para in paragraphs:
            para = para.strip()
            if not para:
                continue
            if len(current) + len(para) + 2 <= MAX_CHUNK_CHARS:
                current = (current + "\n\n" + para).strip()
            else:
                if current:
                    chunks.append(current)
                if len(para) > MAX_CHUNK_CHARS:
                    for i in range(0, len(para), MAX_CHUNK_CHARS):
                        chunks.append(para[i:i + MAX_CHUNK_CHARS])
                else:
                    current = para
        if current:
            chunks.append(current)

    merged: list[str] = []
    for chunk in chunks:
        if len(chunk) < MIN_CHUNK_CHARS and merged:
            merged[-1] = merged[-1] + "\n\n" + chunk
        elif len(chunk) >= MIN_CHUNK_CHARS:
            merged.append(chunk)
    return merged


# ---------------------------------------------------------------------------
# Content cleaners
# ---------------------------------------------------------------------------

_LATEX_INLINE = re.compile(r"\$[^$\n]+\$")
_LATEX_DISPLAY = re.compile(r"\$\$.*?\$\$", re.DOTALL)
_LEADING_NUM = re.compile(r"^\d{1,4}\s+(?=[A-Z])")
_MD_TABLE_SEP = re.compile(r"^\|[-|: ]+\|$")
_MD_TABLE_ROW = re.compile(r"^\|.*\|$")
_HTML_TABLE = re.compile(r"<table[\s>].*?</table>", re.IGNORECASE | re.DOTALL)
_HTML_TAG = re.compile(r"<[a-zA-Z/][^>]*>")
_IMAGE_LINK = re.compile(r"!\[.*?\]\(.*?\)", re.DOTALL)
_BARE_URL = re.compile(r"https?://\S+")


def clean_content_for_store(content: str) -> str:
    """Strip tables/images/bare URLs — keep human-readable prose."""
    content = _LATEX_DISPLAY.sub("", content)
    content = _HTML_TABLE.sub("", content)
    content = _IMAGE_LINK.sub("", content)

    lines: list[str] = []
    for line in content.split("\n"):
        stripped = line.strip()
        if _MD_TABLE_SEP.match(stripped) or _MD_TABLE_ROW.match(stripped):
            continue
        if _HTML_TAG.search(stripped):
            text_only = _HTML_TAG.sub(" ", stripped).strip()
            if len(text_only) < 10:
                continue
            line = text_only
            stripped = line.strip()
        line = _LEADING_NUM.sub("", line)
        stripped = line.strip()
        if _BARE_URL.fullmatch(stripped):
            continue
        if len(stripped) < 3:
            continue
        lines.append(line)
    return "\n".join(lines).strip()


def clean_content_for_embed(content: str) -> str:
    """Extra pass: drop inline `$...$` too. Only used to build embed text."""
    cleaned = clean_content_for_store(content)
    cleaned = _LATEX_INLINE.sub(" ", cleaned)
    return cleaned


def make_embed_text(chunk: dict, prefix: str = "") -> str:
    """Build the text fed to BGE-M3.

    prefix — free-form semantic tag, e.g. `"[BESIII Memo | {section_type}]"`
             or `"[ICLR paper | {section_type}]"`. `{section_type}`,
             `{section_title}`, `{title}`, `{doc_id}` are format-substituted
             from the chunk dict.
    """
    if prefix:
        try:
            prefix_str = prefix.format(**chunk)
        except KeyError:
            prefix_str = prefix
    else:
        prefix_str = ""
    sec_title = chunk.get("section_title", "")
    content = chunk.get("content", "")
    body = clean_content_for_embed(content)
    header = f"{prefix_str} {sec_title}".strip()
    return f"{header}\n{body}".strip() if header else body


# ---------------------------------------------------------------------------
# SecondaryChunker: size-bounded resplit + [TABLE:n]/[FIG:n] extraction
# ---------------------------------------------------------------------------

_SECONDARY_MAX_CHARS: dict[str, int] = {
    "systematic_uncertainty": 600,
    "event_selection": 800,
    "background": 800,
    "results": 800,
    "fit": 800,
    "detector": 900,
    "dataset": 900,
    "abstract": 1200,
    "introduction": 1000,
    "summary": 1000,
    "method": 900,
    "experiments": 900,
    "related_work": 900,
    "conclusion": 900,
}
_OVERFLOW_FACTOR = 1.15


class SecondaryChunker:
    """Section-aware size-bounded resplit + [TABLE:n]/[FIG:n] extraction.

    Instantiate with a custom `max_chars` dict to override per-section budgets
    (defaults are ICLR-paper + BESIII-memo friendly).
    """

    def __init__(self, max_chars: dict[str, int] | None = None) -> None:
        self.max_chars = {**_SECONDARY_MAX_CHARS, **(max_chars or {})}

    def run(self, v1_chunks: list[dict]) -> list[dict]:
        v2: list[dict] = []
        for chunk in v1_chunks:
            sec_type = chunk.get("section_type", "")
            raw = chunk.get("content", "")

            content, assets = self._extract_assets(raw)
            max_c = self.max_chars.get(sec_type)

            if max_c is None or len(content) <= max_c * _OVERFLOW_FACTOR:
                v2.append({**chunk, "content": content, "assets": assets})
                continue

            sub_texts = self._split(content, sec_type, max_c)
            intro = (
                self._extract_intro(content)
                if sec_type == "systematic_uncertainty"
                else ""
            )

            for sub_idx, sub_text in enumerate(sub_texts):
                if intro and sub_idx > 0 and not sub_text.startswith(intro):
                    sub_text = intro + "\n" + sub_text
                sub_assets = [
                    a for a in assets
                    if f"[TABLE:{a['id']}]" in sub_text
                    or f"[FIG:{a['id']}]" in sub_text
                ]
                v2.append({**chunk, "content": sub_text, "assets": sub_assets})
        return v2

    def _extract_assets(self, content: str) -> Tuple[str, list[dict]]:
        assets: list[dict] = []
        ctr = [0]

        def _repl_table(m):
            aid = ctr[0]; ctr[0] += 1
            assets.append({"id": aid, "type": "table", "raw": m.group(0)})
            return f"[TABLE:{aid}]"

        def _repl_fig(m):
            aid = ctr[0]; ctr[0] += 1
            url_m = re.search(r"\((https?://[^)]+)\)", m.group(0))
            assets.append({
                "id": aid, "type": "figure",
                "url": url_m.group(1) if url_m else "",
            })
            return f"[FIG:{aid}]"

        content = re.sub(
            r"<table[\s>].*?</table>", _repl_table,
            content, flags=re.IGNORECASE | re.DOTALL,
        )
        content = re.sub(
            r"!\[.*?\]\(https?://[^)]+\)", _repl_fig,
            content, flags=re.DOTALL,
        )
        return content, assets

    def _split(self, content: str, sec_type: str, max_c: int) -> list[str]:
        use_overlap = sec_type == "event_selection"
        para_only = sec_type in ("detector", "dataset")

        paras = [p.strip() for p in content.split("\n\n") if p.strip()]

        if all(len(p) <= max_c * _OVERFLOW_FACTOR for p in paras):
            return self._merge_paragraphs(paras, max_c, use_overlap)

        result: list[str] = []
        for para in paras:
            if len(para) <= max_c * _OVERFLOW_FACTOR:
                result.append(para)
            elif para_only:
                result.extend(self._hard_split(para, max_c))
            else:
                sents = self._split_sentences(para)
                result.extend(self._merge_sentences(sents, max_c, use_overlap))
        return self._remerge_tiny(result, max_c)

    def _split_sentences(self, text: str) -> list[str]:
        store: list[str] = []

        def _save(m):
            store.append(m.group(0))
            return f"\x00L{len(store) - 1}\x00"

        t = re.sub(r"\$\$.*?\$\$", _save, text, flags=re.DOTALL)
        t = re.sub(r"\$[^$\n]+\$", _save, t)
        parts = re.split(r"(?<=[a-zA-Z\]\)])\.\s+(?=[A-Z\[\(])", t)

        result: list[str] = []
        for part in parts:
            for i, orig in enumerate(store):
                part = part.replace(f"\x00L{i}\x00", orig)
            s = part.strip()
            if s:
                result.append(s)
        return result

    def _merge_paragraphs(self, paras: list[str], max_c: int, overlap: bool) -> list[str]:
        chunks: list[str] = []
        cur = ""
        for para in paras:
            if not cur:
                cur = para
            elif len(cur) + len(para) + 2 <= max_c:
                cur += "\n\n" + para
            else:
                chunks.append(cur)
                cur = (cur.split("\n\n")[-1] + "\n\n" + para) if overlap else para
        if cur:
            chunks.append(cur)
        return chunks

    def _merge_sentences(self, sents: list[str], max_c: int, overlap: bool) -> list[str]:
        chunks: list[str] = []
        cur: list[str] = []
        cur_len = 0
        for sent in sents:
            s_len = len(sent) + 2
            if cur_len + s_len <= max_c * _OVERFLOW_FACTOR or not cur:
                cur.append(sent)
                cur_len += s_len
            else:
                chunks.append(" ".join(cur))
                cur = ([cur[-1], sent] if overlap else [sent])
                cur_len = (len(cur[0]) + s_len + 2) if overlap else s_len
        if cur:
            chunks.append(" ".join(cur))
        return chunks

    def _hard_split(self, text: str, max_c: int) -> list[str]:
        parts: list[str] = []
        while len(text) > max_c * _OVERFLOW_FACTOR:
            cut = max_c
            while cut < len(text):
                if text[:cut].count("$") % 2 == 0:
                    break
                cut += 1
            parts.append(text[:cut])
            text = text[cut:].strip()
        if text:
            parts.append(text)
        return parts

    def _remerge_tiny(self, parts: list[str], max_c: int, min_c: int = 80) -> list[str]:
        result: list[str] = []
        for p in parts:
            if result and len(result[-1]) + len(p) + 2 <= max_c and len(p) < min_c:
                result[-1] += "\n\n" + p
            else:
                result.append(p)
        return result

    def _extract_intro(self, content: str) -> str:
        for line in content.split("\n"):
            s = line.strip()
            if s and not re.match(r"^[•\-\*]|\d+[\.\)]", s):
                return s[:300]
        return ""


# ---------------------------------------------------------------------------
# SectionParser Protocol + PaperSectionChunker (ICLR-flavored)
# ---------------------------------------------------------------------------

@dataclass
class Section:
    section_type: str
    section_number: str
    section_title: str
    body: str
    assets: list[dict] = field(default_factory=list)


class SectionParser(Protocol):
    def parse(self, markdown: str, doc_id: str, **base_meta: Any) -> list[dict]: ...


# ICLR / generic scientific paper section patterns.
# Order matters: earlier entries win when a heading matches multiple patterns.
_PAPER_SECTION_PATTERNS: list[tuple[str, str, str]] = [
    (r"\babstract\b", "abstract", "Abstract"),
    (r"\bintroduction\b", "introduction", "Introduction"),
    (r"\bbackground\b|\bpreliminar", "background", "Background"),
    (r"\brelated\s*work\b|\bprior\s*work\b|\bliterature\s*review\b",
     "related_work", "Related Work"),
    (r"\bmethod(s|ology)?\b|\bapproach\b|\bmodel\b|\barchitecture\b|\bframework\b",
     "method", "Method"),
    (r"\bexperiment(s|al)?\b|\bevaluation\b|\bbenchmarks?\b|\bsetup\b",
     "experiments", "Experiments"),
    (r"\bresults?\b|\bfindings?\b|\banalysis\b|\bablation\b",
     "results", "Results"),
    (r"\bdiscussion\b|\blimitations?\b|\bfuture\s*work\b",
     "discussion", "Discussion"),
    (r"\bconclusion\b|\bsummary\b",
     "conclusion", "Conclusion"),
    (r"\bappendi(x|ces)\b|\bsupplement", "appendix", "Appendix"),
    (r"\breferences?\b|\bbibliograph", "references", "References"),
    (r"\backnowledg(e|ement|ment)s?\b", "acknowledgements", "Acknowledgements"),
]


class PaperSectionChunker:
    """Split an ICLR-style paper markdown into section chunks.

    Heuristic: walk each `#{1,3}` heading, match against the section pattern
    list; the first heading that hits gives that block its `section_type`.
    Everything until the next matching heading becomes the body. Unmatched
    top-level headings default to `section_type='other'`. Unlike
    `BES3MemoParser` there is no title/author/TOC extraction — pass `title` /
    `authors` / `date` via `base_meta`.
    """

    _HEADING_RE = re.compile(r"^(#{1,3})\s+(.+)$", re.MULTILINE)
    _TRAILING_PAGE_NUM = re.compile(r"\s+\d+\s*$")

    def _classify(self, heading: str) -> tuple[str, str, str] | None:
        text = self._TRAILING_PAGE_NUM.sub("", heading).strip()
        num_match = re.match(r"^\s*(\d+(?:\.\d+)*)\s+(.*)$", text)
        sec_num = num_match.group(1) if num_match else ""
        tail = num_match.group(2) if num_match else text
        for pat, sec_type, canonical in _PAPER_SECTION_PATTERNS:
            if re.search(pat, tail, re.IGNORECASE):
                clean = tail.strip() or canonical
                return sec_type, sec_num, clean
        return None

    def parse(self, markdown: str, doc_id: str, **base_meta: Any) -> list[dict]:
        if not markdown.strip():
            return []
        matches = list(self._HEADING_RE.finditer(markdown))
        if not matches:
            return [{
                **base_meta, "doc_id": doc_id,
                "section_number": "", "section_title": "Body",
                "section_type": "other", "content": markdown.strip(),
            }]

        chunks: list[dict] = []
        for i, m in enumerate(matches):
            heading = m.group(2).strip()
            body_start = m.end()
            body_end = matches[i + 1].start() if i + 1 < len(matches) else len(markdown)
            body = markdown[body_start:body_end].strip()
            if not body:
                continue
            cls = self._classify(heading)
            if cls is None:
                sec_type, sec_num, sec_title = "other", "", heading
            else:
                sec_type, sec_num, sec_title = cls
            chunks.append({
                **base_meta,
                "doc_id": doc_id,
                "section_number": sec_num,
                "section_title": sec_title,
                "section_type": sec_type,
                "content": body,
            })
        return chunks
