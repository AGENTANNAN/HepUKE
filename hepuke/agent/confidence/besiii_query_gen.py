"""BESIII calibration-query generator (T7b, upstream of T7 collector).

Turns a directory of BESIII memo markdowns (or a Milvus title dump) into a
``queries.jsonl`` file whose rows feed :func:`collect_trajectories`. Three
query types are produced in the mix demanded by
``entropy/section_method.tex§4.3`` (also §5-hep):

* ``numeric``       — asks for a specific number / process from ONE memo
* ``methodological``— asks how a particular measurement was performed
* ``synthesis``     — asks to compare / summarise across ≥2 memos

The generator is deliberately template-based (no LLM cost) so it can be
regenerated cheaply whenever the memo corpus grows. Downstream, the
collector runs the real agent on each query — that's where the actual
retrieval trajectories and posteriors come from.
"""

from __future__ import annotations

import json
import random
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Literal, Sequence


QueryType = Literal["numeric", "methodological", "synthesis"]


# ---------------------------------------------------------------------------
# memo parsing
# ---------------------------------------------------------------------------

_DOC_ID_RE = re.compile(r"^(BAM-\d{5,6})", re.IGNORECASE)
_TITLE_RE = re.compile(r"^\s*#\s+(.+?)\s*$")


@dataclass
class MemoEntry:
    doc_id: str
    filename: str
    title: str
    keywords: list[str] = field(default_factory=list)


def _clean_title(raw: str) -> str:
    """Strip LaTeX math markers and collapse whitespace."""
    # LaTeX inline math $...$ → keep body without $s so the plain text
    # is still readable in a query.
    t = re.sub(r"\$([^$]+)\$", r"\1", raw)
    t = re.sub(r"[{}\\]", " ", t)  # \psi, ^ { \prime }, etc.
    t = re.sub(r"\s+", " ", t).strip()
    return t


def _extract_keywords(title: str, filename: str) -> list[str]:
    """Cheap keyword harvest: (a) capitalised particle names ``psi``,
    ``chi_c`` etc.; (b) hyphenated slugs in the filename."""
    kws: list[str] = []
    for m in re.finditer(r"[A-Za-z][a-z]{2,}", title):
        w = m.group(0)
        if len(w) >= 4 and w.lower() not in {"decay", "decays", "measurement", "measurements"}:
            kws.append(w)
    slug = Path(filename).stem
    for chunk in slug.split("_"):
        chunk = chunk.strip("-")
        if len(chunk) >= 4 and not chunk.startswith("BAM"):
            kws.append(chunk)
    # De-dup preserving order.
    seen: set[str] = set()
    out: list[str] = []
    for w in kws:
        wl = w.lower()
        if wl in seen:
            continue
        seen.add(wl)
        out.append(w)
    return out[:5]


def parse_memo_head(path: Path) -> MemoEntry | None:
    """Read the first ~40 lines of a memo and return (doc_id, title, kws).

    Returns None when the file has no BAM-* prefix — the generator is
    deliberately conservative and skips anything that doesn't fit the
    BESIII naming convention.
    """
    m = _DOC_ID_RE.match(path.name)
    if not m:
        return None
    doc_id = m.group(1).upper()
    title = ""
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f):
                if i > 40:
                    break
                tm = _TITLE_RE.match(line)
                if tm:
                    title = _clean_title(tm.group(1))
                    break
    except OSError:
        return None
    if not title:
        # Fallback: reconstruct a short title from the filename slug.
        slug = path.stem.split("_", 1)[-1] if "_" in path.stem else path.stem
        title = slug.replace("-", " ").replace("_", " ")
    return MemoEntry(
        doc_id=doc_id,
        filename=path.name,
        title=title,
        keywords=_extract_keywords(title, path.name),
    )


def scan_memo_root(root: str | Path, *, limit: int | None = None) -> list[MemoEntry]:
    """Walk ``root`` for ``BAM-*.md`` files and parse their headers."""
    p = Path(root)
    if not p.exists():
        return []
    entries: list[MemoEntry] = []
    for fp in sorted(p.rglob("BAM-*.md")):
        e = parse_memo_head(fp)
        if e is None:
            continue
        entries.append(e)
        if limit is not None and len(entries) >= limit:
            break
    return entries


# ---------------------------------------------------------------------------
# query templates
# ---------------------------------------------------------------------------

_NUMERIC_TEMPLATES = [
    "在 {doc_id} 中，{title} 的中心值是多少？请给出数值和单位。",
    "{doc_id}（{title}）报道的主要测量结果是什么？只需要给出数值和统计误差。",
    "查阅 {doc_id}，说明 {title} 的最终结果，并附引用编号。",
]

_METHOD_TEMPLATES = [
    "{doc_id} 中 {title} 采用的分析方法有哪些步骤？请概括流程。",
    "在 {doc_id} 里，{title} 的系统误差是如何评估的？",
    "{doc_id} 中 {title} 使用了怎样的数据样本与背景处理策略？",
]

_SYNTHESIS_TEMPLATES = [
    "对比 {doc_id_a} 与 {doc_id_b} 中关于 {topic} 的测量结果，指出主要差异。",
    "综合 {doc_id_a}、{doc_id_b} 两篇 memo，评述 BESIII 在 {topic} 方向的进展。",
    "结合 {doc_id_a} 和 {doc_id_b}，讨论 {topic} 的分析策略是否一致。",
]


def _pick_topic(a: MemoEntry, b: MemoEntry) -> str:
    for kw in a.keywords:
        if kw in b.keywords or kw.lower() in b.title.lower():
            return kw
    # Fallback to a.title first-word bigram.
    words = [w for w in a.title.split() if len(w) >= 3]
    return " ".join(words[:2]) if words else a.title


def render_numeric(entry: MemoEntry, rng: random.Random) -> str:
    tpl = rng.choice(_NUMERIC_TEMPLATES)
    return tpl.format(doc_id=entry.doc_id, title=entry.title)


def render_methodological(entry: MemoEntry, rng: random.Random) -> str:
    tpl = rng.choice(_METHOD_TEMPLATES)
    return tpl.format(doc_id=entry.doc_id, title=entry.title)


def render_synthesis(a: MemoEntry, b: MemoEntry, rng: random.Random) -> str:
    tpl = rng.choice(_SYNTHESIS_TEMPLATES)
    return tpl.format(doc_id_a=a.doc_id, doc_id_b=b.doc_id, topic=_pick_topic(a, b))


# ---------------------------------------------------------------------------
# top-level generator
# ---------------------------------------------------------------------------

@dataclass
class GenSpec:
    """How many of each type to produce.

    Defaults follow the paper's stated calibration mix (numeric-heavy,
    smaller share of synthesis because it costs more per trajectory).
    """
    n_numeric: int = 200
    n_methodological: int = 200
    n_synthesis: int = 100
    seed: int = 0


def generate_queries(
    entries: Sequence[MemoEntry],
    spec: GenSpec = GenSpec(),
) -> Iterator[dict]:
    """Yield ``{id, query, type, doc_ids}`` records.

    ``doc_ids`` lists the memo(s) the template was drawn from — useful as
    a *weak* relevance label for later evaluation (not used by T7 itself,
    which is self-supervised on the posterior stream).
    """
    if not entries:
        return

    rng = random.Random(spec.seed)

    def _emit(qid: str, text: str, qtype: QueryType, doc_ids: list[str]):
        yield {
            "id": qid,
            "query": text,
            "type": qtype,
            "doc_ids": doc_ids,
        }

    numeric = list(entries)
    rng.shuffle(numeric)
    for i in range(spec.n_numeric):
        e = numeric[i % len(numeric)]
        yield from _emit(
            qid=f"num-{i:04d}-{e.doc_id}",
            text=render_numeric(e, rng),
            qtype="numeric",
            doc_ids=[e.doc_id],
        )

    method = list(entries)
    rng.shuffle(method)
    for i in range(spec.n_methodological):
        e = method[i % len(method)]
        yield from _emit(
            qid=f"met-{i:04d}-{e.doc_id}",
            text=render_methodological(e, rng),
            qtype="methodological",
            doc_ids=[e.doc_id],
        )

    if len(entries) >= 2:
        pool = list(entries)
        for i in range(spec.n_synthesis):
            a = rng.choice(pool)
            b = rng.choice(pool)
            attempts = 0
            while b.doc_id == a.doc_id and attempts < 8:
                b = rng.choice(pool)
                attempts += 1
            if b.doc_id == a.doc_id:
                continue
            yield from _emit(
                qid=f"syn-{i:04d}-{a.doc_id}-{b.doc_id}",
                text=render_synthesis(a, b, rng),
                qtype="synthesis",
                doc_ids=[a.doc_id, b.doc_id],
            )


def dump_queries(rows: Iterable[dict], out_path: str | Path) -> int:
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with p.open("w", encoding="utf-8") as f:
        for r in rows:
            f.write(json.dumps(r, ensure_ascii=False))
            f.write("\n")
            n += 1
    return n


__all__ = [
    "MemoEntry",
    "GenSpec",
    "parse_memo_head",
    "scan_memo_root",
    "render_numeric",
    "render_methodological",
    "render_synthesis",
    "generate_queries",
    "dump_queries",
]
