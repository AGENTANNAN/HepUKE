"""Retrieval interface over the ``besiii_dsl_desc_v1`` collection.

``search_dsl(query, top_k=5)`` performs a hybrid search (dense + sparse) over
the DSL description index, then re-ranks with the configured reranker if
enabled. Each hit is resolved to its gold ``.rb`` path via the ``stem``
metadata written by ``build_index.py``.

Only import this module from generate_candidates+ scripts; it does NOT install any gating
side-effects into the main RAG loop.

Returned schema per hit:
    {"stem": str, "score": float, "description": str, "rb_path": str,
     "rb_text": str}
"""
from __future__ import annotations

import json
import sys
import threading
from functools import lru_cache
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hepuke import HepUKE, load_config  # noqa: E402


DEFAULT_COLLECTION = "besiii_dsl_desc_v1"

# Module-level lock protecting BOTH (a) the first-time HepUKE init inside
# _rag() and (b) every retriever.search() call. FlagEmbedding's BGEM3
# tokenizer.pad() is not thread-safe: two workers calling encode() at the
# same time race and one gets `'list' object has no attribute 'keys'`.
# Serialising search here matches the encoding-lock already present in
# hepuke.core.embedding.BgeM3LocalBackend and keeps concurrent generate_candidates workers
# safe without duplicating the 4 GB model per thread.
_search_lock = threading.Lock()


_rag_instance: dict[str, HepUKE] = {}
_init_lock = threading.Lock()


def _rag(collection: str = DEFAULT_COLLECTION) -> HepUKE:
    if collection in _rag_instance:
        return _rag_instance[collection]
    with _init_lock:
        if collection in _rag_instance:
            return _rag_instance[collection]
        cfg = load_config()
        cfg.milvus.default_collection = collection
        rag = HepUKE(cfg)
        rag.connect_collection(collection)
        _rag_instance[collection] = rag
        return rag


def _resolve(hit: dict[str, Any]) -> dict[str, Any] | None:
    try:
        meta = json.loads(hit.get("metadata_json") or "{}")
    except Exception:
        meta = {}
    stem = meta.get("stem") or hit.get("doc_id")
    rb_path = meta.get("rb_path")
    if not stem or not rb_path:
        return None
    rb_file = Path(rb_path)
    if not rb_file.exists():
        return None
    return {
        "stem": stem,
        "score": float(hit.get("score", 0.0) or 0.0),
        "description": hit.get("content") or "",
        "rb_path": rb_path,
        "rb_text": rb_file.read_text(encoding="utf-8"),
    }


def search_dsl(
    query: str,
    top_k: int = 5,
    collection: str = DEFAULT_COLLECTION,
    with_reranker: bool = True,
    exclude_stems: set[str] | None = None,
) -> list[dict[str, Any]]:
    """Return up to ``top_k`` DSL candidates for ``query``.

    ``with_reranker`` defaults True — must be passed to retriever.search()
    every call because HybridRetriever.search() itself defaults to
    ``with_rerank=False`` regardless of config.retrieval.with_reranker.

    ``exclude_stems`` — never return these arXiv stems (leakage guard). Because
    the KB currently contains every gold DSL (868 docs, including all 267
    dev+test targets), callers MUST pass at least ``{qid}`` when scoring a
    query whose own gold is the answer — otherwise the model will retrieve its
    own gold as top-k and the RAG comparison becomes vacuous. We over-fetch by
    ``len(exclude_stems)`` slots to keep ``top_k`` real hits after filtering.
    """
    excl = set(exclude_stems or ())
    rag = _rag(collection)
    retriever = rag.retriever
    fetch_k = top_k + len(excl) if excl else top_k
    with _search_lock:
        hits = retriever.search(query, collection,
                                top_k=fetch_k, with_rerank=with_reranker)
    resolved: list[dict[str, Any]] = []
    for h in hits:
        r = _resolve(h)
        if r is None:
            continue
        if r["stem"] in excl:
            continue
        resolved.append(r)
        if len(resolved) >= top_k:
            break
    return resolved


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("query")
    ap.add_argument("--k", type=int, default=5)
    args = ap.parse_args()
    for i, r in enumerate(search_dsl(args.query, top_k=args.k), 1):
        print(f"[{i}] stem={r['stem']} score={r['score']:.4f} rb={r['rb_path']}")
        print(f"    desc[:200]={r['description'][:200]!r}")
