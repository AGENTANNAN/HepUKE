"""HepUKE — facade over HybridRetriever + MilvusStore + LLM.

Rewritten in P1 as a thin wrapper on top of `HybridRetriever` (BGE-M3 hybrid
+ RRF + small-to-big). llama-index is gone; the public method surface is
preserved so FastAPI / smoke tests / quickstart keep working:

  * `get_collections()` / `create_collection()` / `drop_collection()`
  * `connect_collection(name)`
  * `insert(document, split_chunk)` — chunks + embeds + inserts
  * `retrieve(content, ...)` — returns list of hit dicts (native shape)
  * `query(content, ...)` — retrieve → LLM answer
  * `get_doc(key, value)` — flat scalar lookup
  * `delete(doc_id)` — best-effort delete by `doc_id`
  * `persist()` — no-op (Milvus is durable); kept for API compat
"""
from __future__ import annotations

import json
import logging
from hashlib import sha256
from typing import Any, List, Tuple

from hepuke.config import Config
from hepuke.core.chunker import (
    MAX_CHUNK_CHARS,
    chunk_pdf_markdown,
    clean_content_for_store,
)
from hepuke.core.embedding import build_embedding
from hepuke.core.hybrid_retriever import DEFAULT_OUTPUT_FIELDS, HybridRetriever
from hepuke.core.llm import build_llm, LLMClient
from hepuke.core.query_preprocessor import preprocess
from hepuke.core.query_rewriter import LLMQueryRewriter
from hepuke.core.reranker import build_reranker
from hepuke.core.schema_profile import SchemaProfile, detect_schema_profile
from hepuke.core.vector_store import MilvusStore
from hepuke.document import Document
from hepuke.observability.trace import configure_trace

logger = logging.getLogger(__name__)


class HepUKE:
    def __init__(self, config: Config) -> None:
        self.config = config
        configure_trace(
            enabled=config.trace.enabled,
            file=config.trace.file,
            include_context=config.trace.include_context,
            backend=config.trace.backend,
            db_path=config.trace.db_path,
            console=config.trace.console,
        )
        self.embedder = build_embedding(config.embedding)
        self.reranker = build_reranker(config.reranker)
        self.llm = build_llm(config.llm)
        self.store = MilvusStore(
            uri=config.milvus.uri,
            token=config.milvus.token,
            database=config.milvus.database,
            dim=config.embedding.dim,
            text_search_mode=config.retrieval.text_search_mode,
        )
        self.retriever = HybridRetriever(
            store=self.store,
            embedder=self.embedder,
            reranker=self.reranker,
            preprocess=preprocess if config.retrieval.with_query_preprocess else None,
            rrf_k=config.retrieval.rrf_k,
        )
        # LLM query rewriter (multi-candidate, RRF-fused). When enabled,
        # subagents pass `query_pairs=` to retriever.search so the corpus
        # gets searched under multiple English candidates for each
        # question — this is what closes the CN-question / EN-corpus gap
        # and the entity-vs-term-family recall gap (Zc(3900) ↔
        # charmonium-like). Disabled by default; users opt in via
        # config.retrieval.query_rewrite_llm_enabled. When disabled the
        # attr stays None so downstream `if self.rewriter is not None`
        # gates work uniformly.
        self.rewriter: LLMQueryRewriter | None = None
        if getattr(config.retrieval, "query_rewrite_llm_enabled", False):
            try:
                # Dedicated LLM client for the rewriter — usually a fast/
                # cheap model (glm-5.2, deepseek-v4-flash) rather than the
                # heavy agent LLM. Falls back to self.llm when model is
                # empty to preserve the old single-LLM contract.
                rewriter_model = str(
                    getattr(config.retrieval, "query_rewrite_llm_model", "") or ""
                ).strip()
                if rewriter_model:
                    rewriter_llm = LLMClient(
                        base_url=config.llm.base_url,
                        api_key=config.llm.api_key,
                        model=rewriter_model,
                        temperature=0.0,
                        max_tokens=1024,
                    )
                else:
                    rewriter_llm = self.llm
                self.rewriter = LLMQueryRewriter(
                    llm=rewriter_llm,
                    cache_size=int(
                        getattr(config.retrieval, "query_rewrite_cache_size", 512),
                    ),
                    max_candidates=int(
                        getattr(config.retrieval, "query_rewrite_n_candidates", 4),
                    ),
                    keep_original=bool(
                        getattr(config.retrieval, "query_rewrite_keep_original", True),
                    ),
                    prompt_style=str(
                        getattr(config.retrieval, "query_rewrite_prompt_style", "physics"),
                    ),
                )
            except Exception as e:  # pragma: no cover
                logger.warning(
                    "LLMQueryRewriter init failed (%s); falling back to "
                    "single-query retrieval",
                    e,
                )
                self.rewriter = None
        self.collection_name: str | None = None
        self.schema_profile: SchemaProfile | None = None
        # Preload the configured default collection so the first request
        # doesn't pay the connect + Milvus `load_collection` cost (100ms-2s).
        # Errors are logged but don't crash startup — a bad config or a
        # missing collection should still let the health endpoint answer.
        default_coll = config.milvus.default_collection
        if default_coll:
            try:
                self.connect_collection(default_coll)
            except Exception as e:  # pragma: no cover — depends on live Milvus
                logger.warning(
                    "preload of default collection %r failed: %s (first "
                    "request will retry)", default_coll, e,
                )

    # ---- collections -------------------------------------------------------
    def get_collections(self) -> List[str]:
        return self.store.list_collections()

    def create_collection(self, name: str) -> List[str]:
        self.store.create_hybrid_collection(name)
        return self.get_collections()

    def drop_collection(self, name: str) -> List[str]:
        self.store.drop_collection(name)
        if self.collection_name == name:
            self.collection_name = None
        return self.get_collections()

    def connect_collection(self, name: str) -> bool:
        if not self.store.has_collection(name):
            self.store.create_hybrid_collection(name)
        self.store.load_collection(name)
        self.collection_name = name
        self._apply_schema_profile(name)
        logger.info(
            "connected collection %s (doc_id_field=%s, content_field=%s, hybrid=%s, chunking=%s)",
            name,
            self.schema_profile.doc_id_field,
            self.schema_profile.primary_content_field,
            self.schema_profile.sparse_field is not None,
            self.schema_profile.has_chunking,
        )
        return True

    def _apply_schema_profile(self, name: str) -> None:
        """Introspect the connected collection and reconfigure `self.retriever`
        so it uses the actual field names (`bam_id` vs `doc_id`, `text` vs
        `content`, dense-only vs hybrid) and the actual dense metric type.
        Called once per `connect_collection()`."""
        profile = detect_schema_profile(self.store, name)
        self.schema_profile = profile
        self.retriever.doc_id_field = profile.doc_id_field
        self.retriever.content_field = profile.primary_content_field
        if profile.dense_field:
            self.retriever.dense_field = profile.dense_field
        self.retriever.sparse_field = profile.sparse_field  # may be None
        self.retriever.dense_metric = profile.dense_metric
        # Always use the schema-derived output_fields — the P0 default list
        # references columns (`doc_id`, `metadata_json`, ...) that don't exist
        # on other layouts, and passing a non-existent output_field to
        # hybrid_search is a hard error in Milvus.
        self.retriever.default_output_fields = profile.output_fields

    def _require_collection(self) -> None:
        if self.collection_name is None:
            raise RuntimeError("no collection connected; call connect_collection() first")

    # ---- documents ---------------------------------------------------------
    def insert(
        self,
        document: Document,
        split_chunk: bool = False,
    ) -> Tuple[str, str | None]:
        self._require_collection()
        text = document.text or ""
        if not text.strip():
            return "empty document", None
        if split_chunk and len(text) > MAX_CHUNK_CHARS * 20:
            return f"chunk_size {len(text)} > {MAX_CHUNK_CHARS*20}", None

        doc_id = document.doc_id or _hash_doc_id(text)
        doc_hash = sha256(text.encode("utf-8", "surrogatepass")).hexdigest()
        meta = dict(document.metadata or {})
        meta.setdefault("hash", doc_hash)

        # If a doc with this hash already exists, skip.
        existing = self.get_doc("metadata_hash", doc_hash)
        if existing:
            return "doc with same hash already exists", None

        chunk_texts: list[str] = (
            chunk_pdf_markdown(text) if split_chunk else [text]
        )
        rows: list[dict] = []
        for i, ct in enumerate(chunk_texts):
            rows.append({
                "doc_id": doc_id,
                "chunk_type": "child",
                "parent_id": f"{doc_id}::doc::0",
                "section_type": meta.get("section_type", "body"),
                "section_number": str(i),
                "section_title": meta.get("section_title", ""),
                "title": meta.get("title", ""),
                "authors": meta.get("authors", ""),
                "date": meta.get("date", ""),
                "source_url": meta.get("source_url", ""),
                "content": clean_content_for_store(ct) or ct,
                "assets": "",
                "metadata_json": json.dumps(meta, ensure_ascii=False),
            })
        # Add one parent chunk with the full text (so small-to-big and
        # `get_full_doc` behave predictably).
        rows.append({
            "doc_id": doc_id,
            "chunk_type": "parent",
            "parent_id": f"{doc_id}::doc::0",
            "section_type": meta.get("section_type", "body"),
            "section_number": "0",
            "section_title": meta.get("section_title", ""),
            "title": meta.get("title", ""),
            "authors": meta.get("authors", ""),
            "date": meta.get("date", ""),
            "source_url": meta.get("source_url", ""),
            "content": clean_content_for_store(text) or text,
            "assets": "",
            "metadata_json": json.dumps(meta, ensure_ascii=False),
        })
        n = self.retriever.upsert_chunks(self.collection_name, rows)
        return f"inserted {n} rows", doc_id

    def delete(self, doc_id: str) -> bool | str:
        self._require_collection()
        doc_id_field = (
            self.schema_profile.doc_id_field if self.schema_profile else "doc_id"
        )
        try:
            self.store.load_collection(self.collection_name)
            self.store.client.delete(
                collection_name=self.collection_name,
                filter=f'{doc_id_field} == "{doc_id}"',
            )
            self.store.client.flush(self.collection_name)
            return True
        except Exception as exc:  # pragma: no cover
            logger.warning("delete(%s) failed: %s", doc_id, exc)
            return f"delete failed: {exc}"

    def update(self, document: Document) -> bool:
        """Delete-then-insert by `doc_id`."""
        if document.doc_id:
            self.delete(document.doc_id)
        _, _ = self.insert(document)
        return True

    # ---- retrieval / query -------------------------------------------------
    def retrieve(
        self,
        content: str,
        *,
        similarity_top_k: int | None = None,
        filters: List[dict] | None = None,
        condition: str = "and",
        with_reranker: bool | None = None,
    ) -> List[dict]:
        self._require_collection()
        top_k = similarity_top_k or self.config.retrieval.similarity_top_k
        extra_filter = _filters_to_expr(filters or [], condition)
        use_rerank = (
            self.config.retrieval.with_reranker if with_reranker is None else with_reranker
        )
        hits = self.retriever.search(
            question=content,
            collection=self.collection_name,
            top_k=top_k,
            extra_filter=extra_filter,
            with_rerank=use_rerank,
            with_small_to_big=self.config.retrieval.with_small_to_big,
            per_doc_dedup=self.config.retrieval.per_doc_dedup,
        )
        return hits

    def query(self, content: str, **retrieve_kwargs: Any) -> dict:
        hits = self.retrieve(content, **retrieve_kwargs)
        content_field = (
            self.schema_profile.primary_content_field
            if self.schema_profile
            else "content"
        )
        doc_id_field = (
            self.schema_profile.doc_id_field if self.schema_profile else "doc_id"
        )
        contexts = [h.get(content_field, "") for h in hits]
        answer = self.llm.chat(content, contexts=contexts)
        return {
            "answer": answer,
            "sources": [
                {
                    "score": h.get("score"),
                    "text": h.get(content_field, ""),
                    "metadata": {
                        k: h.get(k) for k in
                        (doc_id_field, "section_type", "section_title", "title")
                        if h.get(k) is not None
                    },
                }
                for h in hits
            ],
        }

    # ---- raw doc lookup ----------------------------------------------------
    def get_doc(self, key: str, value: str | dict) -> List[dict]:
        self._require_collection()
        collection = self.collection_name
        try:
            if key == "metadata":
                if not isinstance(value, dict):
                    return []
                pieces = [f'{k} == "{v}"' for k, v in value.items()]
                expr = " and ".join(pieces)
            elif key == "metadata_hash":
                # metadata_json contains `"hash": "<sha>"` — do a substring match.
                expr = f'metadata_json like "%\\\"hash\\\": \\\"{value}\\\"%"'
            elif key == "search":
                if not isinstance(value, str):
                    return []
                expr = self.store.build_text_expr(value)
            else:
                expr = f'{key} == "{value}"'
            return self.store.query_by_expr(collection, expr, limit=1000)
        except Exception as exc:
            logger.warning("get_doc(%s=%r) failed: %s", key, value, exc)
            return []

    def persist(self) -> None:
        """No-op: Milvus persists on flush(). Kept for API compatibility."""
        if self.collection_name is None:
            return
        try:
            self.store.client.flush(self.collection_name)
        except Exception:
            pass


# Backwards-compat alias.
MyIndex = HepUKE


def _hash_doc_id(text: str) -> str:
    return "doc_" + sha256(text.encode("utf-8", "surrogatepass")).hexdigest()[:16]


def _filters_to_expr(filters: list[dict], condition: str) -> str:
    """Translate the old llama-index filter list into a Milvus expr."""
    if not filters:
        return ""
    pieces: list[str] = []
    for f in filters:
        key = f.get("key")
        value = f.get("value")
        op = f.get("operator", "==")
        if key is None:
            continue
        if isinstance(value, str):
            pieces.append(f'{key} {op} "{value}"')
        else:
            pieces.append(f"{key} {op} {value}")
    joiner = " and " if condition.lower() == "and" else " or "
    return joiner.join(pieces)
