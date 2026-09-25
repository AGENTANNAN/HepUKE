"""Hybrid retriever: BGE-M3 dense+sparse + RRF + rerank + small-to-big.

Uses `pymilvus.MilvusClient` throughout (no ORM `Collection` calls). The
schema-aware knobs (`doc_id_field` / `content_field` / `dense_field` /
`sparse_field` / `dense_metric` / `default_output_fields`) are set by
`HepUKE.connect_collection()` via `SchemaProfile`, so this class stays
schema-agnostic and works against heterogeneous collections (P0 hybrid
chunking, dense-only `test`, sparse-less variants, `bam_id`/`paper_id`/`id`
primary keys, ...).

Depends on:
  * `hepuke.core.embedding.EmbeddingBackend`
  * `hepuke.core.vector_store.MilvusStore`  (MilvusClient wrapper)
  * `hepuke.core.reranker.Reranker`         (optional)
  * `hepuke.observability.trace.emit`
"""
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Iterable, List, Literal, Optional, Sequence

from hepuke.observability.trace import emit


DEFAULT_OUTPUT_FIELDS: list[str] = [
    "pk", "doc_id", "chunk_type", "parent_id",
    "section_type", "section_number", "section_title",
    "title", "authors", "date", "source_url",
    "content", "assets", "metadata_json",
]


def _rrf_merge(hits_lists: Sequence[Sequence[dict]], *, k: int = 60) -> list[dict]:
    """Reciprocal-rank-fusion merge of several ranked lists into one.

    Score = sum over lists of 1/(k + rank). Same primary key across lists
    accumulates; unique keys just get their single-list contribution. The
    fused list is returned sorted by fused score descending, with each
    hit carrying `score` overwritten to the fused value (so downstream
    dedup/rerank sees a consistent scoring).

    Dedup key uses `pk` when present, else the compound
    `(doc_id_field-hit-uses, parent_id, content[:64])` fallback — the
    caller-agnostic bit is that anything without a stable id degrades
    gracefully to per-record.
    """
    scores: dict[Any, float] = {}
    rep: dict[Any, dict] = {}
    for hits in hits_lists:
        for rank, h in enumerate(hits):
            key = h.get("pk")
            if key is None:
                key = (
                    h.get("doc_id"), h.get("parent_id"),
                    (h.get("content") or "")[:64],
                )
            contribution = 1.0 / (k + rank + 1)
            scores[key] = scores.get(key, 0.0) + contribution
            # Keep the first-seen row; RRF cares about ranks not fields.
            rep.setdefault(key, h)
    fused: list[dict] = []
    for key, s in sorted(scores.items(), key=lambda kv: kv[1], reverse=True):
        row = dict(rep[key])
        row["score"] = s
        fused.append(row)
    return fused


class HybridRetriever:
    """BGE-M3 dense+sparse hybrid retriever with RRF fusion.

    All Milvus calls go through `store.client` (a `MilvusClient`). Schema
    knobs are overwritten per collection by `HepUKE._apply_schema_profile()`.
    """

    def __init__(
        self,
        *,
        store,
        embedder,
        reranker=None,
        preprocess: Optional[Callable[[str], tuple[str, str]]] = None,
        doc_id_field: str = "doc_id",
        content_field: str = "content",
        dense_field: str = "dense_vector",
        sparse_field: Optional[str] = "sparse_vector",
        dense_metric: str = "COSINE",
        rrf_k: int = 60,
        multi_query_max_parallel: int = 8,
    ) -> None:
        self.store = store
        self.embedder = embedder
        self.reranker = reranker
        self._preprocess = preprocess
        self.doc_id_field = doc_id_field
        self.content_field = content_field
        self.dense_field = dense_field
        self.sparse_field = sparse_field
        self.dense_metric = dense_metric
        self.rrf_k = rrf_k
        # Wall-clock cap on concurrent Milvus hits when `search(query_pairs=..)`
        # is called with multiple candidates. Milvus tolerates well over 8
        # concurrent client requests; this cap protects the local thread
        # pool, not the server.
        self.multi_query_max_parallel = max(1, int(multi_query_max_parallel))
        # Per-collection output_fields override. When None, callers/`search()`
        # fall back to `DEFAULT_OUTPUT_FIELDS` (the P0 hybrid layout).
        self.default_output_fields: Optional[List[str]] = None

    # ---------------- helpers ----------------------------------------------

    def _schema_field_names(self, collection: str) -> set[str]:
        desc = self.store.describe_collection(collection)
        return {f["name"] for f in desc.get("fields", [])}

    @staticmethod
    def _row_from_hit(hit: dict, collection: str) -> dict:
        """Normalize a MilvusClient hit into `{score, collection, <fields>...}`.

        Client-mode hits look like `{'id': ..., 'distance': 0.42,
        'entity': {'field': value, ...}}`. Older code assumed `.score` and
        `.fields`; we adapt here.
        """
        entity = hit.get("entity") or {}
        record: dict = {
            "score": float(hit.get("distance", 0.0)),
            "collection": collection,
        }
        record.update(entity)
        return record

    # ---------------- core hybrid search ----------------------------------

    def _hybrid_search(
        self,
        collection_name: str,
        top_k: int,
        output_fields: List[str],
        dense_text: str,
        sparse_text: str,
        filter_expr: str = "",
        partition_names: Optional[List[str]] = None,
        channel: Literal["dense", "sparse", "hybrid"] = "hybrid",
    ) -> List[dict]:
        from pymilvus import AnnSearchRequest, RRFRanker

        t0 = time.time()
        self.store.load_collection(collection_name)

        dense, sparse = self.embedder.encode_dual(dense_text, sparse_text)

        if (sparse is None or not self.sparse_field) and channel != "dense":
            # Backend disabled sparse, or the collection has no sparse column
            # (dense-only schemas like `test`) — degrade gracefully.
            channel = "dense"

        if channel == "hybrid":
            req_dense = AnnSearchRequest(
                data=[dense],
                anns_field=self.dense_field,
                param={"metric_type": self.dense_metric,
                       "params": {"ef": max(64, top_k * 2)}},
                limit=top_k,
                expr=filter_expr or None,
            )
            req_sparse = AnnSearchRequest(
                data=[sparse],
                anns_field=self.sparse_field,
                param={"metric_type": "IP",
                       "params": {"drop_ratio_search": 0.1}},
                limit=top_k,
                expr=filter_expr or None,
            )
            kwargs: dict = dict(
                collection_name=collection_name,
                reqs=[req_dense, req_sparse],
                ranker=RRFRanker(k=self.rrf_k),
                limit=top_k,
                output_fields=output_fields,
            )
            if partition_names:
                kwargs["partition_names"] = partition_names
            results = self.store.client.hybrid_search(**kwargs)
        else:
            anns_field = self.dense_field if channel == "dense" else self.sparse_field
            search_params = (
                {"metric_type": self.dense_metric,
                 "params": {"ef": max(64, top_k * 2)}}
                if channel == "dense"
                else {"metric_type": "IP",
                      "params": {"drop_ratio_search": 0.1}}
            )
            data_vec = [dense] if channel == "dense" else [sparse]
            search_kwargs: dict = dict(
                collection_name=collection_name,
                data=data_vec,
                anns_field=anns_field,
                search_params=search_params,
                limit=top_k,
                filter=filter_expr or "",
                output_fields=output_fields,
            )
            if partition_names:
                search_kwargs["partition_names"] = partition_names
            results = self.store.client.search(**search_kwargs)

        raw = results[0] if results else []
        hits = [self._row_from_hit(h, collection_name) for h in raw]

        emit(
            "retrieve",
            collection=collection_name,
            channel=channel,
            top_k=top_k,
            filter_expr=filter_expr,
            partition_names=partition_names,
            n_hits=len(hits),
            top_pks=[h.get("pk") for h in hits[:top_k]],
            top_doc_ids=[h.get(self.doc_id_field) for h in hits[:top_k]],
            top_parent_ids=[h.get("parent_id") for h in hits[:top_k]],
            top_scores=[round(float(h.get("score", 0.0)), 4) for h in hits[:top_k]],
            latency_ms=round((time.time() - t0) * 1000, 1),
        )
        return hits

    # ---------------- high-level search API --------------------------------

    def search(
        self,
        question: str,
        collection: str,
        *,
        top_k: int = 5,
        doc_id: Optional[str] = None,
        section_types: Optional[Sequence[str]] = None,
        extra_filter: str = "",
        with_rerank: bool = False,
        with_small_to_big: bool = True,
        per_doc_dedup: bool = False,
        output_fields: Optional[List[str]] = None,
        channel: Literal["dense", "sparse", "hybrid"] = "hybrid",
        query_pairs: Optional[Sequence[tuple[str, str]]] = None,
    ) -> List[dict]:
        """One-stop hybrid search: preprocess → hybrid → dedup → rerank → small-to-big.

        Parameters
        ----------
        doc_id : if set, hard-filter to a single document (`doc_id_field == "..."`).
        section_types : optional whitelist for `section_type`.
        with_small_to_big : only takes effect when the collection has a
            `chunk_type` column (P0 hybrid schema); silently disabled otherwise.
        per_doc_dedup : when doc_id is not fixed, keep the top-scoring hit per
            `doc_id_field` before rerank to avoid diversity collapse.
        query_pairs : when provided, MULTI-QUERY mode. Each entry is a
            `(dense_text, sparse_text)` pair — typically produced by
            `LLMQueryRewriter.rewrite_multi(question)` and unpacked by
            the caller. Skips `self._preprocess` entirely, runs one
            `_hybrid_search` per pair, RRF-merges the ranked lists,
            then feeds the merged hits into the usual dedup/rerank/
            small-to-big tail. `None` → legacy single-query behaviour.
        """
        if query_pairs:
            pairs = list(query_pairs)
        else:
            dense_text, sparse_text = (
                self._preprocess(question) if self._preprocess else (question, question)
            )
            pairs = [(dense_text, sparse_text)]

        schema_field_names = self._schema_field_names(collection)
        has_chunk_type = "chunk_type" in schema_field_names
        effective_small_to_big = with_small_to_big and has_chunk_type

        filters: list[str] = []
        if effective_small_to_big:
            filters.append('chunk_type == "child"')
        if section_types and "section_type" in schema_field_names:
            quoted = [f'"{s}"' for s in section_types]
            filters.append(f"section_type in [{', '.join(quoted)}]")
        if doc_id:
            filters.append(f'{self.doc_id_field} == "{doc_id}"')
        if extra_filter:
            filters.append(f"({extra_filter})")
        filter_expr = " and ".join(filters)

        # Section types are row metadata, not physical Milvus partitions.
        # Collections may store all section types in the default partition.
        partition_names = None
        fetch_k = top_k if with_rerank else top_k
        of = output_fields or self.default_output_fields or DEFAULT_OUTPUT_FIELDS

        # Single-query fast path preserves EXACTLY the old behaviour —
        # important because most tests / production callers still go
        # through here without touching query_pairs.
        if len(pairs) == 1:
            dense_text, sparse_text = pairs[0]
            hits = self._hybrid_search(
                collection_name=collection,
                top_k=fetch_k,
                output_fields=of,
                dense_text=dense_text,
                sparse_text=sparse_text,
                filter_expr=filter_expr,
                partition_names=partition_names,
                channel=channel,
            )
        else:
            # Multi-query fan-out. Each pair fetches its own top-fetch_k,
            # then we RRF-merge by primary key across the ranked lists.
            # RRF's normalisation makes cross-query score fusion sane
            # even though absolute scores don't line up.
            #
            # The fan-out is CONCURRENT: without a thread pool, three
            # candidates would triple the wall-clock (Milvus RTT is the
            # dominant cost per pair). With N workers = N pairs, wall-
            # clock collapses to max(latencies). The upper worker bound
            # is `multi_query_max_parallel` (defaults to len(pairs) —
            # Milvus can handle it, and pairs are small ints).
            workers = max(1, min(len(pairs), self.multi_query_max_parallel))

            def _run_pair(pair: tuple[str, str]) -> list[dict]:
                d, s = pair
                return self._hybrid_search(
                    collection_name=collection,
                    top_k=fetch_k,
                    output_fields=of,
                    dense_text=d,
                    sparse_text=s,
                    filter_expr=filter_expr,
                    partition_names=partition_names,
                    channel=channel,
                )

            t_fan = time.time()
            hits_lists: list[list[dict]] = [[]] * len(pairs)
            with ThreadPoolExecutor(
                max_workers=workers, thread_name_prefix="rrf",
            ) as pool:
                futs = {pool.submit(_run_pair, p): i for i, p in enumerate(pairs)}
                for fut in futs:
                    hits_lists[futs[fut]] = fut.result()
            hits = _rrf_merge(hits_lists, k=self.rrf_k)
            emit(
                "multi_query_merge",
                n_pairs=len(pairs),
                n_merged=len(hits),
                fanout_latency_ms=round((time.time() - t_fan) * 1000, 1),
                top_pks=[h.get("pk") for h in hits[:top_k]],
            )

        if per_doc_dedup and not doc_id and hits:
            best: dict = {}
            for h in hits:
                key = h.get(self.doc_id_field, "")
                if key not in best or h["score"] > best[key]["score"]:
                    best[key] = h
            hits = sorted(best.values(), key=lambda x: x["score"], reverse=True)

        if with_rerank and hits and self.reranker is not None:
            try:
                hits = self.reranker.rerank(question, hits, top_n=top_k)
            except Exception as e:
                emit("rerank_error", error=str(e), fallback="rrf")
                hits = hits[:top_k]
        else:
            hits = hits[:top_k]

        if effective_small_to_big and hits:
            hits = self._expand_small_to_big(hits, collection=collection)

        return hits

    # ---------------- small-to-big expansion ------------------------------

    def _expand_small_to_big(self, hits: list[dict], *, collection: str) -> list[dict]:
        parent_ids = sorted({h.get("parent_id", "") for h in hits if h.get("parent_id")})
        if not parent_ids:
            return hits
        self.store.load_collection(collection)
        quoted = [f'"{pid}"' for pid in parent_ids]
        parent_rows = self.store.client.query(
            collection_name=collection,
            filter=f'chunk_type == "parent" and parent_id in [{", ".join(quoted)}]',
            output_fields=["parent_id", self.content_field, "assets"],
            limit=len(parent_ids) * 2,
        )
        parent_map = {r["parent_id"]: r for r in parent_rows}
        for h in hits:
            pid = h.get("parent_id", "")
            if pid in parent_map:
                h["child_content"] = h.get(self.content_field, "")
                h[self.content_field] = parent_map[pid].get(self.content_field, "")
                pa = parent_map[pid].get("assets", "")
                if pa:
                    h["assets"] = pa
        return hits

    # ---------------- full document retrieval ------------------------------

    def get_full_doc(
        self,
        doc_id: str,
        *,
        collection: str,
        section_order: Optional[Sequence[str]] = None,
        limit: int = 100,
    ) -> Optional[dict]:
        """Fetch all parent chunks of one document, sorted by section order.

        On collections without a `chunk_type` column (metadata-only layouts
        like `hypernews_hybrid`), we skip the `chunk_type == "parent"` clause
        instead of tripping a Milvus parse error, and return whatever single
        row matches `doc_id_field == doc_id` (its content_field, if any).
        """
        self.store.load_collection(collection)
        schema_fields = self._schema_field_names(collection)
        parts = [f'{self.doc_id_field} == "{doc_id}"']
        if "chunk_type" in schema_fields:
            parts.insert(0, 'chunk_type == "parent"')
        filter_expr = " and ".join(parts)

        wanted = [
            self.doc_id_field, "title", "authors", "date",
            "section_number", "section_title", "section_type",
            self.content_field,
        ]
        output_fields = [f for f in wanted if f in schema_fields or f == self.doc_id_field]

        rows = self.store.client.query(
            collection_name=collection,
            filter=filter_expr,
            output_fields=output_fields,
            limit=limit,
        )
        if not rows:
            return None

        order = list(section_order) if section_order else []

        def _sort_key(r: dict) -> tuple:
            st = r.get("section_type", "")
            idx = order.index(st) if st in order else 99
            try:
                num = float(r.get("section_number") or 0)
            except (ValueError, TypeError):
                num = 0.0
            return (idx, num)

        rows.sort(key=_sort_key)
        meta = rows[0]
        return {
            "doc_id": doc_id,
            "title": meta.get("title", ""),
            "authors": meta.get("authors", ""),
            "date": meta.get("date", ""),
            "sections": rows,
        }

    def list_indexed_docs(
        self,
        collection: str,
        *,
        limit: int = 500,
        marker_section: str = "abstract",
    ) -> list[dict]:
        """Return one entry per doc_id in the collection (by scanning parent-chunks)."""
        self.store.load_collection(collection)
        rows = self.store.client.query(
            collection_name=collection,
            filter=(
                f'chunk_type == "parent" and '
                f'section_type == "{marker_section}"'
            ),
            output_fields=[self.doc_id_field, "title", "authors", "date"],
            limit=limit,
        )
        seen: set[str] = set()
        result: list[dict] = []
        for r in rows:
            key = r.get(self.doc_id_field, "")
            if key and key not in seen:
                seen.add(key)
                result.append({
                    "doc_id": key,
                    "title": r.get("title", ""),
                    "authors": r.get("authors", ""),
                    "date": r.get("date", ""),
                })
        result.sort(key=lambda x: x["doc_id"])
        return result

    # ---------------- upsert ------------------------------------------------

    def upsert_chunks(
        self,
        collection: str,
        chunks: Iterable[dict],
        *,
        embed_text_key: str = "embed_text",
        content_key: str = "content",
        partition: Optional[str] = None,
    ) -> int:
        """Encode each chunk's `embed_text_key` (or `content_key`) and insert.

        Each chunk dict must supply the schema fields (`doc_id`, `chunk_type`,
        `parent_id`, `section_type`, `section_number`, `section_title`,
        `title`, `authors`, `date`, `source_url`, `content`, `assets`,
        `metadata_json`). Missing fields default to empty strings.
        """
        import json

        chunks = list(chunks)
        if not chunks:
            return 0

        texts = [c.get(embed_text_key) or c.get(content_key) or "" for c in chunks]
        dense_vecs, sparse_vecs = self.embedder.encode_batch(texts)

        def _s(d: dict, k: str, default: str = "") -> str:
            v = d.get(k, default)
            if v is None:
                return default
            if isinstance(v, (dict, list)):
                return json.dumps(v, ensure_ascii=False)
            return str(v)

        rows: list[dict] = []
        for c, dv, sv in zip(chunks, dense_vecs, sparse_vecs):
            rows.append({
                "doc_id": _s(c, "doc_id"),
                "chunk_type": _s(c, "chunk_type", "child"),
                "parent_id": _s(c, "parent_id"),
                "section_type": _s(c, "section_type"),
                "section_number": _s(c, "section_number"),
                "section_title": _s(c, "section_title"),
                "title": _s(c, "title"),
                "authors": _s(c, "authors"),
                "date": _s(c, "date"),
                "source_url": _s(c, "source_url"),
                "content": _s(c, "content"),
                "assets": _s(c, "assets"),
                "metadata_json": _s(c, "metadata_json"),
                "dense_vector": dv,
                "sparse_vector": sv if sv is not None else {},
            })
        insert_kwargs: dict = dict(collection_name=collection, data=rows)
        if partition:
            insert_kwargs["partition_name"] = partition
        self.store.client.insert(**insert_kwargs)
        self.store.client.flush(collection)
        return len(rows)
