"""Milvus wrapper: pure `MilvusClient` API (no ORM).

Previously this module mixed `MilvusClient` (admin plane) with pymilvus ORM
(`connections.connect` + `Collection`) for the retrieval plane. Under
pymilvus 2.6 the ORM path is flagged with `PyMilvusDeprecationWarning` and
scheduled for removal in 3.1, so the retrieval plane was migrated to
`MilvusClient` as well.

The generic parent-child schema is still declared via `FieldSchema` /
`CollectionSchema` (those live outside the ORM connection code and are the
only sane way to spell out a heterogeneous schema before creation). Index
setup uses `MilvusClient.prepare_index_params` / `create_index`.
"""
from __future__ import annotations

import logging
from typing import Any, List

from pymilvus import CollectionSchema, DataType, FieldSchema, MilvusClient

logger = logging.getLogger(__name__)


# ---- generic parent-child schema (shared by P1 papers/memos and P2 wiki) ----

def _default_fields(dim: int) -> list:
    """Return the shared parent-child schema field list.

    Fields intentionally use generic names (`doc_id` / `section_type` / ...).
    Callers stuff their own domain-specific extras into `metadata_json`.
    """
    return [
        FieldSchema(name="pk", dtype=DataType.INT64, is_primary=True, auto_id=True),
        FieldSchema(name="doc_id", dtype=DataType.VARCHAR, max_length=128),
        FieldSchema(name="chunk_type", dtype=DataType.VARCHAR, max_length=16),
        FieldSchema(name="parent_id", dtype=DataType.VARCHAR, max_length=128),
        FieldSchema(name="section_type", dtype=DataType.VARCHAR, max_length=64),
        FieldSchema(name="section_number", dtype=DataType.VARCHAR, max_length=16),
        FieldSchema(name="section_title", dtype=DataType.VARCHAR, max_length=256),
        FieldSchema(name="title", dtype=DataType.VARCHAR, max_length=512),
        FieldSchema(name="authors", dtype=DataType.VARCHAR, max_length=1024),
        FieldSchema(name="date", dtype=DataType.VARCHAR, max_length=64),
        FieldSchema(name="source_url", dtype=DataType.VARCHAR, max_length=256),
        FieldSchema(name="content", dtype=DataType.VARCHAR, max_length=65535),
        FieldSchema(name="assets", dtype=DataType.VARCHAR, max_length=65535),
        FieldSchema(name="metadata_json", dtype=DataType.VARCHAR, max_length=65535),
        FieldSchema(name="dense_vector", dtype=DataType.FLOAT_VECTOR, dim=dim),
        FieldSchema(name="sparse_vector", dtype=DataType.SPARSE_FLOAT_VECTOR),
    ]


_DEFAULT_SCALAR_INDEX_FIELDS = ("doc_id", "chunk_type", "parent_id", "section_type")


class MilvusStore:
    def __init__(
        self,
        *,
        uri: str,
        token: str,
        dim: int,
        database: str = "default",
        text_search_mode: str = "like",
    ) -> None:
        self.uri = uri
        self.token = token
        self.database = database
        self.dim = dim
        self.text_search_mode = text_search_mode
        self._client = MilvusClient(uri=uri, token=token, db_name=database)
        # Track collections we've already loaded so we don't call
        # `load_collection` on every search.
        self._loaded: set[str] = set()

    # ---- client ------------------------------------------------------------
    @property
    def client(self) -> MilvusClient:
        return self._client

    def load_collection(self, name: str) -> None:
        """Idempotent load. `MilvusClient.load_collection` is a no-op when the
        collection is already loaded, but we still cache to avoid the RPC."""
        if name in self._loaded:
            return
        self._client.load_collection(name)
        self._loaded.add(name)

    def describe_collection(self, name: str) -> dict:
        return self._client.describe_collection(name)

    # ---- admin plane -------------------------------------------------------
    def list_collections(self) -> List[str]:
        return list(self._client.list_collections())

    def has_collection(self, name: str) -> bool:
        return self._client.has_collection(name)

    def drop_collection(self, name: str) -> None:
        if self._client.has_collection(name):
            self._client.drop_collection(name)
            self._loaded.discard(name)
            logger.info("dropped collection %s", name)

    def query_by_expr(
        self,
        collection: str,
        expr: str,
        limit: int = 1000,
        output_fields: List[str] | None = None,
    ) -> List[dict]:
        self.load_collection(collection)
        return list(
            self._client.query(
                collection_name=collection,
                filter=expr,
                limit=limit,
                output_fields=output_fields or ["*"],
            )
        )

    def get_by_ids(
        self,
        collection: str,
        ids: List[str],
        output_fields: List[str] | None = None,
    ) -> List[dict]:
        self.load_collection(collection)
        return list(
            self._client.get(
                collection_name=collection,
                ids=ids,
                output_fields=output_fields or ["*"],
            )
        )

    def stats(self, collection: str) -> dict[str, Any]:
        return self._client.get_collection_stats(collection)

    def build_text_expr(self, keyword: str) -> str:
        safe = keyword.replace('"', '\\"')
        if self.text_search_mode == "bm25":
            return f'TEXT_MATCH(content, "{safe}")'
        return f'content like "%{safe}%"'

    # ---- creation ----------------------------------------------------------
    def create_hybrid_collection(
        self,
        name: str,
        *,
        fields: list | None = None,
        partitions: list[str] | None = None,
        scalar_index_fields: list[str] | tuple[str, ...] = _DEFAULT_SCALAR_INDEX_FIELDS,
        description: str = "",
    ) -> None:
        """Create a hybrid (dense+sparse) collection with the generic schema.

        Indexes:
          - dense_vector : HNSW, COSINE (M=16, efConstruction=64)
          - sparse_vector: SPARSE_INVERTED_INDEX, IP (drop_ratio_build=0.2)
          - scalar: INVERTED on each of `scalar_index_fields`
        """
        if self._client.has_collection(name):
            logger.info("collection %s already exists", name)
            return

        schema = CollectionSchema(
            fields=fields or _default_fields(self.dim),
            description=description or f"hepuke hybrid collection: {name}",
            enable_dynamic_field=False,
        )
        self._client.create_collection(collection_name=name, schema=schema)

        index_params = self._client.prepare_index_params()
        index_params.add_index(
            field_name="dense_vector",
            index_type="HNSW",
            index_name="dense_vector",
            metric_type="COSINE",
            params={"M": 16, "efConstruction": 64},
        )
        index_params.add_index(
            field_name="sparse_vector",
            index_type="SPARSE_INVERTED_INDEX",
            index_name="sparse_vector",
            metric_type="IP",
            params={"drop_ratio_build": 0.2},
        )
        for f in scalar_index_fields:
            try:
                index_params.add_index(field_name=f, index_type="INVERTED")
            except Exception as e:  # pragma: no cover
                logger.warning("scalar INVERTED index on %s failed: %s", f, e)
        self._client.create_index(collection_name=name, index_params=index_params)

        for p in partitions or []:
            if not self._client.has_partition(name, p):
                self._client.create_partition(name, p)

        logger.info(
            "created hybrid collection %s (dim=%s, partitions=%s)",
            name, self.dim, partitions or [],
        )
