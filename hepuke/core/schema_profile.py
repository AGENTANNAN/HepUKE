"""Auto-detect collection schema so `HepUKE` can talk to heterogeneous Milvus
collections without hard-coding field names.

We inspect the collection's schema (via `MilvusClient.describe_collection`)
once at `connect_collection()` and derive:

- `doc_id_field`  : the per-document identifier (used in `doc_id == "..."`
                    filters, dedup, small-to-big linking).
- `content_fields`: the VARCHAR field(s) that hold retrievable text. Multiple
                    are supported (e.g. `question_text` + `answer_text`).
- `dense_field`   : the FLOAT_VECTOR field.
- `sparse_field`  : the SPARSE_FLOAT_VECTOR field, or None (dense-only).
- `dense_metric`  : metric_type on the dense index (from `describe_index`).
- `has_chunking`  : True iff the collection has both `chunk_type` and
                    `parent_id` (parent-child layout).
- `has_section_type`: True iff `section_type` exists.
- `output_fields` : all non-vector fields (safe default for search output).

Milvus DataType numeric codes we care about (per pymilvus):

    5   = INT64          (usually the auto_id pk `pk`)
    21  = VARCHAR
    101 = FLOAT_VECTOR
    104 = SPARSE_FLOAT_VECTOR
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional


# Field names we treat as "the primary retrievable text" in decreasing priority.
_CONTENT_PREFERENCES = ("content", "text", "chunk_text", "body")

# Extra text fields we'll happily include when nothing above matches.
_SECONDARY_TEXT_HINTS = ("question_text", "answer_text", "title", "summary")

# Names we treat as the "document id" in decreasing priority.
_DOC_ID_PREFERENCES = ("doc_id", "bam_id", "paper_id", "document_id")

_DTYPE_INT64 = 5
_DTYPE_VARCHAR = 21
_DTYPE_FLOAT_VECTOR = 101
_DTYPE_SPARSE_FLOAT_VECTOR = 104


@dataclass
class SchemaProfile:
    """Result of `detect_schema_profile()`. Frozen at connect time."""
    collection: str
    doc_id_field: str
    content_fields: List[str]
    dense_field: Optional[str]
    sparse_field: Optional[str]
    dense_metric: str
    has_chunking: bool
    has_section_type: bool
    has_parent_id: bool
    output_fields: List[str] = field(default_factory=list)

    @property
    def primary_content_field(self) -> str:
        """First content field, used as the canonical text column."""
        return self.content_fields[0] if self.content_fields else "content"


def _dtype_of(field_desc: dict) -> int:
    """`describe_collection` returns pymilvus `DataType` enums; coerce to int."""
    dt = field_desc.get("type")
    return int(dt) if dt is not None else 0


def _detect_dense_metric(store, collection: str, dense_field: Optional[str]) -> str:
    """Look up the metric_type actually configured on the dense vector index.

    Falls back to COSINE (P0 default) when no index info is available.
    """
    if not dense_field:
        return "COSINE"
    try:
        info = store.client.describe_index(
            collection_name=collection, index_name=dense_field
        )
        metric = info.get("metric_type") if isinstance(info, dict) else None
        if metric:
            return str(metric).upper()
    except Exception:
        pass
    return "COSINE"


def detect_schema_profile(store, name: str) -> SchemaProfile:
    """Read the collection schema via `MilvusClient` and derive a SchemaProfile."""
    desc = store.describe_collection(name)
    fields = list(desc.get("fields", []))

    field_names = [f["name"] for f in fields]
    field_by_name = {f["name"]: f for f in fields}

    # -- vector fields ------------------------------------------------------
    dense_field: Optional[str] = None
    sparse_field: Optional[str] = None
    for f in fields:
        dt = _dtype_of(f)
        if dt == _DTYPE_FLOAT_VECTOR and dense_field is None:
            dense_field = f["name"]
        elif dt == _DTYPE_SPARSE_FLOAT_VECTOR and sparse_field is None:
            sparse_field = f["name"]

    # -- doc_id field -------------------------------------------------------
    doc_id_field = ""
    for candidate in _DOC_ID_PREFERENCES:
        if candidate in field_by_name:
            doc_id_field = candidate
            break
    if not doc_id_field:
        # Fall back to primary key. If it's auto_id INT64 (`pk`) it's not
        # really usable as a doc id, but at least filter expressions won't
        # reference a non-existent column.
        for f in fields:
            if f.get("is_primary"):
                doc_id_field = f["name"]
                break

    # -- content field(s) ---------------------------------------------------
    content_fields: List[str] = []
    for candidate in _CONTENT_PREFERENCES:
        f = field_by_name.get(candidate)
        if f is not None and _dtype_of(f) == _DTYPE_VARCHAR:
            content_fields.append(candidate)
    if not content_fields:
        excluded = {doc_id_field, "pk", "id", "assets", "metadata_json"}
        candidates = [
            f for f in fields
            if _dtype_of(f) == _DTYPE_VARCHAR
            and f["name"] not in excluded
            and not f.get("is_primary")
            and not f["name"].endswith("_url")
            and not f["name"].endswith("_id")
        ]
        cand_names = {c["name"] for c in candidates}
        for hint in _SECONDARY_TEXT_HINTS:
            if hint in cand_names:
                content_fields.append(hint)
        remaining = [c for c in candidates if c["name"] not in content_fields]
        remaining.sort(
            key=lambda f: (f.get("params") or {}).get("max_length", 0),
            reverse=True,
        )
        for c in remaining[:2]:
            content_fields.append(c["name"])

    if not content_fields:
        content_fields = ["content"]

    has_chunking = "chunk_type" in field_by_name and "parent_id" in field_by_name
    has_section_type = "section_type" in field_by_name
    has_parent_id = "parent_id" in field_by_name

    output_fields = [
        f["name"] for f in fields
        if _dtype_of(f) not in (_DTYPE_FLOAT_VECTOR, _DTYPE_SPARSE_FLOAT_VECTOR)
    ]

    return SchemaProfile(
        collection=name,
        doc_id_field=doc_id_field or "doc_id",
        content_fields=content_fields,
        dense_field=dense_field,
        sparse_field=sparse_field,
        dense_metric=_detect_dense_metric(store, name, dense_field),
        has_chunking=has_chunking,
        has_section_type=has_section_type,
        has_parent_id=has_parent_id,
        output_fields=output_fields,
    )
