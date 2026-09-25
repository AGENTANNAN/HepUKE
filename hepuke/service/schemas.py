"""Request / response models for the FastAPI service."""
from __future__ import annotations

from typing import Any, List

from pydantic import BaseModel, Field


class DocumentIn(BaseModel):
    text: str
    metadata: dict[str, Any] = Field(default_factory=dict)
    doc_id: str | None = None


class InsertRequest(BaseModel):
    doc: DocumentIn
    split_chunk: bool = False


class InsertResponse(BaseModel):
    msg: str
    doc_id: str | None


class RetrieveRequest(BaseModel):
    content: str
    similarity_top_k: int | None = None
    filters: List[dict] | None = None
    condition: str = "and"
    with_reranker: bool | None = None


class RetrievedNode(BaseModel):
    score: float | None
    text: str
    metadata: dict[str, Any]


class RetrieveResponse(BaseModel):
    nodes: List[RetrievedNode]


class QueryRequest(RetrieveRequest):
    pass


class QueryResponse(BaseModel):
    answer: str
    sources: List[RetrievedNode]


class DocLookupRequest(BaseModel):
    key: str
    value: Any


class SimpleResponse(BaseModel):
    ok: bool
    detail: Any = None
