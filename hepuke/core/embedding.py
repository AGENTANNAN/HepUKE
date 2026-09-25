"""BGE-M3 embedding backends (dense + sparse) for Milvus hybrid retrieval.

Default: `BgeM3ApiBackend` (HepAI `hepai/bge-m3`, no local GPU needed).
Optional (via `pip install hepuke[gpu]`): `BgeM3LocalBackend` (FlagEmbedding + torch).

The Protocol / factory here is the contract that P1 HybridRetriever consumes.
"""
from __future__ import annotations

import threading
from typing import Any, Iterable, Protocol, Tuple, runtime_checkable


SparseVec = dict[int, float]
DenseVec = list[float]


@runtime_checkable
class EmbeddingBackend(Protocol):
    dim: int
    supports_sparse: bool

    def encode(self, text: str) -> Tuple[DenseVec, SparseVec | None]: ...

    def encode_batch(
        self, texts: list[str]
    ) -> Tuple[list[DenseVec], list[SparseVec | None]]: ...

    def encode_dual(
        self, dense_text: str, sparse_text: str
    ) -> Tuple[DenseVec, SparseVec | None]: ...


def _to_int_sparse(raw: dict) -> SparseVec:
    return {int(k): float(v) for k, v in raw.items()}


class BgeM3ApiBackend:
    """HepAI `hepai/bge-m3-hybrid` remote backend (default).

    Ported from rag-besiii `EmbeddingModelAPI` — response shape:
      result['data'][i] = { 'embedding': [1024 floats], 'sparse_embedding': {str: float} }
    """

    supports_sparse = True

    def __init__(
        self,
        *,
        base_url: str,
        api_key: str,
        model: str = "hepai/bge-m3:latest",
        dim: int = 1024,
        enable_sparse: bool = True,
    ) -> None:
        try:
            from hepai import HRModel  # type: ignore
        except ImportError as e:  # pragma: no cover
            raise RuntimeError(
                "hepai package required for BgeM3ApiBackend. "
                "Install with `pip install hepai`."
            ) from e
        self._model = HRModel.connect(name=model, base_url=base_url, api_key=api_key)
        self._model_name = model
        self.dim = dim
        self.enable_sparse = enable_sparse
        self.supports_sparse = enable_sparse

    def encode(self, text: str) -> Tuple[DenseVec, SparseVec | None]:
        result = self._model.embeddings(input=text, model=self._model_name)
        item = result["data"][0]
        sparse = None
        if self.enable_sparse and "sparse_embedding" in item:
            sparse = _to_int_sparse(item["sparse_embedding"])
        return item["embedding"], sparse

    def encode_batch(
        self, texts: list[str]
    ) -> Tuple[list[DenseVec], list[SparseVec | None]]:
        result = self._model.embeddings(input=texts, model=self._model_name)
        items = sorted(result["data"], key=lambda x: x["index"])
        dense = [it["embedding"] for it in items]
        if self.enable_sparse and items and "sparse_embedding" in items[0]:
            sparse: list[SparseVec | None] = [
                _to_int_sparse(it["sparse_embedding"]) for it in items
            ]
        else:
            sparse = [None] * len(items)
        return dense, sparse

    def encode_dual(
        self, dense_text: str, sparse_text: str
    ) -> Tuple[DenseVec, SparseVec | None]:
        # HepAI API doesn't accept two texts per call; fall back to dense_text.
        return self.encode(dense_text)


class BgeM3LocalBackend:
    """Local BGE-M3 backend (requires `pip install hepuke[gpu]` + weights)."""

    supports_sparse = True

    def __init__(
        self,
        *,
        model_path: str,
        dim: int = 1024,
        enable_sparse: bool = True,
        device: str = "cuda",
        use_fp16: bool = True,
    ) -> None:
        try:
            from FlagEmbedding import BGEM3FlagModel  # type: ignore
        except ImportError as e:
            raise RuntimeError(
                "FlagEmbedding not installed. `pip install hepuke[gpu]` "
                "or install FlagEmbedding + torch manually."
            ) from e
        self._model = BGEM3FlagModel(model_path, use_fp16=use_fp16, device=device)
        # BGEM3FlagModel.encode is not thread-safe: it shares tokenizer + torch
        # buffers, and concurrent forwards from parallel tool dispatch have
        # been observed to segfault libcuda on driver 580.95.05 + Blackwell
        # (see storage/collector_smoke_w2.log). Serialize forwards here.
        self._encode_lock = threading.Lock()
        self.dim = dim
        self.enable_sparse = enable_sparse
        self.supports_sparse = enable_sparse

    def encode(self, text: str) -> Tuple[DenseVec, SparseVec | None]:
        with self._encode_lock:
            out = self._model.encode(text, return_dense=True, return_sparse=self.enable_sparse)
        dense = out["dense_vecs"].tolist()
        sparse = (
            _to_int_sparse(out["lexical_weights"]) if self.enable_sparse else None
        )
        return dense, sparse

    def encode_batch(
        self, texts: list[str]
    ) -> Tuple[list[DenseVec], list[SparseVec | None]]:
        with self._encode_lock:
            out = self._model.encode(
                texts,
                return_dense=True,
                return_sparse=self.enable_sparse,
                batch_size=32,
            )
        dense = [v.tolist() for v in out["dense_vecs"]]
        if self.enable_sparse:
            sparse: list[SparseVec | None] = [
                _to_int_sparse(w) for w in out["lexical_weights"]
            ]
        else:
            sparse = [None] * len(texts)
        return dense, sparse

    def encode_dual(
        self, dense_text: str, sparse_text: str
    ) -> Tuple[DenseVec, SparseVec | None]:
        if dense_text == sparse_text:
            return self.encode(dense_text)
        with self._encode_lock:
            dense_out = self._model.encode(dense_text, return_dense=True, return_sparse=False)
            if not self.enable_sparse:
                return dense_out["dense_vecs"].tolist(), None
            sparse_out = self._model.encode(sparse_text, return_dense=False, return_sparse=True)
        return dense_out["dense_vecs"].tolist(), _to_int_sparse(sparse_out["lexical_weights"])


def build_embedding(cfg) -> EmbeddingBackend:
    """Factory: EmbeddingConfig -> EmbeddingBackend."""
    provider = cfg.provider
    if provider == "bge_m3_api":
        return BgeM3ApiBackend(
            base_url=cfg.base_url,
            api_key=cfg.api_key,
            model=cfg.model or "hepai/bge-m3:latest",
            dim=cfg.dim,
            enable_sparse=cfg.enable_sparse,
        )
    if provider == "bge_m3_local":
        if not cfg.local_model_path:
            raise ValueError(
                "EmbeddingConfig.local_model_path must be set for provider=bge_m3_local"
            )
        return BgeM3LocalBackend(
            model_path=cfg.local_model_path,
            dim=cfg.dim,
            enable_sparse=cfg.enable_sparse,
        )
    raise ValueError(f"Unknown embedding provider: {provider!r}")
