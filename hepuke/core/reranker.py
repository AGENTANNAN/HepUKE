"""HTTP reranker client (OpenAI-compatible /rerank endpoint).

Operates on plain hit dicts (must contain a text field, default `content`) —
no llama-index dependency. Each hit gets a `rerank_score` field; the returned
list is sorted by that score descending. Ported from rag-besiii
`rag/reranker.py::CrossEncoderReranker` with retry/backoff.
"""
from __future__ import annotations

import logging
import time
from typing import Any, Iterable, List, Optional

import requests

logger = logging.getLogger(__name__)


class Reranker:
    def __init__(
        self,
        *,
        base_url: str,
        api_key: str,
        model: str,
        top_n: int = 5,
        timeout: int = 60,
        max_retries: int = 3,
    ) -> None:
        self.url = base_url.rstrip("/") + "/rerank"
        self.api_key = api_key
        self.model = model
        self.top_n = top_n
        self.timeout = timeout
        self.max_retries = max_retries

    def rerank(
        self,
        query: str,
        hits: List[dict],
        *,
        text_field: str = "content",
        top_n: Optional[int] = None,
    ) -> List[dict]:
        if not hits:
            return hits
        top_n = top_n or self.top_n
        documents = [str(h.get(text_field) or "") for h in hits]
        payload = {
            "model": self.model,
            "query": query,
            "top_n": top_n,
            "documents": documents,
        }
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
        }
        last_err: Exception | None = None
        for attempt in range(self.max_retries):
            try:
                resp = requests.post(
                    self.url, headers=headers, json=payload, timeout=self.timeout,
                )
                resp.raise_for_status()
                break
            except requests.HTTPError as e:
                last_err = e
                if resp.status_code in (429, 503) and attempt < self.max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                raise
            except requests.RequestException as e:
                last_err = e
                if attempt < self.max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                raise
        data = resp.json()

        results = data.get("results") or data.get("data") or []

        # DEBUG: log what the reranker actually returned so we can tell the
        # difference between "collector logs a stale doc" vs "reranker really
        # returned this doc". Prints the query (clipped), the top-N reranked
        # (doc_id + score), and the input pool's doc_ids in original order.
        try:
            in_ids = [str(h.get("doc_id", "?")) for h in hits]
            out_preview = []
            for item in (results or [])[:5]:
                idx = item.get("index")
                sc = item.get("relevance_score", item.get("score"))
                did = str(hits[idx].get("doc_id", "?")) if isinstance(idx, int) and 0 <= idx < len(hits) else "?"
                out_preview.append(f"{did}({sc})")
            logger.info(
                "[rerank-debug] q=%r  in_ids(%d)=%s  out_top5=%s",
                (query[:80] + "…") if len(query) > 80 else query,
                len(in_ids),
                in_ids[:10],
                out_preview,
            )
        except Exception:  # never let debug logging break rerank
            pass

        if not results:
            logger.warning("rerank response empty, keeping RRF order")
            return hits[:top_n]

        scored: list[dict] = []
        for item in results:
            idx = item.get("index")
            score = item.get("relevance_score", item.get("score"))
            if idx is None or idx >= len(hits):
                continue
            row = dict(hits[idx])
            if score is not None:
                # Promote the cross-encoder relevance to the primary
                # `score` so downstream consumers (agent `_compact`,
                # offline evals) see a value consistent with the reranked
                # ORDER. The original RRF fusion score is preserved under
                # `rrf_score` for diagnostics — before this, `_compact`
                # dropped `rerank_score` and surfaced the stale RRF score,
                # which no longer matched the row's rank (5ade9c9c probe).
                if "score" in row and "rrf_score" not in row:
                    row["rrf_score"] = row["score"]
                row["rerank_score"] = float(score)
                row["score"] = float(score)
            scored.append(row)
        scored.sort(key=lambda x: x.get("rerank_score", 0.0), reverse=True)
        return scored[:top_n] if scored else hits[:top_n]


def build_reranker(cfg) -> Reranker | None:
    if not cfg.enabled:
        return None
    return Reranker(
        base_url=cfg.base_url,
        api_key=cfg.api_key,
        model=cfg.model,
        top_n=cfg.top_n,
    )
