"""Config loader: YAML + ${ENV_VAR} expansion + HEPUKE__* env overrides."""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, Field

try:  # dotenv is required at runtime but keep the import guarded so importing
      # this module never crashes if the package isn't installed yet.
    from dotenv import load_dotenv as _load_dotenv
except Exception:  # pragma: no cover
    def _load_dotenv(*args, **kwargs):  # type: ignore
        return False


_ENV_PATTERN = re.compile(r"\$\{([^}:]+)(?::-([^}]*))?\}")


def _expand_env(value: Any) -> Any:
    if isinstance(value, str):
        def repl(m: re.Match) -> str:
            var, default = m.group(1), m.group(2)
            return os.environ.get(var, default if default is not None else "")
        return _ENV_PATTERN.sub(repl, value)
    if isinstance(value, dict):
        return {k: _expand_env(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand_env(v) for v in value]
    return value


def _apply_env_overrides(data: dict[str, Any]) -> dict[str, Any]:
    """Env like HEPUKE__MILVUS__URI overrides data['milvus']['uri']."""
    for key, val in os.environ.items():
        if not key.startswith("HEPUKE__"):
            continue
        parts = [p.lower() for p in key[len("HEPUKE__"):].split("__")]
        cursor = data
        for p in parts[:-1]:
            cursor = cursor.setdefault(p, {})
        cursor[parts[-1]] = val
    return data


def _migrate_legacy_provider(data: dict[str, Any]) -> dict[str, Any]:
    """Rewrite legacy embedding.provider values to the new Literal set.

    Historically we auto-upgraded ``hepai/bge-m3`` → ``hepai/bge-m3-hybrid``
    because the -hybrid worker returned both dense + sparse. As of
    2026-09 the hybrid worker started returning 500 INTERNAL_ERROR while
    the plain ``hepai/bge-m3`` endpoint stays up, so the auto-rewrite is
    disabled — whatever the user writes in ``config.yaml`` is what we
    call.
    """
    emb = data.get("embedding")
    if isinstance(emb, dict):
        p = str(emb.get("provider", "")).strip().lower()
        if p in ("", "hepai", "hepai_api", "bge_m3", "bge-m3"):
            emb["provider"] = "bge_m3_api"
        elif p in ("local", "bge_m3_local", "gpu"):
            emb["provider"] = "bge_m3_local"
    return data


class MilvusConfig(BaseModel):
    uri: str
    token: str = ""
    database: str = "default"
    default_collection: str = "hepuke_default"


class EmbeddingConfig(BaseModel):
    provider: Literal["bge_m3_api", "bge_m3_local"] = "bge_m3_api"
    base_url: str
    api_key: str
    model: str = "hepai/bge-m3:latest"
    dim: int = 1024
    local_model_path: str | None = None
    enable_sparse: bool = True


class RerankerConfig(BaseModel):
    enabled: bool = False
    provider: str = "hepai"
    base_url: str = ""
    api_key: str = ""
    model: str = ""
    top_n: int = 5


class LLMConfig(BaseModel):
    provider: str = "hepai"
    base_url: str
    api_key: str
    model: str
    temperature: float = 0.1
    max_tokens: int = 1024


class RetrievalConfig(BaseModel):
    chunk_size: int = 500
    chunk_overlap: int = 0
    similarity_top_k: int = 5
    with_reranker: bool = False
    text_search_mode: str = Field(default="like", pattern="^(like|bm25)$")
    hybrid_enabled: bool = True
    rrf_k: int = 60
    with_small_to_big: bool = True
    with_query_preprocess: bool = True
    per_doc_dedup: bool = False
    # LLM query rewriter (root solution for the CN-question / EN-corpus
    # gap and for term-family recall). When enabled, retrievers can call
    # `LLMQueryRewriter.rewrite_multi(question)` to produce multiple
    # English candidate queries, fan out one hybrid search per candidate,
    # and RRF-merge the results. Falls back to the classic dictionary-
    # based Preprocessor on LLM failure — never worse than today.
    query_rewrite_llm_enabled: bool = False
    query_rewrite_n_candidates: int = 4
    query_rewrite_cache_size: int = 512
    query_rewrite_llm_model: str = ""
    query_rewrite_keep_original: bool = True
    query_rewrite_prompt_style: str = "physics"


class AgentConfig(BaseModel):
    enabled: bool = False
    protocol: Literal["openai_tools"] = "openai_tools"
    max_steps: int = 16
    tool_timeout_s: int = 30
    use_subagents: bool = False
    hypernews_db_path: str = ""
    parallel_tool_calls: bool = True
    max_parallel: int = 8
    route: bool = False
    decompose_max_parallel: int = 4
    # Warn the LLM when it re-issues the same (tool_name, args) more than N
    # times in one run — the observed failure mode is a subagent that keeps
    # re-searching a corpus it already knows is empty. 0 disables.
    repeated_tool_call_threshold: int = 3
    # Cap total output (completion) tokens spent inside one RagAgent.run().
    # Once exceeded, the next step is skipped and the loop returns a
    # `terminated=token_budget` result. Backstop against runaway subagent
    # fan-outs — cheaper than max_steps because it counts actual $ spent.
    # Env overrides: HEPUKE_TOKEN_BUDGET / HEPUKE_TOKEN_BUDGET_ENABLED.
    token_budget: int = 20000
    token_budget_enabled: bool = True
    # Per-turn progress summary: when True, at every step-review pause an
    # extra tool-less LLM probe fires to summarise "what we have / what's
    # missing / next" in 1–2 sentences. Surfaced on the review payload as
    # `progress`. Off saves one LLM hop per step at the cost of the human
    # only seeing subagent answer previews.
    enable_perturn_summary: bool = True


class TraceConfig(BaseModel):
    enabled: bool = False
    backend: Literal["jsonl", "sqlite", "both"] = "sqlite"
    file: str = "./trace/rag_trace.jsonl"
    db_path: str = "./storage/traces.db"
    include_context: bool = True
    console: bool = False


class WikiConfig(BaseModel):
    enabled: bool = False
    root: str = ""
    memo_md_root: str = ""
    memo_pdf_root: str = ""
    # Navigation aids for the memo corpus. Both live one level up from
    # memo_md_root by convention. INDEX.md is a single table of every doc's
    # id/title/abstract-lede; summaries/<id>.md is a ~1KB abstract + heading
    # list. These are the corpus's "table of contents" — cheap first hops
    # for the memo_search subagent before it grep-scans the 1150 body files.
    memo_index_path: str = ""
    memo_summaries_root: str = ""
    pages_collection: str = "wiki_pages"
    edges_db: str = "./storage/wiki_edges.sqlite"
    hops: int = 1


class ConfidenceConfig(BaseModel):
    """Controller-based confidence for RAG loop (see entropy/section_method.tex §4).

    When ``enabled=False`` the RAG loop keeps its fixed ``max_steps`` behavior;
    all other fields are inert until the controller is wired in.
    """

    enabled: bool = False
    # Look-ahead: KL-gain predictor checkpoint path (MLP + isotonic calibrator).
    predictor_ckpt: str = ""
    # Look-back: BOCPD geometric hazard 1/λ_h.
    bocpd_hazard: float = 0.1
    # Cold-start: force RETRIEVE for the first r_min rounds.
    r_min: int = 2
    # Answer posterior: top-B logprobs used to build p_t.
    top_b: int = 8
    # Fusion classifier parameters (α, b) checkpoint.
    fusion_ckpt: str = ""
    # Fusion operating threshold z* (controller-stop when z_t > z_star).
    z_star: float = 0.5
    # When set, overrides the fusion ckpt's baked-in z_star (ablation/sweep).
    # None = honour the ckpt's operating point (normal behaviour).
    z_star_override: float | None = None
    # Look-ahead threshold ε for the stand-alone stopping rule (used when
    # fusion_ckpt is absent).
    epsilon: float = 0.01
    # Look-back threshold δ for the stand-alone triggering rule (used when
    # fusion_ckpt is absent).
    delta: float = 0.8
    # Denominator used to normalise the `exhaustion` φ feature: 1 - n_new/denom.
    # Should match the retrieval top_k at collection AND inference time.
    expected_docs_per_round: int = 5


class ServiceConfig(BaseModel):
    host: str = "0.0.0.0"
    port: int = 42796
    route_prefix: str = "/apiv2"
    api_key: str = ""
    persist_dir: str = "./storage"


class Config(BaseModel):
    version: int = 1
    milvus: MilvusConfig
    embedding: EmbeddingConfig
    reranker: RerankerConfig = RerankerConfig()
    llm: LLMConfig
    retrieval: RetrievalConfig = RetrievalConfig()
    service: ServiceConfig = ServiceConfig()
    agent: AgentConfig = AgentConfig()
    trace: TraceConfig = TraceConfig()
    wiki: WikiConfig = WikiConfig()
    confidence: ConfidenceConfig = ConfidenceConfig()


def load_config(path: str | Path | None = None) -> Config:
    """Load config from YAML, expand ${VARS}, apply HEPUKE__* env overrides.

    Reads .env from the CWD (and parents) via python-dotenv before resolving,
    so ${HEPAI_API_KEY} etc. work out of the box.
    """
    _load_dotenv(override=False)

    if path is None:
        env_path = os.environ.get("HEPUKE_CONFIG")
        path = env_path or Path(__file__).resolve().parent.parent / "config.yaml"
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}

    raw = _expand_env(raw)
    raw = _apply_env_overrides(raw)
    raw = _migrate_legacy_provider(raw)
    return Config.model_validate(raw)
