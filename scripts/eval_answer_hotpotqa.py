"""HotpotQA official Ans/Sup evaluation harness (A2).

Runs the RAG agent on the same 500-query dev sample used by the T7 replay
collector, then picks supporting sentences with a *separate* small LLM
(deepseek-v4-flash) so the Sup/Joint F1 numbers are directly comparable
to the official hotpot_evaluate_v1 metric.

Design notes
------------
* Sentence boundaries come from ``storage/hotpot_title2sents.json`` (A1),
  so ``(title, sent_idx)`` predictions align 1:1 with the gold JSON.
* Milvus ``doc_id`` = ``"-".join(title.split())[:200]`` — a slug. We
  reverse via ``storage/hotpot_slug2title.json`` (A1, no collisions).
* Agent runs the SAME loop shape as the collector (no route, no
  subagents, ``hotpotqa_distractor`` collection). The only variable
  we manipulate is ``confidence.enabled`` via ``--ctrl on|off``.
* Answer selection LLM (v4-flash) is *only* used for sp picking; the
  agent still generates the answer with the config-set model.
* Output: one JSONL row per query, ready for ``hotpot_score.py`` (A3).
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Iterable, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root

from hepuke import HepUKE, load_config
from hepuke.agent.confidence.hotpotqa_queries import (
    load_hotpotqa_json,
    sample_examples,
)
from hepuke.agent.defaults import build_default_tools
from hepuke.agent.loop import RagAgent
from hepuke.agent.prompts import build_system_prompt
from hepuke.agent.tools import ToolRegistry
from hepuke.config import ConfidenceConfig
from hepuke.core.llm import LLMClient


DEFAULT_HOTPOT_JSON = "data/HotpotQA/raw/hotpot_dev_distractor_v1.json"
DEFAULT_SLUG2TITLE = "storage/hotpot_slug2title.json"
DEFAULT_TITLE2SENTS = "storage/hotpot_title2sents.json"


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hotpot-json", type=Path, default=Path(DEFAULT_HOTPOT_JSON))
    p.add_argument("--slug2title", type=Path, default=Path(DEFAULT_SLUG2TITLE))
    p.add_argument("--title2sents", type=Path, default=Path(DEFAULT_TITLE2SENTS))
    p.add_argument("--out", required=True, type=Path,
                   help="Output JSONL: one row per qid.")
    p.add_argument("--n", type=int, default=500,
                   help="Sample size (must match the collector run: 500).")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--T-max", type=int, default=8)
    p.add_argument("--top-b", type=int, default=8)

    p.add_argument("--ctrl", choices=("on", "off"), required=True,
                   help="Override cfg.confidence.enabled for this run. "
                        "The single independent variable of the ablation.")

    p.add_argument("--collection", default="hotpotqa_distractor")
    p.add_argument("--model", default=None,
                   help="Override cfg.llm.model for the AGENT run "
                        "(e.g. hepai/deepseek-v4-pro). config.yaml stays "
                        "untouched. Sp picker model is separately fixed.")
    #p.add_argument("--sp-model", default="deepseek-ai/deepseek-v4-flash",
    p.add_argument("--sp-model", default="zhipu/glm-5.3",
                   help="LLM used to pick supporting sentences from "
                        "retrieved candidates. Default v4-flash matches the "
                        "original design (eval header comment). Historical "
                        "note: was temporarily switched to v4-pro because "
                        "flash sometimes returned {} / raw long answer under "
                        "load — that TPM contention has since been fixed.")
    p.add_argument("--short-model", default="zhipu/glm-5.3",
                   help="LLM used to compress the agent's long answer "
                        "to a HotpotQA-style short span. Default v4-flash "
                        "matches --sp-model. Runs AFTER the agent so it does "
                        "not affect controller behaviour.")
    p.add_argument("--posterior-model", default="openai/gpt-4.1-mini",
                   help="LLM used only for the shadow posterior probe "
                        "(cheap top-B logprob harvesting). Independent of "
                        "the main agent model; default gpt-4.1-mini keeps "
                        "the primary loop on deepseek without pushing it to "
                        "return logprobs (which HepAI/deepseek 400s on).")
    p.add_argument("--sp-max-titles", type=int, default=10,
                   help="Cap on distinct titles fed to the sp picker. "
                        "Order preserved from observed_hits (retrieval-rank).")

    p.add_argument("--shard", default=None,
                   help='"i/N" — process only sampled indices where i%%N==shard_i.')
    p.add_argument("--skip-existing", type=Path, nargs="+", default=None,
                   help="One or more previous run JSONLs whose qids to skip. "
                        "Pass THIS run's --out on resume to make it idempotent.")
    p.add_argument("--limit", type=int, default=None,
                   help="Debug knob: stop after this many qids (post-shard).")
    p.add_argument("--workers", type=int, default=1,
                   help="Concurrent queries in-process (threadpool). Each "
                       "worker mints its own RagAgent per query, but shares "
                       "the loaded slug2title / title2sents / rag connection. "
                       "Upstream HepAI rate is the real bottleneck — 4-8 is "
                       "usually the sweet spot.")
    p.add_argument("--score-every", type=int, default=0,
                   help="If >0, run hotpot_score.score_file on --out every "
                        "N completions and print a compact metric line. "
                        "Zero disables in-run scoring.")
    p.add_argument("--auto-resume", action="store_true",
                   help="Treat --out as skip-existing too. Rerunning the "
                        "same command then picks up wherever it left off — "
                        "the canonical way to extend 500 → 7405 later.")
    p.add_argument("--predictor-ckpt", type=Path, default=None,
                   help="Override cfg.confidence.predictor_ckpt for this "
                        "run only. Config.yaml stays untouched — needed "
                        "when running multiple seeds in parallel.")
    p.add_argument("--fusion-ckpt", type=Path, default=None,
                   help="Override cfg.confidence.fusion_ckpt for this run "
                        "only. Config.yaml stays untouched.")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# Agent factory (mirrors collector; controller wired via --ctrl)              #
# --------------------------------------------------------------------------- #

# HotpotQA-friendly sparse extractor. The default preprocess (built for
# BESIII physics) leaves English full-sentence questions untouched, so
# BGE-M3 sparse ends up scoring stopwords like "who", "the", "is" —
# noise that dominates BM25. This extractor:
#   * keeps the dense text intact (BGE-M3 dense handles full questions well)
#   * strips a small English stopword list from sparse text
#   * uppercase-first / all-caps tokens (proper nouns, acronyms, particle
#     symbols) are preserved verbatim and duplicated to raise their sparse
#     weight vs common English words
# It is scoped to this eval script — the module-level `preprocess` used by
# other collections is untouched.
_HOTPOT_STOPWORDS = frozenset("""
a an the and or but if of in on at to for from by with without into onto
is are was were be been being am do does did doing done have has had having
not no nor so than then that this these those it its as at
who whom whose what which where when why how
i you he she we they me him her us them my your his hers our their mine yours
which what where when why how many much more most less least any some all
about above below between through during before after over under again
here there just very really quite too own same other another
can could may might must shall should will would ought
""".split())

_HOTPOT_TOKEN_RE = __import__("re").compile(r"[A-Za-z0-9][A-Za-z0-9\-/]*")


def _hotpot_preprocess(question: str) -> tuple[str, str]:
    """Return ``(dense_text, sparse_text)`` for HotpotQA-style English."""
    dense = question.strip()
    toks = _HOTPOT_TOKEN_RE.findall(dense)
    kept: list[str] = []
    for i, t in enumerate(toks):
        if not t:
            continue
        # Sentence-initial capitalization is not evidence of a proper noun.
        # Check stopword membership case-insensitively first — this drops
        # "Who/What/Where/The/Both/..." while still keeping later Capitalized
        # tokens (real named entities).
        if t.lower() in _HOTPOT_STOPWORDS:
            continue
        # Keep everything else: proper nouns, acronyms, numerals, content
        # words. No further filtering — recall matters more than concision.
        kept.append(t)
    if not kept:
        return dense, dense
    return dense, " ".join(kept)


def _build_agent_factory(cfg, rag: HepUKE, *, T_max: int, top_b: int, ctrl_on: bool, posterior_llm: LLMClient | None = None, predictor_ckpt_override=None, fusion_ckpt_override=None):
    registry = ToolRegistry()
    # HotpotQA only exposes the Milvus collection; grep/read_memo/list_memos
    # etc. point at cfg.wiki.memo_md_root, which is empty for this eval.
    # Leaving them enabled makes the agent burn rounds on tools that can't
    # help — measured 39/59 tool calls in the earlier 20-qid probe were
    # non-search, and every one of the OFF max_steps=8 rows was dominated
    # by grep/read_memo loops. Restrict the surface to search only.
    all_tools = build_default_tools(
        rag.retriever,
        memo_root=cfg.wiki.memo_md_root,
        default_collection=cfg.milvus.default_collection,
        schema_profile=rag.schema_profile,
        use_subagents=cfg.agent.use_subagents,
        llm=rag.llm,
        hypernews_db_path=cfg.agent.hypernews_db_path,
        default_with_rerank=cfg.retrieval.with_reranker,
        rewriter=rag.rewriter,
        # Narrow each search to the single highest-reranked doc so the
        # agent reads one strongly-relevant full doc per lane instead of
        # skimming N shallow hits. Breadth comes from round-2 fan-out.
        search_return_top=1,
    )
    # Keep only `search` and (if wired) `rewrite`. The rewrite tool is the
    # first-round planner; search is single-query retrieval. All other
    # tools point at empty memo/hypernews roots for HotpotQA and would
    # only burn rounds.
    keep = {"search", "rewrite"}
    tools_kept = [t for t in all_tools if t.name in keep]
    if not any(t.name == "search" for t in tools_kept):
        raise RuntimeError("search tool missing from build_default_tools output")
    registry.register_all(tools_kept)
    has_rewrite = any(t.name == "rewrite" for t in tools_kept)
    system_prompt = build_system_prompt(
        rag.schema_profile, use_subagents=cfg.agent.use_subagents,
    )

    conf = ConfidenceConfig(
        enabled=ctrl_on,
        top_b=top_b,
        predictor_ckpt=(str(predictor_ckpt_override) if predictor_ckpt_override
                        else cfg.confidence.predictor_ckpt),
        fusion_ckpt=(str(fusion_ckpt_override) if fusion_ckpt_override
                     else cfg.confidence.fusion_ckpt),
        bocpd_hazard=cfg.confidence.bocpd_hazard,
        r_min=cfg.confidence.r_min,
        z_star=cfg.confidence.z_star,
        epsilon=cfg.confidence.epsilon,
        delta=cfg.confidence.delta,
    )
    # collect_logprobs_only only matters when ctrl is OFF — controller pulls
    # logprobs on its own. Set True either way so posteriors are captured
    # for later diagnostics, at zero controller cost.
    collect_only = not ctrl_on

    def factory() -> RagAgent:
        return RagAgent(
            llm=rag.llm,
            tools=registry,
            max_steps=T_max,
            system_prompt=system_prompt,
            parallel_tool_calls=cfg.agent.parallel_tool_calls,
            max_parallel=cfg.agent.max_parallel,
            route=False,
            confidence=conf,
            collect_logprobs_only=collect_only,
            force_first_tool="rewrite" if has_rewrite else None,
            posterior_llm=posterior_llm,
        )

    return factory


# --------------------------------------------------------------------------- #
# Candidate sentence assembly                                                 #
# --------------------------------------------------------------------------- #
def _unique_doc_ids(observed_hits: list[dict]) -> list[str]:
    """Retrieval-rank-ordered unique doc_ids observed across all rounds.

    First appearance wins so downstream truncation preserves the earliest
    (highest-scoring) titles.
    """
    seen: set[str] = set()
    out: list[str] = []
    for h in observed_hits:
        if not isinstance(h, dict):
            continue
        did = h.get("doc_id")
        if not did or did in seen:
            continue
        seen.add(str(did))
        out.append(str(did))
    return out


def _assemble_candidates(
    doc_ids: list[str],
    slug2title: dict[str, str],
    title2sents: dict[str, list[str]],
    *,
    max_titles: int,
) -> tuple[list[tuple[str, int, str]], list[str]]:
    """Return (labeled_sentences, missing_slugs).

    labeled_sentences: list of ``(title, sent_idx, sentence_text)``, cap at
    ``max_titles`` distinct titles. missing_slugs: doc_ids that didn't map
    to any known title (logged but not fatal — some tools observe non-
    Wikipedia doc_ids).
    """
    labeled: list[tuple[str, int, str]] = []
    missing: list[str] = []
    used_titles: set[str] = set()
    for slug in doc_ids:
        if len(used_titles) >= max_titles:
            break
        title = slug2title.get(slug)
        if title is None:
            missing.append(slug)
            continue
        if title in used_titles:
            continue
        sents = title2sents.get(title)
        if not sents:
            continue
        used_titles.add(title)
        for k, sent in enumerate(sents):
            labeled.append((title, k, sent))
    return labeled, missing


# --------------------------------------------------------------------------- #
# Sp picker (v4-flash)                                                        #
# --------------------------------------------------------------------------- #
_SP_SYSTEM = (
    "You select supporting sentences from Wikipedia snippets that JUSTIFY "
    "a given short answer to a multi-hop question. "
    'Return ONLY a JSON object of the exact shape '
    '{"supporting_facts": [[title, sent_idx], ...]}. '
    "Rules: "
    "(1) titles and sent_idx must appear verbatim in the candidate list; "
    "(2) HotpotQA is multi-hop — the gold supporting-fact set for most "
    "questions spans 2 sentences drawn from 2 DIFFERENT titles (a "
    '"bridge" and a "final" entity). When one candidate title names the '
    "bridge entity and another names the final answer, BOTH belong in "
    "the output. When in doubt, prefer higher recall over higher "
    "precision — a missed required title zeroes out sp-recall, while "
    "one extra sentence only costs some sp-precision. "
    "(3) It is OK (and often required) to include 2 sentences from the "
    "SAME title when both jointly carry the justification, but never "
    "output identical (title, sent_idx) pairs twice. "
    "(4) no prose, no code fences, no comments — the response body must "
    'start with `{` and end with `}`.'
)


# --------------------------------------------------------------------------- #
# Short-answer extractor (v4-flash) — fixes the length mismatch between the   #
# agent's explanatory long answer and HotpotQA's 1-6 token gold span.         #
# Runs AFTER the agent, on the already-generated long answer, so agent        #
# behaviour and controller variables stay unchanged.                          #
# --------------------------------------------------------------------------- #
_SHORT_SYSTEM = (
    "You extract the shortest possible answer span from a long answer to a "
    "question. HotpotQA gold answers are typically 1-6 tokens (a name, "
    'date, number, place, or "yes"/"no"). Rules: (1) output ONLY tokens '
    "that appear in the long answer — no paraphrase; (2) prefer the exact "
    'shortest span that answers the question; (3) if the question is a '
    'yes/no question, output exactly "yes" or "no"; (4) return ONLY a JSON '
    'object of shape {"short": "..."} — no prose, no fences.'
)


def _extract_short_ans(
    short_llm: LLMClient,
    *,
    question: str,
    long_answer: str,
) -> tuple[str, Optional[str]]:
    """Compress a long generative answer down to the HotpotQA gold-length
    span. Returns (short, error). On any failure, short defaults to the
    long answer (so downstream scoring degrades gracefully — worst case
    matches today's numbers)."""
    if not long_answer.strip():
        return "", "empty_long_answer"
    prompt = (
        f"Question: {question}\n"
        f"Long answer: {long_answer}\n\n"
        'Return {"short": "<shortest exact span>"}.'
    )
    try:
        obj = short_llm.ask_json(prompt, system_prompt=_SHORT_SYSTEM)
    except Exception as exc:  # noqa: BLE001
        return long_answer, f"llm_error:{type(exc).__name__}:{exc}"
    val = obj.get("short") if isinstance(obj, dict) else None
    if not isinstance(val, str) or not val.strip():
        return long_answer, "bad_shape"
    return val.strip(), None


def _pick_sp(
    sp_llm: LLMClient,
    *,
    question: str,
    pred_answer: str,
    labeled: list[tuple[str, int, str]],
) -> tuple[list[list], Optional[str]]:
    """Call the sp picker LLM. Returns (supporting_facts, error).

    Robust to model misbehaviour: any invalid entry (unknown title, out-
    of-range idx, malformed shape) is silently dropped rather than crashing
    the batch.
    """
    if not labeled:
        return [], "no_candidates"

    # Compact rendering: `TITLE [k] sentence` per line.
    rendered = "\n".join(
        f"{title} [{k}] {sent}"
        for title, k, sent in labeled
    )
    prompt = (
        f"Question: {question}\n"
        f"Predicted answer: {pred_answer}\n\n"
        "Candidate sentences (one per line, format `TITLE [sent_idx] text`):\n"
        f"{rendered}\n\n"
        "Example — a bridge question whose gold justification spans TWO "
        "sentences from TWO DIFFERENT titles (the common HotpotQA pattern):\n"
        '  Q: "Vasily Agapkin\'s most well-known march was written in '
        'honor of what event?"\n'
        '  A: "the Slavic women accompanying their husbands in the '
        'First Balkan War"\n'
        "  Candidates include:\n"
        "    `Vasily Agapkin [0] Vasily Agapkin ... composed the march "
        "'Farewell of Slavianka' in 1912.`\n"
        "    `Farewell of Slavianka [0] Farewell of Slavianka is a "
        "Russian patriotic march, written by Vasily Agapkin in honour "
        "of the Slavic women accompanying their husbands in the First "
        "Balkan War.`\n"
        '  Correct output: {"supporting_facts": [["Vasily Agapkin", 0], '
        '["Farewell of Slavianka", 0]]}\n'
        "  (Note: BOTH titles are required — the bridge title identifies "
        "who wrote the march, the final title identifies what event it "
        "honored. Returning only one zeroes out sp-recall.)\n\n"
        'Now return supporting_facts as JSON, e.g. {"supporting_facts": '
        '[["Barack Obama", 0], ["Honolulu", 2]]}'
    )

    try:
        obj = sp_llm.ask_json(prompt, system_prompt=_SP_SYSTEM)
    except Exception as exc:  # noqa: BLE001 — one bad row shouldn't kill the run
        return [], f"llm_error:{type(exc).__name__}:{exc}"

    raw = obj.get("supporting_facts") if isinstance(obj, dict) else None
    if not isinstance(raw, list):
        preview = json.dumps(obj)[:200] if isinstance(obj, dict) else repr(obj)[:200]
        return [], f"bad_shape:{preview}"

    # Whitelist against the labeled set so metrics never see hallucinated
    # (title, k) pairs that would just look like precision errors.
    valid: set[tuple[str, int]] = {(t, k) for t, k, _ in labeled}
    out: list[list] = []
    seen: set[tuple[str, int]] = set()
    for row in raw:
        if not isinstance(row, (list, tuple)) or len(row) != 2:
            continue
        t, k = row
        try:
            k = int(k)
        except Exception:
            continue
        if not isinstance(t, str):
            continue
        key = (t, k)
        if key in valid and key not in seen:
            seen.add(key)
            out.append([t, k])
    return out, None


# --------------------------------------------------------------------------- #
# Per-query worker                                                            #
# --------------------------------------------------------------------------- #
def _run_one_query(
    ex,
    *,
    agent_factory,
    sp_llm: LLMClient,
    short_llm: LLMClient,
    slug2title: dict[str, str],
    title2sents: dict[str, list[str]],
    sp_max_titles: int,
    ctrl_on: bool,
    log: logging.Logger,
    round_observer=None,
    session_id: Optional[str] = None,
) -> dict:
    row: dict = {
        "query_id": ex.id,
        "question": ex.question,
        "gold_ans": ex.answer,
        "gold_sp": [list(sp) for sp in ex.supporting_facts],
        "ctrl_on": ctrl_on,
    }
    t0 = time.time()
    result = None
    last_exc: Optional[Exception] = None
    # Retry-once on transient upstream errors (HepAI occasionally returns
    # 400 INTERNAL_ERROR mid-batch). Two attempts covers the observed
    # flake rate; deeper failures land in `row["error"]` for later resume.
    for attempt in (1, 2):
        try:
            agent = agent_factory()
            # collector-side hook: install the TrajectoryCollector as
            # round_observer BEFORE run() so it sees every posterior/
            # round update. `session_id` gets forwarded to agent.run so the
            # SQLite sink can join predictor rows with trace events later.
            if round_observer is not None:
                agent.round_observer = round_observer
            if session_id is not None:
                result = agent.run(ex.question, session_id=session_id)
            else:
                result = agent.run(ex.question)
            last_exc = None
            break
        except Exception as exc:  # noqa: BLE001
            last_exc = exc
            log.warning("  %s attempt %d failed: %s", ex.id, attempt, exc)
            time.sleep(1.0)
    if result is None:
        row["error"] = f"agent:{type(last_exc).__name__}:{last_exc}"
        row["latency_s"] = round(time.time() - t0, 2)
        return row

    row["pred_ans"] = (result.get("answer") or "").strip()
    row["rounds"] = int(result.get("steps") or 0)
    row["terminated"] = result.get("terminated")

    # Precise search-call count from the trace. For HotpotQA the only
    # retrieval tool the agent is allowed to fire is ``search`` (over the
    # ``hotpotqa_distractor`` collection), so this equals "how many Milvus
    # RPCs the agent made for this question" — a cleaner unit than
    # ``rounds`` (which also increments on read_memo / grep turns).
    trace = result.get("trace") or []
    row["search_calls"] = sum(1 for e in trace if e.get("tool") == "search")
    # Break down every tool the agent actually called — used to verify
    # HotpotQA runs stay search-only (no leakage into grep/read_memo/etc).
    tool_counts: dict[str, int] = {}
    for e in trace:
        t = e.get("tool")
        if t:
            tool_counts[t] = tool_counts.get(t, 0) + 1
    row["tool_counts"] = tool_counts
    # Verify search-only collection targeting: any search args_summary
    # mentioning a collection other than hotpotqa_distractor is a bug.
    non_hotpot_search = 0
    for e in trace:
        if e.get("tool") != "search":
            continue
        args = e.get("args_summary") or {}
        coll = args.get("collection") if isinstance(args, dict) else None
        if coll and coll != "hotpotqa_distractor":
            non_hotpot_search += 1
    if non_hotpot_search:
        row["non_hotpot_search"] = non_hotpot_search

    # Short-answer extraction: agent generates explanatory long answers
    # (~60 tokens median) but HotpotQA gold is 1-6 tokens. Without this
    # step, Ans-F1 precision is diluted by ~20× (see hotpot_score.py).
    short_ans, short_err = _extract_short_ans(
        short_llm, question=ex.question, long_answer=row["pred_ans"],
    )
    row["short_ans"] = short_ans
    if short_err:
        row["short_error"] = short_err

    doc_ids = _unique_doc_ids(result.get("observed_hits") or [])
    row["doc_ids"] = doc_ids

    labeled, missing = _assemble_candidates(
        doc_ids, slug2title, title2sents, max_titles=sp_max_titles,
    )
    if missing:
        row["missing_slugs_count"] = len(missing)

    pred_sp, sp_err = _pick_sp(
        sp_llm,
        question=ex.question,
        pred_answer=short_ans or row["pred_ans"],
        labeled=labeled,
    )
    row["pred_sp"] = pred_sp
    if sp_err:
        row["sp_error"] = sp_err

    row["latency_s"] = round(time.time() - t0, 2)
    return row


# --------------------------------------------------------------------------- #
# Skip-existing helper                                                        #
# --------------------------------------------------------------------------- #
def _load_done_qids(paths: Iterable[Path], log: logging.Logger) -> set[str]:
    done: set[str] = set()
    for p in paths:
        if not p.exists():
            log.warning("skip-existing %s missing — ignored", p)
            continue
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    qid = json.loads(line).get("query_id")
                except Exception:
                    continue
                if qid:
                    done.add(str(qid))
    return done


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    log = logging.getLogger("eval_hotpot")

    cfg = load_config()
    cfg.milvus.default_collection = args.collection
    if args.model:
        log.info("agent LLM model override: %s -> %s", cfg.llm.model, args.model)
        cfg.llm.model = args.model
    # Match the collector: hotpot corpus doesn't need memo_search, and
    # subagents would burn the step budget on empty fan-outs.
    cfg.agent.use_subagents = False
    ctrl_on = (args.ctrl == "on")
    log.info("controller = %s (cfg.confidence.enabled override)", ctrl_on)

    slug2title = json.loads(args.slug2title.read_text())
    title2sents = json.loads(args.title2sents.read_text())
    log.info("loaded %d slugs, %d titles", len(slug2title), len(title2sents))

    rag = HepUKE(cfg)
    rag.connect_collection(cfg.milvus.default_collection)

    # Swap the physics preprocessor for a HotpotQA-friendly one. Physics
    # preprocess leaves English questions untouched (dense=sparse=raw
    # question), so BGE-M3 sparse ends up scoring "who/the/is" noise. See
    # _hotpot_preprocess above for details. We only overwrite the private
    # attribute on this run's retriever — cfg.retrieval.with_query_preprocess
    # already turned preprocessing on, so the hook exists.
    if getattr(rag.retriever, "_preprocess", None) is not None:
        rag.retriever._preprocess = _hotpot_preprocess
        log.info("installed hotpot-aware sparse preprocessor")
    else:
        log.warning("retriever._preprocess is None; cfg.retrieval.with_query_preprocess "
                    "may be off. Skipping hotpot preprocessor install.")

    # Sp picker: separate LLMClient so its model override doesn't affect
    # the agent's LLM. Reuse base_url / api_key from cfg.llm (both live on
    # HepAI).
    sp_llm = LLMClient(
        base_url=cfg.llm.base_url,
        api_key=cfg.llm.api_key,
        model=args.sp_model,
        temperature=0.0,
        max_tokens=2048,
    )
    log.info("sp picker LLM = %s (temp=0)", args.sp_model)

    short_llm = LLMClient(
        base_url=cfg.llm.base_url,
        api_key=cfg.llm.api_key,
        model=args.short_model,
        temperature=0.0,
        max_tokens=512,
    )
    log.info("short-ans LLM = %s (temp=0)", args.short_model)

    # Dedicated posterior-probe client. Never used for tool-choice=auto —
    # only for the shadow probe with tool_choice="none" + logprobs=True.
    # Reuses the same base_url / api_key as the main llm so the HepAI/openai
    # gateway routing is identical. `max_tokens` here is only an upper bound;
    # the loop caps the shadow probe at 8 tokens anyway.
    posterior_llm = LLMClient(
        base_url=cfg.llm.base_url,
        api_key=cfg.llm.api_key,
        model=args.posterior_model,
        temperature=0.0,
        max_tokens=8,
    )
    log.info("posterior LLM = %s (shadow probes only)", args.posterior_model)

    examples = load_hotpotqa_json(args.hotpot_json)
    sub = sample_examples(examples, n=args.n, seed=args.seed)
    log.info("sampled %d queries (seed=%d) from %d", len(sub), args.seed, len(examples))

    if args.shard:
        try:
            i_str, n_str = args.shard.split("/")
            shard_i, shard_n = int(i_str), int(n_str)
        except (ValueError, AttributeError):
            log.error("bad --shard %r; expected 'i/N'", args.shard); return 2
        if not (0 <= shard_i < shard_n) or shard_n <= 0:
            log.error("bad --shard indices: i=%d N=%d", shard_i, shard_n); return 2
        sub = [ex for k, ex in enumerate(sub) if k % shard_n == shard_i]
        log.info("shard %d/%d → %d queries", shard_i, shard_n, len(sub))

    if args.skip_existing:
        done = _load_done_qids(args.skip_existing, log)
        before = len(sub)
        sub = [ex for ex in sub if ex.id not in done]
        log.info("skip-existing: %d already covered → %d remaining",
                 before - len(sub), len(sub))

    if args.auto_resume and args.out.exists():
        done = _load_done_qids([args.out], log)
        before = len(sub)
        sub = [ex for ex in sub if ex.id not in done]
        log.info("auto-resume from %s: %d already covered → %d remaining",
                 args.out, before - len(sub), len(sub))

    if args.limit is not None:
        sub = sub[: max(0, int(args.limit))]
        log.info("--limit → %d queries", len(sub))

    factory = _build_agent_factory(
        cfg, rag, T_max=args.T_max, top_b=args.top_b, ctrl_on=ctrl_on,
        posterior_llm=posterior_llm,
        predictor_ckpt_override=args.predictor_ckpt,
        fusion_ckpt_override=args.fusion_ckpt,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    workers = max(1, int(args.workers))
    log.info("workers = %d", workers)

    total = len(sub)
    n_ok = 0
    n_err = 0
    n_sp_err = 0
    # Single writer serialised via a lock — one open file handle, one
    # flush per row, no interleaving. Faster than reopening per row and
    # keeps the JSONL parseable even if the process is killed mid-batch.
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed

    write_lock = threading.Lock()
    done_count = 0
    f_out = args.out.open("a", encoding="utf-8")

    def _handle(row: dict, ex) -> None:
        nonlocal n_ok, n_err, n_sp_err, done_count
        with write_lock:
            f_out.write(json.dumps(row, ensure_ascii=False))
            f_out.write("\n")
            f_out.flush()
            done_count += 1
            i = done_count
            if "error" in row:
                n_err += 1
                log.warning("  [%d/%d] %s FAILED %s (%.1fs)",
                            i, total, ex.id, row["error"], row["latency_s"])
            else:
                n_ok += 1
                if row.get("sp_error"):
                    n_sp_err += 1
                log.info("  [%d/%d] %s rounds=%d docs=%d pred_sp=%d (%.1fs)%s",
                         i, total, ex.id, row["rounds"],
                         len(row["doc_ids"]), len(row["pred_sp"]),
                         row["latency_s"],
                         f" [sp_err={row.get('sp_error')}]" if row.get("sp_error") else "")

            # In-run scoring: reuse hotpot_score.score_file so numbers are
            # produced by the exact same code A6 uses. Held under the lock
            # so the read sees a consistent JSONL snapshot.
            if args.score_every > 0 and i % args.score_every == 0:
                try:
                    from scripts.hotpot_score import score_file as _score  # type: ignore
                except Exception:
                    # scripts/ isn't a package by default — fall back to path load.
                    import importlib.util as _iu
                    spec = _iu.spec_from_file_location(
                        "_hotpot_score",
                        Path(__file__).with_name("hotpot_score.py"),
                    )
                    mod = _iu.module_from_spec(spec)
                    spec.loader.exec_module(mod)  # type: ignore
                    _score = mod.score_file
                res = _score(args.out)
                m = res["metrics"]
                log.info(
                    "  [score @ %d] ans_em=%.2f ans_f1=%.2f sp_em=%.2f "
                    "sp_f1=%.2f joint_em=%.2f joint_f1=%.2f  "
                    "(n_ok=%d n_err=%d rounds=%.2f)",
                    i,
                    100 * m["em"], 100 * m["f1"],
                    100 * m["sp_em"], 100 * m["sp_f1"],
                    100 * m["joint_em"], 100 * m["joint_f1"],
                    res["n_ok"], res["n_error"], res["avg_rounds"],
                )

    def _worker(ex) -> None:
        row = _run_one_query(
            ex,
            agent_factory=factory,
            sp_llm=sp_llm,
            short_llm=short_llm,
            slug2title=slug2title,
            title2sents=title2sents,
            sp_max_titles=args.sp_max_titles,
            ctrl_on=ctrl_on,
            log=log,
        )
        _handle(row, ex)

    try:
        if workers == 1:
            for ex in sub:
                _worker(ex)
        else:
            with ThreadPoolExecutor(
                max_workers=workers, thread_name_prefix="eval",
            ) as pool:
                futures = [pool.submit(_worker, ex) for ex in sub]
                for _ in as_completed(futures):
                    pass  # completion side effects handled inside _worker
    finally:
        f_out.close()

    log.info("done: %d ok / %d err / %d sp_err → %s", n_ok, n_err, n_sp_err, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
