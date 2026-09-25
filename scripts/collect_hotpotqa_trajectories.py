"""Smoke-run T7 replay collector on a small HotpotQA slice.

Usage:
    # 1) download the distractor dev split once (~50 MB) into `data/`:
    #      curl -L -o data/hotpot_dev_distractor_v1.json \
    #        https://hotpotqa.github.io/data/hotpot_dev_distractor_v1.json
    # 2) run this script:
    #      python scripts/collect_hotpotqa_trajectories.py \
    #        --hotpot-json data/hotpot_dev_distractor_v1.json \
    #        --out storage/predictor_train/hotpotqa_dev_smoke.jsonl \
    #        --n 20 --T-max 8

The confidence controller stays OFF (T7 requires unshortened trajectories).
Every parent LLM call carries `logprobs=True, top_logprobs=cfg.top_b` so
each round records a real posterior on `AgentSession.posteriors`.

This is a SMOKE script — 20 queries × ~4 rounds ≈ 80 LLM calls, ~5-10 min
on hepai/deepseek-v4-pro. Scale up to n=500-1000 once the JSONL shape has
been eyeballed.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Make `scripts.*` importable when this file is run as a top-level script
# (`python scripts/collect_hotpotqa_trajectories.py …`). Without this the
# `from scripts.eval_answer_hotpotqa import _build_agent_factory` below
# fails with ModuleNotFoundError. Same trick probe_8round_qid.py uses.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from hepuke import HepUKE, load_config
from hepuke.agent.confidence.collect_predictor_data import collect_trajectories
from hepuke.agent.confidence.hotpotqa_queries import (
    load_hotpotqa_json,
    sample_examples,
)
from hepuke.agent.loop import RagAgent  # noqa: F401 — kept for backwards-compat imports
from hepuke.config import ConfidenceConfig  # noqa: F401 — see above
from hepuke.core.llm import LLMClient
# Reuse the eval-side factory + preprocessor so collector, eval, and probe
# share ONE agent-construction path. Any drift (tool set, force_first_tool,
# search_return_top, sparse preprocessor) would silently skew the collected
# trajectories vs. the runtime the controller is later evaluated on.
from scripts.eval_answer_hotpotqa import (  # type: ignore
    DEFAULT_SLUG2TITLE,
    DEFAULT_TITLE2SENTS,
    _build_agent_factory,
    _hotpot_preprocess,
    _run_one_query,
)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hotpot-json", required=True, type=Path,
                   help="Path to hotpot_dev_distractor_v1.json (or a sample).")
    p.add_argument("--out", required=True, type=Path,
                   help="Output JSONL of per-round records.")
    p.add_argument("--n", type=int, default=20,
                   help="How many queries to sample (default 20 for a smoke).")
    p.add_argument("--T-max", type=int, default=8,
                   help="Fixed round budget per query (paper default: 8).")
    p.add_argument("--top-b", type=int, default=8,
                   help="Top-B answer tokens sent as logprobs (default 8).")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--collection", default=None,
                   help="Milvus collection to connect. Defaults to cfg.milvus.default_collection.")
    p.add_argument("--model", default=None,
                   help="Override cfg.llm.model for this run (e.g. openai/gpt-4o-mini). "
                        "Leaves config.yaml untouched so other tasks are unaffected.")
    p.add_argument("--no-subagents", action="store_true",
                   help="Force cfg.agent.use_subagents=False for this run. Needed on "
                        "collections like hotpotqa_distractor where the memo_search "
                        "subagent is irrelevant and eats the whole step budget.")
    p.add_argument("--shard", default=None,
                   help='"i/N" — process only qids where index%%N==i in the '
                        'sampled list. Lets multiple processes cover the same '
                        '--n sample without overlap. Combine with --skip-existing '
                        'on resume runs. Example: 0/4 through 3/4.')
    p.add_argument("--skip-existing", type=Path, default=None, nargs="+",
                   help="One or more JSONL files whose query_ids should be "
                        "skipped. Use to resume after a crash — pass the "
                        "previous run's output(s). Multiple paths allowed.")
    p.add_argument("--posterior-model", default="openai/gpt-4.1-mini",
                   help="LLM used only for the shadow posterior probe. "
                        "Independent of the main agent model; must match "
                        "the eval-side --posterior-model so predictor "
                        "training / inference distributions align.")
    p.add_argument("--workers", type=int, default=1,
                   help="Concurrent queries in-process (threadpool). "
                        "workers=1 (default) preserves the original serial "
                        "loop; higher values dispatch queries in parallel. "
                        "Each worker mints its own RagAgent per query.")
    p.add_argument("--auto-resume", action="store_true",
                   help="On start, read `--out` (if present) and skip any "
                        "qids already recorded there. Idempotent — safe to "
                        "re-run the exact same command after a crash.")
    # Eval-shape scoring fields — same LLMs / files as eval_answer_hotpotqa.
    p.add_argument("--sp-model", default="zhipu/glm-5.3",
                   help="LLM used by the sp picker (supporting-fact selector). "
                        "Default glm-5.3 matches the eval-side default and "
                        "has proven far more reliable at strict-JSON output "
                        "than v4-flash (which produced {} under thinking-"
                        "mode-tight max_tokens budgets).")
    p.add_argument("--short-model", default="zhipu/glm-5.3",
                   help="LLM used to compress the agent's long answer into "
                        "a HotpotQA-length short span. Default glm-5.3 "
                        "matches sp-model (same rationale).")
    p.add_argument("--sp-max-titles", type=int, default=10,
                   help="Cap on distinct titles fed to the sp picker.")
    p.add_argument("--slug2title", type=Path, default=Path(DEFAULT_SLUG2TITLE),
                   help="doc_id (slug) → Wikipedia title map, same file as eval.")
    p.add_argument("--title2sents", type=Path, default=Path(DEFAULT_TITLE2SENTS),
                   help="Wikipedia title → [sentences] map, same file as eval.")
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    log = logging.getLogger("collect_hotpotqa")

    cfg = load_config()
    if args.collection:
        cfg.milvus.default_collection = args.collection
    if args.model:
        log.info("overriding cfg.llm.model: %s -> %s", cfg.llm.model, args.model)
        cfg.llm.model = args.model
    # Match the eval script: hotpot corpus doesn't need memo_search, and
    # subagents would burn the step budget on empty fan-outs. Force it OFF
    # unconditionally — the --no-subagents flag stays as a no-op alias.
    if cfg.agent.use_subagents:
        log.info("forcing cfg.agent.use_subagents: %s -> False (hotpot)",
                 cfg.agent.use_subagents)
        cfg.agent.use_subagents = False
    rag = HepUKE(cfg)
    rag.connect_collection(cfg.milvus.default_collection)

    # Install the hotpot-aware sparse preprocessor the eval script uses.
    # Physics preprocess leaves English questions untouched; BGE-M3 sparse
    # then scores "who/the/is" noise. Same rationale + guard as eval.
    if getattr(rag.retriever, "_preprocess", None) is not None:
        rag.retriever._preprocess = _hotpot_preprocess
        log.info("installed hotpot-aware sparse preprocessor")
    else:
        log.warning("retriever._preprocess is None; cfg.retrieval.with_query_preprocess "
                    "may be off. Skipping hotpot preprocessor install.")

    examples = load_hotpotqa_json(args.hotpot_json)
    log.info("loaded %d hotpot examples from %s", len(examples), args.hotpot_json)
    sub = sample_examples(examples, n=args.n, seed=args.seed)
    log.info("sampled %d queries (seed=%d)", len(sub), args.seed)

    # Optional: shard the sample across parallel processes.
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

    # Optional: skip qids already covered by previous run(s).
    if args.skip_existing:
        import json as _json
        done: set[str] = set()
        for p in args.skip_existing:
            if not p.exists():
                log.warning("skip-existing file %s does not exist — ignoring", p); continue
            with p.open("r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line: continue
                    try:
                        qid = _json.loads(line).get("query_id")
                    except Exception:
                        continue
                    if qid:
                        done.add(str(qid))
        before = len(sub)
        sub = [ex for ex in sub if getattr(ex, "id", None) not in done]
        log.info("skip-existing: %d already covered → %d remaining",
                 before - len(sub), len(sub))

    # Auto-resume: read args.out (if exists) and add its qids to the skip set.
    # Deliberately runs AFTER --skip-existing so both stack; idempotent
    # rerun of the same command becomes safe.
    if args.auto_resume and args.out.exists():
        import json as _json
        done: set[str] = set()
        with args.out.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                try:
                    qid = _json.loads(line).get("query_id")
                except Exception:
                    continue
                if qid:
                    done.add(str(qid))
        before = len(sub)
        sub = [ex for ex in sub if getattr(ex, "id", None) not in done]
        log.info("auto-resume from %s: %d already covered → %d remaining",
                 args.out, before - len(sub), len(sub))

    # Dedicated posterior client. `openai/gpt-4.1-mini` is cheap and
    # deterministic under `temperature=0`; the loop caps the shadow probe
    # at 8 tokens so cost stays low even on T_max=8 trajectories.
    posterior_llm = LLMClient(
        base_url=cfg.llm.base_url,
        api_key=cfg.llm.api_key,
        model=args.posterior_model,
        temperature=0.0,
        max_tokens=8,
    )
    log.info("posterior LLM = %s (shadow probes only)", args.posterior_model)
    log.info("agent LLM = %s (main tool-calling loop, from cfg.llm.model)", cfg.llm.model)

    # Sp picker + short-answer LLMs — same models eval uses. Their outputs
    # (pred_sp, short_ans) get merged into every predictor training row so
    # downstream training can filter by gold-hit / answer F1 without re-
    # running the agent.
    import json as _json_mod
    sp_llm = LLMClient(
        base_url=cfg.llm.base_url,
        api_key=cfg.llm.api_key,
        model=args.sp_model,
        temperature=0.0,
        max_tokens=1024,
    )
    short_llm = LLMClient(
        base_url=cfg.llm.base_url,
        api_key=cfg.llm.api_key,
        model=args.short_model,
        temperature=0.0,
        # 2048 is a "safe headroom" ceiling, not a per-call cost. glm-5.3
        # has a thinking-mode that spends tokens on internal reasoning
        # before writing the final JSON — max_tokens=64 was too tight
        # and produced empty {} objects. Real completion_tokens per call
        # remain ~30-80 (just the {"short": "..."} envelope).
        max_tokens=2048,
    )
    log.info("sp LLM = %s, short-ans LLM = %s", args.sp_model, args.short_model)
    slug2title = _json_mod.loads(args.slug2title.read_text())
    title2sents = _json_mod.loads(args.title2sents.read_text())
    log.info("loaded %d slugs, %d titles", len(slug2title), len(title2sents))

    # Use the eval-side factory verbatim (same tool set = {search, rewrite},
    # same search_return_top=1) so the collected trajectories match what the
    # controller will see at eval time.
    # ctrl_on=False keeps the confidence controller from firing and
    # shortening trajectories (T7 unshortened-trajectory requirement); the
    # factory sets collect_logprobs_only=True when ctrl_on is False, which
    # is exactly what this collector needs.
    _base_factory = _build_agent_factory(
        cfg, rag,
        T_max=args.T_max, top_b=args.top_b, ctrl_on=False,
        posterior_llm=posterior_llm,
    )

    # zhipu/glm-5.2 in "Thinking mode" rejects function-name tool_choice
    # (Upstream 400001: "Thinking mode does not support this tool_choice").
    # The eval factory sets force_first_tool="rewrite" → step-0 sends
    # tool_choice={"type":"function","function":{"name":"rewrite"}}, which
    # is exactly the unsupported form. Neutralise it on the collector-built
    # agents only; eval/probe paths untouched.
    def factory():
        agent = _base_factory()
        agent.force_first_tool = None
        # Serialize tool dispatch: same-turn parallel `search` calls run three
        # M3 encodes on the local GPU in three threads simultaneously, which
        # trips a libcuda segfault on driver 580.95.05 + Blackwell (see
        # storage/collector_smoke_w2.log). Serial dispatch removes the race.
        # agent.parallel_tool_calls = False
        return agent

    workers = max(1, int(args.workers))
    log.info("workers = %d", workers)

    # Shared helper: run one example through _run_one_query with a
    # TrajectoryCollector installed, then merge eval-shape scoring fields
    # (pred_ans / short_ans / pred_sp / gold_ans / gold_sp / doc_ids /
    # rounds / tool_counts / search_calls / latency_s) into EVERY
    # predictor training row for that query. This way one collector run
    # produces both:
    #   * predictor training rows (per-round φ + posterior + g_next)
    #   * eval-shape scoring (agent long answer + short answer + sp),
    #     ready to be fed to scripts/hotpot_score.py for EM/F1.
    import json as _json
    import uuid as _uuid
    from hepuke.agent.confidence.collect_predictor_data import TrajectoryCollector

    _EVAL_ROW_KEYS = (
        "pred_ans", "short_ans", "short_error",
        "pred_sp", "sp_error", "missing_slugs_count",
        "gold_ans", "gold_sp",
        "doc_ids_eval",   # renamed to avoid clashing with per-round `doc_ids`
        "rounds", "terminated", "tool_counts",
        "search_calls_eval",  # renamed: eval counts search_calls too
        "non_hotpot_search",
        "latency_s", "error",
    )

    def _score_from_result(ex, result: dict, rows: list[dict]) -> list[dict]:
        """Score using RagAgent.run()'s return dict. Idempotent: mutates
        rows in-place, returns them. `result` must have `answer` and
        `observed_hits` and `trace` fields, per RagAgent's contract."""
        if not result:
            return rows
        pred_ans = (result.get("answer") or "").strip()
        observed = list(result.get("observed_hits") or [])
        trace = list(result.get("trace") or [])

        # doc_ids in retrieval order, deduped
        seen: set[str] = set()
        doc_ids_eval: list[str] = []
        for h in observed:
            did = h.get("doc_id") if isinstance(h, dict) else None
            if did and did not in seen:
                seen.add(did)
                doc_ids_eval.append(did)

        from scripts.eval_answer_hotpotqa import (  # local import to avoid cycles
            _assemble_candidates,
            _extract_short_ans,
            _pick_sp,
        )
        short_ans, short_err = _extract_short_ans(
            short_llm, question=ex.question, long_answer=pred_ans,
        )
        labeled, missing = _assemble_candidates(
            doc_ids_eval, slug2title, title2sents,
            max_titles=args.sp_max_titles,
        )
        pred_sp, sp_err = _pick_sp(
            sp_llm, question=ex.question,
            pred_answer=short_ans or pred_ans, labeled=labeled,
        )

        tool_counts: dict[str, int] = {}
        for e in trace:
            t = e.get("tool")
            if t:
                tool_counts[t] = tool_counts.get(t, 0) + 1
        search_calls_eval = sum(1 for e in trace if e.get("tool") == "search")

        extras = {
            "pred_ans": pred_ans,
            "short_ans": short_ans,
            "pred_sp": pred_sp,
            "gold_ans": ex.answer,
            "gold_sp": [list(sp) for sp in ex.supporting_facts],
            "doc_ids_eval": doc_ids_eval,
            "rounds": int(result.get("steps") or 0),
            "terminated": result.get("terminated"),
            "tool_counts": tool_counts,
            "search_calls_eval": search_calls_eval,
        }
        if short_err:
            extras["short_error"] = short_err
        if sp_err:
            extras["sp_error"] = sp_err
        if missing:
            extras["missing_slugs_count"] = len(missing)
        for r in rows:
            r.update(extras)
        return rows

    if workers == 1:
        # Serial path.
        args.out.parent.mkdir(parents=True, exist_ok=True)
        total = 0
        n_done = 0
        n_total = len(sub)
        with args.out.open("a", encoding="utf-8") as f_out:
            for ex in sub:
                qid = str(getattr(ex, "id", "") or "")
                qtext = str(getattr(ex, "question", "") or "")
                if not qtext:
                    log.info("  [%d/%d] %s → 0 rows (empty question)",
                             n_done + 1, n_total, qid)
                    n_done += 1
                    continue
                agent = factory()
                col = TrajectoryCollector(query_id=qid, query=qtext, T_max=args.T_max)
                agent.round_observer = col
                sid = _uuid.uuid4().hex
                err: str | None = None
                result: dict | None = None
                try:
                    result = agent.run(qtext, max_steps=args.T_max, session_id=sid)
                except Exception as e:  # noqa: BLE001
                    err = f"{type(e).__name__}: {e}"
                finally:
                    col.finalize()
                rows = list(col.records)
                if err is not None:
                    for r in rows:
                        r["agent_error"] = err
                elif result is not None:
                    try:
                        rows = _score_from_result(ex, result, rows)
                    except Exception as e:  # noqa: BLE001
                        for r in rows:
                            r["score_error"] = f"{type(e).__name__}: {e}"
                n_rows = 0
                for row in rows:
                    f_out.write(_json.dumps(row, ensure_ascii=False))
                    f_out.write("\n")
                    n_rows += 1
                f_out.flush()
                total += n_rows
                n_done += 1
                log.info("  [%d/%d] %s → %d rows%s",
                         n_done, n_total, qid, n_rows,
                         f" ERROR {err}" if err else "")
        log.info("wrote %d rows → %s", total, args.out)
        return 0

    # Parallel path.
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_lock = threading.Lock()
    total = 0
    n_done = 0
    n_total = len(sub)
    f_out = args.out.open("a", encoding="utf-8")

    def _run_one(ex) -> tuple[str, int, str | None]:
        nonlocal total, n_done
        qid = str(getattr(ex, "id", "") or "")
        qtext = str(getattr(ex, "question", "") or "")
        if not qtext:
            return qid, 0, "empty-question"
        agent = factory()
        if getattr(agent, "_confidence_controller", None) is not None:
            return qid, 0, "controller-on"
        col = TrajectoryCollector(query_id=qid, query=qtext, T_max=args.T_max)
        agent.round_observer = col
        sid = _uuid.uuid4().hex
        err: str | None = None
        result: dict | None = None
        try:
            result = agent.run(qtext, max_steps=args.T_max, session_id=sid)
        except Exception as e:  # noqa: BLE001
            err = f"{type(e).__name__}: {e}"
        finally:
            col.finalize()
        rows = list(col.records)
        if err is not None:
            for r in rows:
                r["agent_error"] = err
        elif result is not None:
            try:
                rows = _score_from_result(ex, result, rows)
            except Exception as e:  # noqa: BLE001 — scoring must not kill the row
                for r in rows:
                    r["score_error"] = f"{type(e).__name__}: {e}"
        n_rows = 0
        with write_lock:
            for row in rows:
                f_out.write(_json.dumps(row, ensure_ascii=False))
                f_out.write("\n")
                n_rows += 1
            f_out.flush()
            total += n_rows
            n_done += 1
            log.info("  [%d/%d] %s → %d rows%s",
                     n_done, n_total, qid, n_rows,
                     f" ERROR {err}" if err else "")
        return qid, n_rows, err

    try:
        with ThreadPoolExecutor(
            max_workers=workers, thread_name_prefix="collect",
        ) as pool:
            futures = [pool.submit(_run_one, ex) for ex in sub]
            for _ in as_completed(futures):
                pass
    finally:
        f_out.close()

    log.info("wrote %d rows → %s", total, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
