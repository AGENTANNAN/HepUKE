"""BESIII QA under the HotpotQA protocol — parity runner v1.

Why
---
`storage/submission_60h/protocol_parity_hotpot_besiii_v1.md` audits the frozen
QA120 protocol against `scripts/eval_answer_hotpotqa.py` and finds nine of
fourteen dimensions differ. This runner removes those differences so the BESIII
test measures the same RAG loop HotpotQA measures.

Aligned to HotpotQA (each line is a change from QA120):

* corpus  — `besiii_dsl_block_v1`: one document per DSL logical block, one
  child + one identical parent row, like HotpotQA's one-paragraph documents
  (built by `scripts/index_besiii_dsl_blocks_v1.py`).
* tools   — `{search, rewrite}` with `force_first_tool="rewrite"`.
* depth   — `search_return_top=1`: one complete evidence unit per lane.
* clock   — `T_max=8` with native stopping. No `ForcedSearchLLM`; the agent
  decides each round whether to search, exactly as in HotpotQA.
* arms    — `ctrl_off` / `ctrl_on`, not fixed2/fixed4/adaptive.
* prompt  — generic `build_system_prompt`, not a BESIII-specific system prompt.
  HotpotQA keeps its domain knowledge in the extractor, so we do too.
* answer  — free-form long answer, then a short-answer extractor LLM, the
  direct analogue of HotpotQA's `_SHORT_SYSTEM` compression step.
* support — a separate supporting-block picker LLM over the retrieved blocks,
  the analogue of HotpotQA's `_SP_SYSTEM` call, scored as Sup-EM / Sup-F1 /
  Joint with `scripts/hotpot_score.py`'s own set rules.
* metrics — HotpotQA token EM/F1 **and** the BESIII typed strict scorer, so the
  result is readable next to both benchmarks.

Held constant against the BESIII references (`besiii_paired_dev_retrieverfix`,
native .3000 / hotpot_controller .2639): main model
`deepseek-ai/deepseek-v4.1-flash` at temperature 0 / 2048 tokens, posterior
probe `openai/gpt-4.1-mini` with `top_b=8`, and the same predictor/fusion
checkpoints. Flash reasoning controls are frozen into every call, because the
v4 smoke burned its completion budget on hidden reasoning and returned empty
content.

Boundaries
----------
* The block corpus excludes paper markdown (all 1,619/1,619 reviewed-QA gold
  refs have `source=dsl`). That makes this corpus easier than QA120's, so a
  score here is **not** comparable to the frozen `.2034565/.2524603/.2518254`
  headline. Attributing any change to protocol rather than corpus needs a third
  arm running this protocol on `besiii_dsl_v1`.
* The 12 smoke questions are previously exposed development questions. This is
  a development smoke, not a held-out benchmark.
* HotpotQA's own main-generator model is not recorded in its run logs, so this
  runner pins the BESIII model instead of claiming model parity with HotpotQA.
* Writes only into its own versioned directory. Never touches a frozen QA120 /
  DSL20 manifest, trace, score, or an existing Milvus collection.

Usage
-----
    python scripts/run_besiii_parity_hotpot_v1.py --prepare
    python scripts/run_besiii_parity_hotpot_v1.py --validate
    python scripts/run_besiii_parity_hotpot_v1.py --run --allow-endpoint
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import threading
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

COLLECTION = "besiii_dsl_block_v1"
QA_PATH = ROOT / "data" / "besiii_physicsqa" / "qa.jsonl"
# Produced by scripts/index_besiii_dsl_blocks_v1.py from the shipped
# data/besiii_dsl_index/chunks.jsonl + data/besiii_physicsqa/qa.jsonl.
GOLD_BLOCKS = ROOT / "data" / "besiii_block_index" / "gold_blocks.json"
BLOCK_MANIFEST = ROOT / "data" / "besiii_block_index" / "manifest.json"
OUT_DIR = ROOT / "storage" / "besiii_parity_hotpot_v1"
# Optional: only read when --qid-source qa120 is used. Not shipped.
QA120_SELECTION = ROOT / "storage" / "qa120_selection.json"

# The 12 development questions used by besiii_native_vs_hotpot_controller_dev_v1
# so this run has native .3000 / hotpot_controller .2639 as reference points.
SMOKE_QIDS = ["besiii-A-4199", "besiii-A-0042", "besiii-C-0020",
              "besiii-A-0745", "besiii-A-1688", "besiii-B-0068",
              "besiii-A-0098", "besiii-A-2286", "besiii-C-0022",
              "besiii-A-1866", "besiii-A-2102", "besiii-C-0013"]


def qa120_qids() -> list[str]:
    """The frozen QA120 selection, so the expansion pairs against the same
    question set the frozen headline used."""
    return [r["qid"] for r in
            json.loads(QA120_SELECTION.read_text(encoding="utf-8"))["rows"]]


def full_qids() -> list[str]:
    """All 980 reviewed QA questions — the full BESIII multihop v2 set."""
    result: list[str] = []
    for line in QA_PATH.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            result.append(json.loads(line)["qid"])
    return result


# Pre-registered before the expansion ran. The target claim is a within-corpus,
# within-protocol non-inferiority statement: ctrl_on holds the answer score
# within delta while spending less. Fixed here so it cannot be chosen after
# seeing the result.
PRE_REGISTRATION = {
    "primary_claim": "ctrl_on is non-inferior to ctrl_off on typed_score "
                     "within delta, at strictly lower cost",
    "primary_metric": "typed_score",
    "equivalence_margin_delta": 0.10,
    "test": "one-sided paired bootstrap on (ctrl_on - ctrl_off); claim holds "
            "iff the 95% lower bound exceeds -delta",
    "bootstrap_B": 10000,
    "bootstrap_seed": 20260921,
    "cost_metrics": ["total_tokens", "rounds", "elapsed_s", "search_calls"],
    "cost_criterion": "two-sided 95% paired bootstrap interval excludes zero",
    "secondary_metrics": ["ans_em", "ans_f1", "sp_em", "sp_f1", "joint_f1"],
    "secondary_note": "Ans-EM cannot support delta<=0.05 at the smoke's point "
                      "estimate; reported descriptively, not as a claim.",
    "sizing_basis": "smoke typed paired sd 0.3321 implies N~58 for delta=0.10",
    "no_arm_added_after_seeing_results": True,
    "declared_at": "before the 120-question expansion started",
}

ARMS = ["ctrl_off", "ctrl_on"]
T_MAX = 8
TOP_B = 8
SP_MAX_BLOCKS = 10

# Block kinds in besiii_dsl_block_v1, i.e. the values of the Milvus
# `section_type` column that the `search` tool's `section_types` whitelist
# filters on. Ordered by share of the 1,619 reviewed-QA gold references.
BLOCK_KINDS = ["exclusive_mc", "selection", "decay_card", "tag_analysis",
               "dataset", "execute", "header", "algorithm", "note"]

# Two configurations. `parity` is the faithful HotpotQA reading and is what
# besiii_parity_hotpot_v1 / _qa120_v1 ran. `dsl_tuned` applies the three
# retrieval fixes chosen after the QA120 diagnosis (gold-block recall 0.222,
# 60% of failures from answering off the wrong block). Each fix is a deliberate
# DEVIATION from HotpotQA and is recorded as such, so `dsl_tuned` must be
# reported as its own arm and never substituted for the parity numbers.
VARIANTS = {
    "parity": {
        "search_return_top": 1,
        "rewrite_prompt_style": None,      # inherit config.yaml ("generic")
        "kind_filter_hint": False,
        "deviations_from_hotpotqa": [],
    },
    "dsl_tuned": {
        "search_return_top": 3,
        "rewrite_prompt_style": "dsl",
        "kind_filter_hint": True,
        "deviations_from_hotpotqa": [
            "search_return_top 1 -> 3: BESIII averages 1.29 gold blocks per "
            "question and `header` blocks crowd rank 1, so top-1 is brittle. "
            "HotpotQA uses 1.",
            "rewrite prompt style generic -> dsl: the generic prompt targets "
            "an 'encyclopedia-style corpus' with a Casablanca few-shot and was "
            "measured emitting paper-register queries against Ruby source. "
            "HotpotQA's corpus really is encyclopedic, so this fix has no "
            "HotpotQA counterpart.",
            "kind filter: the agent is told the nine block kinds and to pass "
            "`section_types`. This is a structural prior over the corpus that "
            "HotpotQA's flat paragraph index does not have.",
        ],
    },
}

KIND_FILTER_HINT = (
    "\n\nThe corpus is BOSS Ruby DSL source code. Every document is one "
    "logical block of one analysis program, and its `section_type` is exactly "
    "one of: " + ", ".join(BLOCK_KINDS) + ". The `search` tool accepts a "
    "`section_types` whitelist — pass the kind that would hold the setting "
    "you need (`selection` for track/photon/PID/kinematic-fit cuts, `dataset` "
    "for sample identifiers, `exclusive_mc` for generated event counts and "
    "sample names, `decay_card` for decay models such as PHSP or ETA_DALITZ, "
    "`tag_analysis` for D-tag/Lambda_c-tag quantities, `execute` for job and "
    "output wiring). `header` holds only leading comments and almost never "
    "holds a setting, so do not filter to it. If a filtered search returns "
    "nothing useful, retry without `section_types`. Dataset and decay-model "
    "names are code identifiers: report their literal underscore spelling "
    "from the source, not a descriptive label."
)
def _model_id(frozen: str) -> str:
    """Frozen paper identifier by default; override with $DSL_MODEL or llm.model."""
    import os as _os
    env = _os.environ.get("DSL_MODEL")
    if env:
        return env
    try:
        import yaml as _yaml
        m = str(((_yaml.safe_load((ROOT / "config.yaml").read_text()) or {})
                 .get("llm") or {}).get("model") or "")
        if m and not m.startswith("${"):
            return m
    except Exception:
        pass
    return frozen


MAIN_MODEL = _model_id("deepseek-flash")
POSTERIOR_MODEL = _model_id("deepseek-flash")
AUX_MODEL = "deepseek-flash"          # short-ans + sp extraction

SOURCE_FILES = [
    "scripts/run_besiii_parity_hotpot_v1.py",
    "scripts/index_besiii_dsl_blocks_v1.py",
    "scripts/besiii_submission_scoring.py",
    "scripts/hotpot_score.py",
    "hepuke/agent/loop.py",
    "hepuke/agent/defaults.py",
    "hepuke/core/hybrid_retriever.py",
    "hepuke/core/llm.py",
]

# --------------------------------------------------------------------------- #
# Prompts — BESIII content, HotpotQA shape                                    #
# --------------------------------------------------------------------------- #
_SHORT_SYSTEM = (
    "You extract the shortest possible answer span from a long answer to a "
    "question about a particle-physics analysis program. Rules: "
    "(1) output ONLY tokens that appear in the long answer — no paraphrase; "
    "(2) preserve comparison operators, signs, units, dataset identifiers and "
    "particle identities exactly as written; "
    "(3) a dataset or decay-model name is a code identifier — keep its literal "
    "spelling, never substitute a descriptive label or a paper ID; "
    "(4) for a question asking for a set, list every item separated by "
    "semicolons; "
    "(5) if the long answer does not determine an answer, output UNKNOWN; "
    '(6) return ONLY a JSON object of shape {"short": "..."} — no prose, no '
    "fences."
)

_SP_SYSTEM = (
    "You select the supporting evidence blocks that JUSTIFY a given short "
    "answer to a question about a particle-physics analysis program. "
    'Return ONLY a JSON object of the exact shape {"supporting_blocks": '
    '["<block_id>", ...]}. Rules: '
    "(1) every block_id must appear verbatim in the candidate list; "
    "(2) most questions are justified by ONE block, some by two blocks of "
    "different kinds (for example a dataset block plus a selection block). "
    "When two blocks jointly carry the justification, include BOTH; "
    "(3) prefer recall over precision — a missing required block zeroes "
    "recall, while one extra block only costs some precision; "
    "(4) never repeat the same block_id; "
    "(5) no prose, no code fences — the body must start with `{` and end "
    "with `}`."
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prepare", action="store_true",
                   help="freeze the manifest; no endpoint contact")
    p.add_argument("--validate", action="store_true",
                   help="re-check hashes, collection and gold blocks offline")
    p.add_argument("--run", action="store_true")
    p.add_argument("--allow-endpoint", action="store_true",
                   help="required with --run; guards API spend")
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    p.add_argument("--variant", choices=tuple(VARIANTS), default="parity",
                   help="parity = faithful HotpotQA reading; dsl_tuned = the "
                        "three retrieval fixes, recorded as deviations")
    p.add_argument("--qid-source", choices=("smoke", "qa120", "full", "file"),
                   default="smoke",
                   help="smoke = the 12 development questions; qa120 = the "
                        "frozen 120-question selection; full = all 980; "
                        "file = --qids-file")
    p.add_argument("--qids-file", type=Path, default=None,
                   help="JSON list of qids, required with --qid-source file. "
                        "Used to validate a fix on exactly the subset another "
                        "variant already completed, so the comparison is paired.")
    p.add_argument("--limit", type=int, default=None,
                   help="only the first N qids (debug)")
    p.add_argument("--shard-id", type=int, default=None,
                   help="0-based shard index (requires --num-shards)")
    p.add_argument("--num-shards", type=int, default=None,
                   help="total number of shards (requires --shard-id)")
    p.add_argument("--merge", action="store_true",
                   help="merge shard results under --out-dir into one summary")
    p.add_argument("--z-star", type=float, default=None,
                   help="override confidence.z_star (fusion stop threshold); "
                        "None = use config.yaml")
    p.add_argument("--arms", type=str, default=None,
                   help="comma-separated arms to run, e.g. 'ctrl_on' (default: "
                        "both)")
    return p.parse_args()


def _save(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False),
                    encoding="utf-8")


def _public_config(cfg) -> dict:
    return {
        "collection": COLLECTION,
        "database": cfg.milvus.database,
        "embedding_model": cfg.embedding.model,
        "embedding_provider": cfg.embedding.provider,
        "reranker_model": cfg.reranker.model,
        "retrieval": {
            "similarity_top_k": cfg.retrieval.similarity_top_k,
            "with_reranker": cfg.retrieval.with_reranker,
            "hybrid_enabled": getattr(cfg.retrieval, "hybrid_enabled", None),
            "rrf_k": getattr(cfg.retrieval, "rrf_k", None),
            "with_small_to_big": getattr(cfg.retrieval, "with_small_to_big", None),
            "query_rewrite_llm_enabled": getattr(
                cfg.retrieval, "query_rewrite_llm_enabled", None),
            "query_rewrite_llm_model": getattr(
                cfg.retrieval, "query_rewrite_llm_model", None),
        },
        "model": MAIN_MODEL,
        "temperature": 0.0,
        "max_tokens": 2048,
    }


def build_manifest(qids: list[str], qid_source: str,
                   variant: str = "parity") -> dict:
    from hepuke import load_config
    cfg = load_config()
    cfg.llm.model = MAIN_MODEL
    v = VARIANTS[variant]
    if v["rewrite_prompt_style"]:
        cfg.retrieval.query_rewrite_prompt_style = v["rewrite_prompt_style"]
    is_smoke = qid_source == "smoke"
    return {
        "protocol": (("besiii_parity_hotpot_v1" if is_smoke
                      else "besiii_parity_hotpot_qa120_v1") if variant == "parity"
                     else ("besiii_dsl_tuned_smoke_v1" if is_smoke
                           else "besiii_dsl_tuned_qa120_v1")),
        "role": ("development_smoke_protocol_parity_with_hotpotqa" if is_smoke
                 else "qa120_selection_expansion_noninferiority_at_lower_cost"),
        "qid_source": qid_source,
        "n_qids": len(qids),
        "pre_registration": PRE_REGISTRATION,
        "smoke_reference": (None if is_smoke else
                            "storage/submission_60h/besiii_parity_hotpot_v1/"
                            "report.md"),
        "third_arm_deferred": (
            "A 'parity protocol on besiii_dsl_v1' arm would separate protocol "
            "from corpus. Deliberately not run: the primary claim is a "
            "within-corpus ctrl_on vs ctrl_off comparison and does not need it. "
            "It can be added later as a separate arm without redoing this run."),
        "reference_audit": "storage/submission_60h/"
                           "protocol_parity_hotpot_besiii_v1.md",
        "qids": qids,
        "arms": ARMS,
        "t_max": T_MAX,
        "top_b": TOP_B,
        "variant": variant,
        "search_return_top": v["search_return_top"],
        "rewrite_prompt_style": (v["rewrite_prompt_style"]
                                 or cfg.retrieval.query_rewrite_prompt_style),
        "kind_filter_hint": v["kind_filter_hint"],
        "block_kinds": BLOCK_KINDS,
        "deviations_from_hotpotqa": v["deviations_from_hotpotqa"],
        "sp_max_blocks": SP_MAX_BLOCKS,
        "native_stopping": True,
        "forced_search": False,
        "force_first_tool": "rewrite",
        "tools_kept": ["search", "rewrite"],
        "posterior_model": POSTERIOR_MODEL,
        "aux_model": AUX_MODEL,
        "aux_roles": ["short_answer_extractor", "supporting_block_picker"],
        "flash_controls": {"reasoning_effort": "none",
                           "thinking": {"type": "disabled"}},
        "public_config": _public_config(cfg),
        "dataset_path": str(QA_PATH.relative_to(ROOT)),
        "dataset_sha256": _sha256(QA_PATH),
        "qa120_selection_sha256": (None if is_smoke
                                   else _sha256(QA120_SELECTION)),
        "runner_lineage": (
            "This runner was extended after the 12-question smoke to add "
            "--qid-source and the pre-registration block. The smoke manifest's "
            "recorded hash for scripts/run_besiii_parity_hotpot_v1.py therefore "
            "refers to the pre-extension file; the smoke's raw results and "
            "report are unchanged and must not be re-validated against this "
            "version. No scoring, prompt, model, tool, or retrieval parameter "
            "was altered."),
        "gold_blocks_sha256": _sha256(GOLD_BLOCKS),
        "block_index_manifest_sha256": _sha256(BLOCK_MANIFEST),
        "predictor_ckpt": cfg.confidence.predictor_ckpt,
        "fusion_ckpt": cfg.confidence.fusion_ckpt,
        "predictor_sha256": _sha256(Path(cfg.confidence.predictor_ckpt)),
        "fusion_sha256": _sha256(Path(cfg.confidence.fusion_ckpt)),
        "source_hashes": {f: _sha256(ROOT / f) for f in SOURCE_FILES},
        "reference_points": {
            "besiii_paired_dev_retrieverfix": {
                "native": {"score": 0.3000, "rounds": 4.0833},
                "hotpot_controller": {"score": 0.2639, "rounds": 2.0},
                "note": "same 12 qids, besiii_dsl_v1 corpus, "
                        "search_return_top=3, max_steps=6",
            },
            "besiii_fixed120_v1": {
                "fixed2": 0.2034565, "fixed4": 0.2524603,
                "adaptive": 0.2518254,
                "note": "different corpus and different 120 questions; "
                        "not a paired comparison",
            },
        },
        "boundaries": [
            "Block corpus excludes paper markdown; easier than QA120's corpus, "
            "so scores are not comparable to the frozen headline.",
            "12 previously exposed development questions; not held out.",
            "HotpotQA's own generator model is unrecorded; model parity is not "
            "claimed.",
            "Native stopping means rounds are an outcome, not a budget.",
        ],
    }


# --------------------------------------------------------------------------- #
# Endpoint helpers                                                            #
# --------------------------------------------------------------------------- #
def attach_usage(llm, role: str, events: list, lock: threading.Lock) -> None:
    """Record every completion's usage and freeze the flash controls."""
    original = llm.client.chat.completions.create

    def wrapped(**kw):
        t0 = time.monotonic()
        if "flash" in (kw.get("model") or ""):
            kw["reasoning_effort"] = "none"
            body = dict(kw.get("extra_body") or {})
            body["thinking"] = {"type": "disabled"}
            kw["extra_body"] = body
        if kw.get("stream"):
            kw["stream_options"] = {"include_usage": True}
        try:
            response = original(**kw)
        except Exception as exc:
            with lock:
                events.append({"role": role, "status": "error",
                               "error_type": type(exc).__name__,
                               "elapsed_s": time.monotonic() - t0})
            raise
        with lock:
            usage = getattr(response, "usage", None)
            events.append({"role": role, "status": "ok",
                           "usage": usage.model_dump() if usage else None,
                           "elapsed_s": time.monotonic() - t0})
        return response

    llm.client.chat.completions.create = wrapped


def extract_short(llm, question: str, long_answer: str) -> tuple[str, str | None]:
    """HotpotQA `_extract_short_ans` analogue. Degrades to the long answer."""
    if not (long_answer or "").strip():
        return "", "empty_long_answer"
    prompt = (f"Question: {question}\nLong answer: {long_answer}\n\n"
              'Return {"short": "<shortest exact span>"}.')
    try:
        obj = llm.ask_json(prompt, system_prompt=_SHORT_SYSTEM)
    except Exception as exc:
        return long_answer, f"llm_error:{type(exc).__name__}"
    val = obj.get("short") if isinstance(obj, dict) else None
    if not isinstance(val, str) or not val.strip():
        return long_answer, "bad_shape"
    return val.strip(), None


def pick_supporting(llm, question: str, short_ans: str,
                    candidates: list[dict]) -> tuple[list[str], str | None]:
    """HotpotQA `_predict_sp` analogue over retrieved blocks."""
    if not candidates:
        return [], "no_candidates"
    listing = "\n\n".join(
        f'block_id: {c["block_id"]}\n{c["content"]}' for c in candidates)
    prompt = (f"Question: {question}\nShort answer: {short_ans}\n\n"
              f"Candidate blocks:\n{listing}\n\n"
              'Return {"supporting_blocks": ["<block_id>", ...]}.')
    allowed = {c["block_id"] for c in candidates}
    try:
        obj = llm.ask_json(prompt, system_prompt=_SP_SYSTEM)
    except Exception as exc:
        return [], f"llm_error:{type(exc).__name__}"
    vals = obj.get("supporting_blocks") if isinstance(obj, dict) else None
    if not isinstance(vals, list):
        return [], "bad_shape"
    out: list[str] = []
    for v in vals:
        if isinstance(v, str) and v in allowed and v not in out:
            out.append(v)
    return out, None if out else "no_valid_block_id"


def fetch_blocks(store, doc_ids: list[str]) -> dict[str, str]:
    """Full block text by doc_id — HotpotQA rebuilds candidates from the corpus
    rather than from truncated retrieval snippets, so we do the same."""
    if not doc_ids:
        return {}
    expr = "doc_id in [" + ",".join(json.dumps(d) for d in doc_ids) + "]"
    rows = store.client.query(collection_name=COLLECTION, filter=expr,
                              output_fields=["doc_id", "chunk_type", "content"],
                              limit=4 * len(doc_ids))
    out: dict[str, str] = {}
    for r in rows:
        if r.get("chunk_type") == "parent" or r["doc_id"] not in out:
            out[r["doc_id"]] = r.get("content") or ""
    return out


# --------------------------------------------------------------------------- #
# Scoring                                                                     #
# --------------------------------------------------------------------------- #
def _contain_norm(s: str) -> str:
    import re
    return re.sub(r'[\s"\'`,]+', '', str(s or "").lower())


def contain_acc(gold: str, response: str) -> float:
    """Contain-Acc: does the gold answer appear verbatim in the response?

    Lenient by design — it ignores surrounding variable names, prose and
    ordering, so it separates "the model did not know / did not retrieve it"
    from "the strict scorer or the extractor threw it away". For set answers
    (semicolon-separated gold) every item must appear.

    Two things it deliberately does NOT absorb: unit-convention differences
    (gold `0.025` vs a response saying `25 MeV` counts as a miss) and the
    redundant-corpus alias problem (gold `psi3686_data` vs `psip_data` counts
    as a miss), because both need a declared semantic rule, not a substring
    test.
    """
    import re
    g, r = _contain_norm(gold), _contain_norm(response)
    if not g:
        return 0.0
    items = [x for x in re.split(r"[;]", str(gold)) if _contain_norm(x)]
    if len(items) > 1:
        return float(all(_contain_norm(x) in r for x in items))
    return float(g in r)


def score_row(short_ans: str, qa_row: dict, pred_sp: list[str],
              gold_sp: list[str], long_answer: str = "") -> dict:
    from besiii_submission_scoring import score_answer
    from hotpot_score import exact_match_score, f1_score

    typed = score_answer(short_ans, qa_row)
    em = float(exact_match_score(short_ans, qa_row.get("answer") or ""))
    f1, prec, rec = f1_score(short_ans, qa_row.get("answer") or "")

    cur, au = set(pred_sp or []), set(gold_sp or [])
    tp = len(cur & au)
    fp = len(cur - au)
    fn = len(au - cur)
    sp_prec = tp / (tp + fp) if (tp + fp) else 0.0
    sp_rec = tp / (tp + fn) if (tp + fn) else 0.0
    sp_f1 = (2 * sp_prec * sp_rec / (sp_prec + sp_rec)
             if (sp_prec + sp_rec) else 0.0)
    sp_em = 1.0 if (fp + fn == 0 and au) else 0.0

    j_prec, j_rec = prec * sp_prec, rec * sp_rec
    j_f1 = (2 * j_prec * j_rec / (j_prec + j_rec)) if (j_prec + j_rec) else 0.0

    gold_ans = qa_row.get("answer") or ""
    ca_long = contain_acc(gold_ans, long_answer)
    ca_short = contain_acc(gold_ans, short_ans)
    return {
        "contain_acc": ca_long,
        "contain_acc_short": ca_short,
        # 1 when the long answer held the gold string but the short-answer
        # extractor dropped it — isolates extractor damage from model error.
        "extractor_lost_answer": float(ca_long > 0 and ca_short == 0),
        "typed_score": typed["score"], "typed_em": typed["em"],
        "scorer_used": typed["scorer_used"],
        "ans_em": em, "ans_f1": f1, "ans_prec": prec, "ans_recall": rec,
        "sp_em": sp_em, "sp_f1": sp_f1, "sp_prec": sp_prec,
        "sp_recall": sp_rec,
        "joint_em": em * sp_em, "joint_f1": j_f1,
    }


# --------------------------------------------------------------------------- #
# Validate                                                                    #
# --------------------------------------------------------------------------- #
def validate(out_dir: Path) -> int:
    manifest = json.loads((out_dir / "manifest.json").read_text())
    problems: list[str] = []
    variant = manifest.get("variant", "parity")
    if variant not in VARIANTS:
        problems.append(f"unknown variant: {variant}")

    for name, want in manifest["source_hashes"].items():
        if _sha256(ROOT / name) != want:
            problems.append(f"source changed: {name}")
    for key, path in (("dataset_sha256", QA_PATH),
                      ("gold_blocks_sha256", GOLD_BLOCKS),
                      ("block_index_manifest_sha256", BLOCK_MANIFEST)):
        if _sha256(path) != manifest[key]:
            problems.append(f"changed after freeze: {path.name}")

    gold = json.loads(GOLD_BLOCKS.read_text())
    for qid in manifest["qids"]:
        if not gold.get(qid):
            problems.append(f"no gold block for {qid}")

    from hepuke import HepUKE, load_config
    cfg = load_config()
    cfg.llm.model = MAIN_MODEL
    cfg.milvus.default_collection = COLLECTION
    # Must match run(): the rewriter's prompt style is chosen at HepUKE
    # construction, so setting it afterwards has no effect.
    cfg.retrieval.query_rewrite_prompt_style = manifest["rewrite_prompt_style"]
    rag = HepUKE(cfg)
    if COLLECTION not in set(rag.get_collections()):
        problems.append(f"collection missing: {COLLECTION}")
    else:
        rag.connect_collection(COLLECTION)
        stats = rag.store.stats(COLLECTION)
        want_rows = 2 * json.loads(
            BLOCK_MANIFEST.read_text())["stats"]["n_blocks"]
        got = stats.get("row_count") if isinstance(stats, dict) else None
        if got != want_rows:
            problems.append(f"row_count {got} != expected {want_rows}")
        from hepuke.agent.defaults import build_default_tools
        tools = build_default_tools(
            rag.retriever, memo_root="", default_collection=COLLECTION,
            schema_profile=rag.schema_profile, use_subagents=False,
            llm=rag.llm, default_with_rerank=cfg.retrieval.with_reranker,
            rewriter=rag.rewriter,
            search_return_top=manifest["search_return_top"])
        names = {t.name for t in tools}
        if "search" not in names:
            problems.append("search tool missing")
        if "rewrite" not in names:
            problems.append("rewrite tool missing — parity requires it")
        want_style = manifest["rewrite_prompt_style"]
        got = getattr(getattr(rag, "rewriter", None), "_system_prompt", "") or ""
        if want_style == "dsl" and "BOSS Ruby DSL source" not in got:
            problems.append("rewriter did not load the dsl prompt style")
        if want_style == "generic" and "encyclopedia-style" not in got:
            problems.append("rewriter did not load the generic prompt style")

    report = {"ok": not problems, "problems": problems,
              "checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    _save(out_dir / "validation.json", report)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if not problems else 3


# --------------------------------------------------------------------------- #
# Run                                                                         #
# --------------------------------------------------------------------------- #
def run(out_dir: Path, limit: int | None,
        shard_id: int | None = None, num_shards: int | None = None,
        z_star: float | None = None, arms: list[str] | None = None) -> int:
    from hepuke import HepUKE, load_config
    from hepuke.agent.defaults import build_default_tools
    from hepuke.agent.loop import RagAgent
    from hepuke.agent.prompts import build_system_prompt
    from hepuke.agent.tools import ToolRegistry
    from hepuke.config import ConfidenceConfig
    from hepuke.core.llm import LLMClient
    from openai import OpenAI

    manifest = json.loads((out_dir / "manifest.json").read_text())
    qids = manifest["qids"]
    # Shard slicing — deterministic, same qid ordering as the manifest
    if shard_id is not None and num_shards is not None:
        qids = [q for i, q in enumerate(qids) if i % num_shards == shard_id]
        out_dir = out_dir / f"shard_{shard_id:02d}_of_{num_shards:02d}"
    if limit:
        qids = qids[:limit]
    gold_blocks = json.loads(GOLD_BLOCKS.read_text())
    qa = {}
    for line in QA_PATH.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            qa[row["qid"]] = row

    cfg = load_config()
    cfg.llm.model = MAIN_MODEL
    cfg.milvus.default_collection = COLLECTION
    cfg.retrieval.query_rewrite_prompt_style = manifest["rewrite_prompt_style"]
    rag = HepUKE(cfg)
    rag.connect_collection(COLLECTION)

    probe = LLMClient(base_url=cfg.llm.base_url, api_key=cfg.llm.api_key,
                      model=POSTERIOR_MODEL, temperature=0.0, max_tokens=8)
    aux = LLMClient(base_url=cfg.llm.base_url, api_key=cfg.llm.api_key,
                    model=AUX_MODEL, temperature=0.0, max_tokens=512)

    events: list = []
    lock = threading.Lock()
    clients = [("main", rag.llm), ("probe", probe), ("aux", aux)]
    if rag.rewriter:
        clients.append(("rewrite", rag.rewriter.llm))
    for role, llm in clients:
        llm.client = OpenAI(base_url=cfg.llm.base_url, api_key=cfg.llm.api_key,
                            timeout=120, max_retries=0)
        # deepseek official models are reasoning-first; disable the thinking
        # phase so `content` and content-token logprobs populate normally
        # (the posterior probe reads the FIRST content token, which a
        # reasoning model otherwise never emits inside an 8-token budget).
        llm.extra_body = {"reasoning_effort": "none"}
        attach_usage(llm, role, events, lock)

    all_tools = build_default_tools(
        rag.retriever, memo_root="", default_collection=COLLECTION,
        schema_profile=rag.schema_profile, use_subagents=False, llm=rag.llm,
        default_with_rerank=cfg.retrieval.with_reranker,
        rewriter=rag.rewriter,
        search_return_top=manifest["search_return_top"])
    kept = [t for t in all_tools if t.name in {"search", "rewrite"}]
    if not any(t.name == "search" for t in kept):
        raise RuntimeError("search tool missing")
    has_rewrite = any(t.name == "rewrite" for t in kept)
    registry = ToolRegistry()
    registry.register_all(kept)
    system_prompt = build_system_prompt(rag.schema_profile,
                                        use_subagents=False)
    if manifest.get("kind_filter_hint"):
        system_prompt += KIND_FILTER_HINT

    (out_dir / "results").mkdir(parents=True, exist_ok=True)
    _save(out_dir / "status.json", {"state": "running", "qids": len(qids)})

    arms = arms or ARMS
    for qi, qid in enumerate(qids, 1):
        qa_row = qa[qid]
        for arm in arms:
            dest = out_dir / "results" / f"{qid}__{arm}.json"
            if dest.exists():
                continue
            events.clear()
            ctrl_on = arm == "ctrl_on"
            force_search = arm == "full_rounds"
            conf = ConfidenceConfig(
                enabled=ctrl_on, top_b=TOP_B,
                predictor_ckpt=cfg.confidence.predictor_ckpt,
                fusion_ckpt=cfg.confidence.fusion_ckpt,
                bocpd_hazard=cfg.confidence.bocpd_hazard,
                r_min=cfg.confidence.r_min, z_star=cfg.confidence.z_star,
                z_star_override=z_star,
                epsilon=cfg.confidence.epsilon, delta=cfg.confidence.delta)
            agent = RagAgent(
                llm=rag.llm, tools=registry, max_steps=T_MAX,
                system_prompt=system_prompt,
                parallel_tool_calls=cfg.agent.parallel_tool_calls,
                max_parallel=cfg.agent.max_parallel, route=False,
                confidence=conf, collect_logprobs_only=not ctrl_on,
                force_first_tool="rewrite" if has_rewrite else None,
                force_search_after_first=force_search,
                posterior_llm=probe)

            record = {"qid": qid, "arm": arm, "ctrl_on": ctrl_on,
                      "question": qa_row["question"],
                      "gold_ans": qa_row["answer"],
                      "scorer": qa_row["scorer"],
                      "n_hops": qa_row["n_hops"],
                      "task_class": qa_row["task_class"],
                      "gold_sp": gold_blocks.get(qid) or []}
            print(f"START {qi}/{len(qids)} {qid} {arm}", flush=True)
            t0 = time.monotonic()
            try:
                result = agent.run(qa_row["question"], max_steps=T_MAX,
                                   session_id=f"{qid}-{arm}")
                long_answer = result.get("answer") or ""
                doc_ids: list[str] = []
                for hit in result.get("observed_hits") or []:
                    did = hit.get("doc_id")
                    if did and did not in doc_ids:
                        doc_ids.append(did)
                contents = fetch_blocks(rag.store, doc_ids[:SP_MAX_BLOCKS])
                candidates = [{"block_id": d, "content": contents[d]}
                              for d in doc_ids[:SP_MAX_BLOCKS] if d in contents]

                short_ans, short_err = extract_short(
                    aux, qa_row["question"], long_answer)
                pred_sp, sp_err = pick_supporting(
                    aux, qa_row["question"], short_ans, candidates)

                record.update(
                    status="ok", pred_ans=long_answer, short_ans=short_ans,
                    short_error=short_err, pred_sp=pred_sp, sp_error=sp_err,
                    rounds=result.get("steps"),
                    terminated=result.get("terminated"),
                    doc_ids=doc_ids, n_candidates=len(candidates),
                    visible_gold_block_hit=bool(
                        set(doc_ids) & set(record["gold_sp"])),
                    **score_row(short_ans, qa_row, pred_sp, record["gold_sp"],
                                long_answer=long_answer))
            except Exception as exc:
                record.update(status="error", error=f"{type(exc).__name__}: {exc}",
                              pred_ans="", short_ans="", pred_sp=[])
            record["elapsed_s"] = time.monotonic() - t0
            record["usage_events"] = list(events)
            _save(dest, record)

    _save(out_dir / "status.json", {"state": "complete", "qids": len(qids)})
    return summarize(out_dir)


def summarize(out_dir: Path) -> int:
    # Merge shard subdirectories (shard_XX_of_YY/results) if present, else read
    # the top-level results dir. Deduplicate by filename so a re-run merge is
    # idempotent.
    result_dirs = [out_dir / "results"]
    shard_dirs = sorted((out_dir).glob("shard_*_of_*"))
    if shard_dirs:
        result_dirs = [d / "results" for d in shard_dirs if (d / "results").is_dir()]
    seen: set[str] = set()
    rows: list[dict] = []
    for rd in result_dirs:
        for p in sorted(rd.glob("*.json")):
            if p.name in seen:
                continue
            seen.add(p.name)
            rows.append(json.loads(p.read_text()))
    ok = [r for r in rows if r.get("status") == "ok"]
    arms_present = sorted({r["arm"] for r in ok})
    metrics = ("contain_acc", "contain_acc_short", "extractor_lost_answer",
               "typed_score", "typed_em", "ans_em", "ans_f1", "sp_em", "sp_f1",
               "joint_em", "joint_f1", "rounds", "elapsed_s")

    def agg(subset: list[dict]) -> dict:
        n = len(subset)
        if not n:
            return {"n": 0}
        out = {"n": n}
        for m in metrics:
            vals = [float(r.get(m) or 0) for r in subset]
            out[f"mean_{m}"] = sum(vals) / n
        out["visible_gold_block_hit_rate"] = (
            sum(1 for r in subset if r.get("visible_gold_block_hit")) / n)
        out["terminated"] = dict(Counter(str(r.get("terminated"))
                                         for r in subset))
        tok = 0
        for r in subset:
            for e in r.get("usage_events") or []:
                tok += ((e.get("usage") or {}).get("total_tokens") or 0)
        out["mean_total_tokens"] = tok / n
        return out

    summary = {
        "protocol": "besiii_parity_hotpot_v1",
        "n_records": len(rows),
        "n_ok": len(ok),
        "statuses": dict(Counter(r.get("status") for r in rows)),
        "arms": {a: agg([r for r in ok if r["arm"] == a]) for a in arms_present},
        "by_scorer": {s: {a: agg([r for r in ok if r["arm"] == a
                                  and r["scorer"] == s]) for a in arms_present}
                      for s in sorted({r["scorer"] for r in ok})},
    }
    _save(out_dir / "summary.json", summary)
    print("\n=== besiii_parity_hotpot_v1 ===")
    for a in arms_present:
        m = summary["arms"][a]
        if not m.get("n"):
            continue
        print(f"{a:<9} n={m['n']:<3} ContainAcc={m['mean_contain_acc']:.4f} "
              f"(short={m['mean_contain_acc_short']:.4f}, "
              f"extractorLost={m['mean_extractor_lost_answer']:.4f}) "
              f"typed={m['mean_typed_score']:.4f} "
              f"AnsEM={100*m['mean_ans_em']:.2f} AnsF1={100*m['mean_ans_f1']:.2f} "
              f"SupEM={100*m['mean_sp_em']:.2f} SupF1={100*m['mean_sp_f1']:.2f} "
              f"JointF1={100*m['mean_joint_f1']:.2f} "
              f"rounds={m['mean_rounds']:.2f} tok={m['mean_total_tokens']:.0f}")
    return 0


def main() -> int:
    args = _parse_args()
    out_dir: Path = args.out_dir
    if args.merge:
        return summarize(out_dir)
    if args.prepare:
        if args.qid_source == "file":
            if not args.qids_file:
                print("--qid-source file requires --qids-file", file=sys.stderr)
                return 2
            qids = json.loads(args.qids_file.read_text(encoding="utf-8"))
        elif args.qid_source == "qa120":
            qids = qa120_qids()
        elif args.qid_source == "full":
            qids = full_qids()
        else:
            qids = SMOKE_QIDS
        _save(out_dir / "manifest.json",
              build_manifest(qids, args.qid_source, args.variant))
        print(f"prepared {out_dir / 'manifest.json'} with {len(qids)} qids "
              f"({args.qid_source}, variant={args.variant})")
        return 0
    if args.validate:
        return validate(out_dir)
    if args.run:
        if not args.allow_endpoint:
            print("--run requires --allow-endpoint", file=sys.stderr)
            return 2
        if (args.shard_id is None) != (args.num_shards is None):
            print("--shard-id and --num-shards must be given together",
                  file=sys.stderr)
            return 2
        arms = [a.strip() for a in args.arms.split(",")] if args.arms else None
        return run(out_dir, args.limit, args.shard_id, args.num_shards,
                   args.z_star, arms)
    print("nothing to do: pass --prepare, --validate, --run or --merge",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
