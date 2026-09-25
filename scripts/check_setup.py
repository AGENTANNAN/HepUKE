"""Preflight check — what is installed, what data is present, what to run next.

    python scripts/check_setup.py            # offline checks only
    python scripts/check_setup.py --live     # also probe Milvus and the LLM endpoint

Exit code is 0 when every REQUIRED item passes, 1 otherwise. Optional items and
pipeline stages that have simply not been run yet never fail the check; they are
reported with the command that produces them.
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

OK, WARN, BAD = "  ok  ", " todo ", " MISS "
_failed = False


def line(state: str, name: str, detail: str = "") -> None:
    global _failed
    if state is BAD:
        _failed = True
    print(f"[{state}] {name:<52} {detail}".rstrip())


def section(title: str) -> None:
    print(f"\n── {title} " + "─" * max(0, 58 - len(title)))


# ---------------------------------------------------------------- packages
REQUIRED = [("pymilvus", "vector store"), ("openai", "OpenAI-compatible client"),
            ("pydantic", ""), ("yaml", "pyyaml"), ("numpy", ""), ("scipy", "BOCPD"),
            ("sklearn", "scikit-learn, predictor training"), ("httpx", ""),
            ("requests", ""), ("tenacity", ""), ("networkx", ""), ("dotenv", "python-dotenv")]
OPTIONAL = [("fastapi", "only for the REST service"),
            ("uvicorn", "only for the REST service"),
            ("matplotlib", "only for figure scripts, pip install .[figures]"),
            ("FlagEmbedding", "local BGE-M3 embedding, pip install .[gpu]"),
            ("torch", "local BGE-M3 embedding, pip install .[gpu]"),
            ("hepai", "remote embedding gateway, pip install .[gateway]")]


def check_packages() -> None:
    section("Python packages")
    v = sys.version_info
    line(OK if v >= (3, 10) else BAD, f"python {v.major}.{v.minor}.{v.micro}",
         "" if v >= (3, 10) else "requires >= 3.10")
    for mod, note in REQUIRED:
        try:
            importlib.import_module(mod); line(OK, mod, note)
        except ImportError:
            line(BAD, mod, f"missing — pip install -e .  ({note})" if note else "missing — pip install -e .")
    for mod, note in OPTIONAL:
        try:
            importlib.import_module(mod); line(OK, mod, note)
        except ImportError:
            line(WARN, mod, note)


# ------------------------------------------------------------------- config
def check_config() -> dict:
    section("Configuration")
    cfg_path = ROOT / "config.yaml"
    if not cfg_path.exists():
        line(BAD, "config.yaml", "missing — cp config.example.yaml config.yaml")
        return {}
    import yaml
    cfg = yaml.safe_load(cfg_path.read_text()) or {}
    line(OK, "config.yaml", str(cfg_path))
    if (ROOT / ".env").exists():
        line(OK, ".env", "")
    else:
        line(WARN, ".env", "cp .env.example .env  (or export the variables)")

    import os
    from dotenv import dotenv_values
    env = {**dotenv_values(ROOT / ".env"), **os.environ} if (ROOT / ".env").exists() else dict(os.environ)
    for var, why in [("MILVUS_URI", "vector store"), ("LLM_BASE_URL", "generation"),
                     ("LLM_API_KEY", "generation"), ("LLM_MODEL", "generation")]:
        line(OK if env.get(var) else WARN, var, "" if env.get(var) else f"unset — needed for {why}")
    if not env.get("BGE_M3_PATH"):
        line(WARN, "BGE_M3_PATH", "unset — needed by the local embedding backend")
    return cfg


# --------------------------------------------------------------------- data
def jsonl_len(p: Path) -> int:
    with p.open(encoding="utf-8") as fh:
        return sum(1 for _ in fh)


def check_data() -> None:
    section("Shipped data")
    qa = ROOT / "data" / "besiii_physicsqa" / "qa.jsonl"
    if qa.exists():
        mani = json.loads((qa.parent / "manifest.json").read_text())
        n = jsonl_len(qa)
        line(OK if n == mani["n"] else BAD, "besiii_physicsqa/qa.jsonl",
             f"{n} items" + ("" if n == mani["n"] else f" — manifest says {mani['n']}"))
    else:
        line(BAD, "besiii_physicsqa/qa.jsonl", "missing")

    idx = ROOT / "data" / "besiii_dsl_index"
    for name, expect in (("chunks.jsonl", 66360), ("facts.jsonl", 48865)):
        p = idx / name
        if p.exists():
            n = jsonl_len(p)
            line(OK, f"besiii_dsl_index/{name}", f"{n:,} rows" + ("" if n == expect else f" (expected {expect:,})"))
        else:
            line(BAD, f"besiii_dsl_index/{name}", "missing")

    corpus = ROOT / "data" / "corpus"
    for sub, ext in (("dsl", "*.rb"), ("descriptions", "*.txt"), ("generated_dsl", "*.rb")):
        n = len(list((corpus / sub).glob(ext))) if (corpus / sub).is_dir() else 0
        line(OK if n == 700 else BAD, f"corpus/{sub}", f"{n} files" + ("" if n == 700 else " — expected 700"))


# --------------------------------------------------------- derived artefacts
STAGES = [
    ("data/corpus/splits.json", "python token_conf/dsl_gate/splits.py",
     "train/dev/test + perturbation subsets"),
    ("data/besiii_block_index/gold_blocks.json",
     "python scripts/index_besiii_dsl_blocks_v1.py --dry-run",
     "gold block map, needed by run_besiii_parity_hotpot_v1.py"),
    ("token_conf/dsl_gate/runs_gate/all.jsonl",
     "python token_conf/dsl_gate/generate_candidates.py --split dev",
     "candidate generations, consumed by strategy_matrix / sweet_spot / three_tier_tables / interval_gate"),
    ("token_conf/dsl_gate/runs_gate/all_baseline.jsonl",
     "python token_conf/dsl_gate/generate_baseline.py --split dev",
     "No-RAG baseline generations"),
    ("storage/predictor.pkl", "python scripts/train_predictor.py --data ... --out storage/predictor.pkl",
     "look-ahead predictor for the controller"),
    ("storage/fusion.pkl", "python scripts/grid_search_fusion.py --out storage/fusion.pkl",
     "fusion checkpoint for the controller"),
]


def check_stages() -> None:
    section("Pipeline stages (produced locally, not shipped)")
    for rel, cmd, why in STAGES:
        p = ROOT / rel
        if p.exists():
            line(OK, rel, why)
        else:
            line(WARN, rel, f"not built yet → {cmd}")


# --------------------------------------------------------------- live probes
def check_live(cfg: dict) -> None:
    section("Live endpoints")
    try:
        from hepuke import load_config
        c = load_config()
    except Exception as e:                                    # pragma: no cover
        line(BAD, "load_config()", f"{type(e).__name__}: {e}")
        return
    try:
        from pymilvus import MilvusClient
        cl = MilvusClient(uri=c.milvus.uri, token=c.milvus.token or None,
                          db_name=getattr(c.milvus, "database", None) or "default")
        cols = cl.list_collections()
        line(OK, "milvus", f"{len(cols)} collections: {', '.join(cols[:4])}")
    except Exception as e:
        line(BAD, "milvus", f"{type(e).__name__}: {str(e)[:90]}")
    try:
        from openai import OpenAI
        r = OpenAI(base_url=c.llm.base_url, api_key=c.llm.api_key).chat.completions.create(
            model=c.llm.model, max_tokens=4, messages=[{"role": "user", "content": "ping"}])
        line(OK, "llm", f"{c.llm.model} responded ({r.usage.total_tokens} tokens)")
    except Exception as e:
        line(BAD, "llm", f"{type(e).__name__}: {str(e)[:90]}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--live", action="store_true", help="probe Milvus and the LLM endpoint")
    args = ap.parse_args()

    print(f"repo: {ROOT}")
    check_packages()
    cfg = check_config()
    check_data()
    check_stages()
    if args.live:
        check_live(cfg)

    print()
    if _failed:
        print("FAILED — fix the [ MISS ] rows above before running anything.")
        return 1
    print("OK — required packages and shipped data are in place.")
    print("[ todo ] rows are pipeline stages you have not built yet; that is expected on a fresh clone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
