"""`hepuke` command-line interface.

Subcommands:
    hepuke ingest --source memo_md --root PATH --collection NAME [--limit N]
    hepuke search "query" --collection NAME [--top-k 5]
    hepuke grep "pattern" [--root PATH] [--doc-id ID]
    hepuke read DOC_ID [--collection NAME]
    hepuke ask "question" [--collection NAME] [--max-steps N]
    hepuke serve  [--host 0.0.0.0 --port 42796]
    hepuke collections list|create|drop [NAME]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional

from hepuke.config import load_config
from hepuke.document import Document


def _rag_from_env(cfg_path: Optional[str] = None):
    from hepuke.core.index import HepUKE
    cfg = load_config(cfg_path or os.environ.get("HEPUKE_CONFIG"))
    return cfg, HepUKE(cfg)


# ---------------------------------------------------------------------------
# ingest
# ---------------------------------------------------------------------------

def _cmd_ingest(args: argparse.Namespace) -> int:
    cfg, rag = _rag_from_env(args.config)
    rag.connect_collection(args.collection)

    from hepuke.core.chunker import (
        PaperSectionChunker, SecondaryChunker, make_embed_text,
        clean_content_for_store,
    )
    parser = PaperSectionChunker()
    secondary = SecondaryChunker()

    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        print(f"[ingest] root not found: {root}", file=sys.stderr)
        return 2

    if args.source == "memo_md":
        files = sorted(root.rglob("*.md"))
    else:
        print(f"[ingest] unsupported --source {args.source!r}", file=sys.stderr)
        return 2

    if args.limit:
        files = files[: args.limit]
    print(f"[ingest] {len(files)} file(s), collection={args.collection}")

    total_rows = 0
    for i, path in enumerate(files, 1):
        text = path.read_text(encoding="utf-8", errors="ignore")
        doc_id = path.stem
        sections = parser.parse(text, doc_id=doc_id, title=doc_id)
        if not sections:
            print(f"[ingest] {path.name}: no sections; falling back to whole-doc")
            sections = [{
                "doc_id": doc_id, "section_type": "body", "section_number": "",
                "section_title": doc_id, "content": text, "title": doc_id,
            }]
        v2 = secondary.run(sections)

        parent_rows: dict[str, dict] = {}
        child_rows: list[dict] = []
        for j, ch in enumerate(v2):
            parent_id = f"{doc_id}::{ch.get('section_type', 'body')}::{ch.get('section_number', '0') or '0'}"
            child = {
                "doc_id": doc_id,
                "chunk_type": "child",
                "parent_id": parent_id,
                "section_type": ch.get("section_type", "body"),
                "section_number": str(ch.get("section_number", "")),
                "section_title": ch.get("section_title", ""),
                "title": ch.get("title", doc_id),
                "authors": ch.get("authors", ""),
                "date": ch.get("date", ""),
                "source_url": str(path.relative_to(root)),
                "content": clean_content_for_store(ch.get("content", "")) or ch.get("content", ""),
                "assets": json.dumps(ch.get("assets", []), ensure_ascii=False)
                if ch.get("assets") else "",
                "metadata_json": "",
                "embed_text": make_embed_text(
                    ch, prefix=args.embed_prefix.format(**ch) if args.embed_prefix else "",
                ),
            }
            child_rows.append(child)
            if parent_id not in parent_rows:
                parent_rows[parent_id] = {
                    **child,
                    "chunk_type": "parent",
                    "content": ch.get("content", ""),
                    "embed_text": make_embed_text(ch, prefix=args.embed_prefix),
                }
        rows = child_rows + list(parent_rows.values())
        n = rag.retriever.upsert_chunks(args.collection, rows)
        total_rows += n
        if i % 10 == 0 or i == len(files):
            print(f"[ingest] {i}/{len(files)} — {path.name} → {n} rows (total {total_rows})")

    print(f"[ingest] done — {total_rows} rows across {len(files)} docs")
    return 0


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

def _cmd_search(args: argparse.Namespace) -> int:
    cfg, rag = _rag_from_env(args.config)
    rag.connect_collection(args.collection)
    hits = rag.retrieve(args.query, similarity_top_k=args.top_k)
    for i, h in enumerate(hits, 1):
        print(f"[{i}] score={h.get('score'):.4f}  doc_id={h.get('doc_id','?')}  "
              f"section={h.get('section_type','?')}  {h.get('section_title','')}")
        content = (h.get("content") or "").replace("\n", " ")
        print(f"     {content[:180]}")
    return 0


# ---------------------------------------------------------------------------
# grep
# ---------------------------------------------------------------------------

def _cmd_grep(args: argparse.Namespace) -> int:
    from hepuke.tools import memo_tools
    cfg = load_config(args.config or os.environ.get("HEPUKE_CONFIG"))
    root = args.root or cfg.wiki.memo_md_root
    doc_ids = [args.doc_id] if args.doc_id else None
    hits = memo_tools.grep(
        pattern=args.pattern,
        path_root=root,
        doc_ids=doc_ids,
        context=args.context,
        ignore_case=args.ignore_case,
        max_results=args.max_results,
    )
    for h in hits:
        tag = "M" if h["is_match"] else " "
        print(f"{tag} {h['file']}:{h['line']}: {h['text']}")
    return 0


# ---------------------------------------------------------------------------
# read
# ---------------------------------------------------------------------------

def _cmd_read(args: argparse.Namespace) -> int:
    from hepuke.tools import memo_tools
    cfg = load_config(args.config or os.environ.get("HEPUKE_CONFIG"))
    if args.collection:
        _, rag = _rag_from_env(args.config)
        rag.connect_collection(args.collection)
        r = memo_tools.read_memo(
            args.doc_id, retriever=rag.retriever, collection=args.collection,
            path_root=cfg.wiki.memo_md_root,
        )
    else:
        r = memo_tools.read_memo(args.doc_id, path_root=cfg.wiki.memo_md_root)
    print(f"[source={r.get('source')}]")
    print(r.get("text", ""))
    return 0


# ---------------------------------------------------------------------------
# ask (agent loop)
# ---------------------------------------------------------------------------

def _cmd_ask(args: argparse.Namespace) -> int:
    cfg, rag = _rag_from_env(args.config)
    collection = args.collection or cfg.milvus.default_collection
    rag.connect_collection(collection)

    from hepuke.agent import RagAgent, ToolRegistry
    from hepuke.agent.defaults import build_default_tools

    registry = ToolRegistry()
    registry.register_all(
        build_default_tools(
            rag.retriever, memo_root=cfg.wiki.memo_md_root,
            default_collection=collection,
            default_with_rerank=cfg.retrieval.with_reranker,
        )
    )
    agent = RagAgent(
        llm=rag.llm, tools=registry,
        max_steps=args.max_steps or cfg.agent.max_steps,
    )
    result = agent.run(args.question)
    print("=" * 60)
    print("ANSWER:")
    print(result["answer"])
    print("=" * 60)
    print(f"steps: {result['steps']}")
    for t in result["trace"]:
        print(f"  step={t['step']} tool={t['tool']} latency_ms={t['latency_ms']}"
              f"{' [ERR: ' + t['error'] + ']' if t.get('error') else ''}")
    return 0


# ---------------------------------------------------------------------------
# serve
# ---------------------------------------------------------------------------

def _cmd_serve(args: argparse.Namespace) -> int:
    import uvicorn
    from hepuke.service.fastapi_app import create_app
    cfg = load_config(args.config or os.environ.get("HEPUKE_CONFIG"))
    app = create_app(args.config)
    uvicorn.run(
        app, host=args.host or cfg.service.host,
        port=args.port or cfg.service.port,
    )
    return 0


# ---------------------------------------------------------------------------
# collections
# ---------------------------------------------------------------------------

def _cmd_collections(args: argparse.Namespace) -> int:
    _, rag = _rag_from_env(args.config)
    if args.action == "list":
        for c in rag.get_collections():
            print(c)
    elif args.action == "create":
        rag.create_collection(args.name)
        print(f"created {args.name}")
    elif args.action == "drop":
        rag.drop_collection(args.name)
        print(f"dropped {args.name}")
    return 0


# ---------------------------------------------------------------------------
# bench
# ---------------------------------------------------------------------------

def _cmd_bench(args: argparse.Namespace) -> int:
    from benchmarks import loaders as bench_loaders
    from benchmarks.runner import replay, run_suite

    if args.action == "list":
        available = bench_loaders.list_available(args.data_root)
        if not available:
            print(f"[bench] no suites in {args.data_root}", file=sys.stderr)
            return 1
        for name in available:
            print(name)
        return 0

    if args.action == "replay":
        replay(
            trace_file=args.trace,
            gold_file=args.gold,
            out_dir=args.out,
            ks=tuple(args.ks),
            suite=args.suite or "",
            config=args.name or "replay",
            run_id=args.run_id or None,
        )
        return 0

    # default: run
    config_path = args.config_file
    if not Path(config_path).exists():
        # allow bare name → benchmarks/configs/<name>.yaml
        alt = Path("benchmarks/configs") / f"{config_path}.yaml"
        if alt.exists():
            config_path = str(alt)
        else:
            print(f"[bench] config not found: {args.config_file}", file=sys.stderr)
            return 2
    run_suite(
        suite=args.suite,
        config_file=config_path,
        out_dir=args.out,
        limit=args.limit or None,
        mode=args.mode,
        data_root=args.data_root,
        hai_config_path=args.config,
        trace_backend=args.trace_backend,
        trace_db_path=args.trace_db or None,
    )
    return 0


# ---------------------------------------------------------------------------
# trace (inspect SQLite trace DB)
# ---------------------------------------------------------------------------

def _cmd_trace(args: argparse.Namespace) -> int:
    from hepuke.observability.sinks import (
        is_sqlite_path,
        iter_sqlite_events,
        known_run_ids,
    )

    cfg = load_config(args.config or os.environ.get("HEPUKE_CONFIG"))
    db_path = args.db or cfg.trace.db_path
    if not is_sqlite_path(db_path) or not Path(db_path).exists():
        print(f"[trace] not a SQLite DB or missing: {db_path}", file=sys.stderr)
        return 2

    if args.action == "runs":
        runs = known_run_ids(db_path)
        for r in runs:
            print(r)
        return 0

    if args.action == "tail":
        it = iter_sqlite_events(
            db_path,
            run_id=args.run_id or None,
            kind=args.kind or None,
            case_id=args.case_id or None,
        )
        rows = list(it)
        for rec in rows[-args.limit:] if args.limit else rows:
            print(json.dumps(rec, ensure_ascii=False))
        return 0

    return 2


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hepuke", description="HepUKE CLI")
    p.add_argument("--config", help="path to config.yaml (or $HEPUKE_CONFIG)")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("ingest", help="chunk + embed a corpus into a collection")
    s.add_argument("--source", default="memo_md",
                   choices=["memo_md"], help="ingestion source type")
    s.add_argument("--root", required=True, help="root directory of the corpus")
    s.add_argument("--collection", required=True)
    s.add_argument("--limit", type=int, default=0,
                   help="only ingest the first N files (0 = all)")
    s.add_argument("--embed-prefix", default="",
                   help="format string, e.g. '[BESIII Memo | {section_type}]'")
    s.set_defaults(func=_cmd_ingest)

    s = sub.add_parser("search", help="hybrid search a collection")
    s.add_argument("query")
    s.add_argument("--collection", required=True)
    s.add_argument("--top-k", type=int, default=5)
    s.set_defaults(func=_cmd_search)

    s = sub.add_parser("grep", help="regex search over the source markdown")
    s.add_argument("pattern")
    s.add_argument("--root", help="path root (defaults to cfg.wiki.memo_md_root)")
    s.add_argument("--doc-id", help="restrict to one document")
    s.add_argument("--context", type=int, default=3)
    s.add_argument("-i", "--ignore-case", action="store_true")
    s.add_argument("--max-results", type=int, default=200)
    s.set_defaults(func=_cmd_grep)

    s = sub.add_parser("read", help="fetch full text of one document")
    s.add_argument("doc_id")
    s.add_argument("--collection", help="fetch from Milvus rather than disk")
    s.set_defaults(func=_cmd_read)

    s = sub.add_parser("ask", help="agent loop: grep + RAG + LLM")
    s.add_argument("question")
    s.add_argument("--collection", help="default collection for search")
    s.add_argument("--max-steps", type=int, default=0)
    s.set_defaults(func=_cmd_ask)

    s = sub.add_parser("serve", help="run FastAPI server")
    s.add_argument("--host", default="")
    s.add_argument("--port", type=int, default=0)
    s.set_defaults(func=_cmd_serve)

    s = sub.add_parser("collections", help="manage Milvus collections")
    s.add_argument("action", choices=["list", "create", "drop"])
    s.add_argument("name", nargs="?", default="")
    s.set_defaults(func=_cmd_collections)

    s = sub.add_parser("bench", help="benchmark harness: run / replay / list suites")
    bench_sub = s.add_subparsers(dest="action", required=True)

    b_run = bench_sub.add_parser("run", help="run one suite × one config")
    b_run.add_argument("suite", help="suite name (matches benchmarks/data/<name>.jsonl)")
    b_run.add_argument("--config-file", default="baseline",
                       help="YAML path OR bare name → benchmarks/configs/<name>.yaml")
    b_run.add_argument("--limit", type=int, default=0)
    b_run.add_argument("--mode", choices=["retrieval", "answer", "both"], default="retrieval")
    b_run.add_argument("--data-root", default="benchmarks/data")
    b_run.add_argument("--out", default="benchmarks/reports")
    b_run.add_argument("--trace-backend", choices=["sqlite", "jsonl", "both"],
                       default="sqlite",
                       help="where trace records go (default: sqlite persistent DB)")
    b_run.add_argument("--trace-db", default="",
                       help="override SQLite DB path (else uses config.trace.db_path)")
    b_run.set_defaults(action="run", func=_cmd_bench)

    b_list = bench_sub.add_parser("list", help="list available suites")
    b_list.add_argument("--data-root", default="benchmarks/data")
    b_list.set_defaults(action="list", func=_cmd_bench)

    b_rep = bench_sub.add_parser("replay", help="recompute metrics from trace + gold")
    b_rep.add_argument("--trace", required=True,
                       help="JSONL file or SQLite DB (.db/.sqlite)")
    b_rep.add_argument("--gold", required=True)
    b_rep.add_argument("--out", default="benchmarks/reports")
    b_rep.add_argument("--ks", type=int, nargs="+", default=[1, 3, 5, 10])
    b_rep.add_argument("--suite", default="")
    b_rep.add_argument("--name", default="replay")
    b_rep.add_argument("--run-id", default="",
                       help="filter SQLite DB by run_id (ignored for JSONL)")
    b_rep.set_defaults(action="replay", func=_cmd_bench)

    s = sub.add_parser("trace", help="inspect the persistent SQLite trace DB")
    trace_sub = s.add_subparsers(dest="action", required=True)

    t_runs = trace_sub.add_parser("runs", help="list run_ids stored in the DB")
    t_runs.add_argument("--db", default="",
                        help="override DB path (else uses config.trace.db_path)")
    t_runs.set_defaults(action="runs", func=_cmd_trace)

    t_tail = trace_sub.add_parser("tail", help="print records (optionally filtered)")
    t_tail.add_argument("--db", default="")
    t_tail.add_argument("--run-id", default="")
    t_tail.add_argument("--kind", default="")
    t_tail.add_argument("--case-id", default="")
    t_tail.add_argument("--limit", type=int, default=50,
                        help="print at most the last N matching records")
    t_tail.set_defaults(action="tail", func=_cmd_trace)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cmd == "collections" and args.action in ("create", "drop") and not args.name:
        parser.error(f"{args.action} requires a collection name")
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
