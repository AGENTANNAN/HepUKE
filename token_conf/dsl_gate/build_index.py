"""Build a Milvus collection of BESIII DSL descriptions for the entropy-gating experiment.

Each row = one gold description in ``token_conf/batch/descriptions/<stem>.txt``,
paired with the path to its gold ``token_conf/batch/runs/<stem>.rb``. Rows are
embedded with the local BGE-M3 model (dense 1024d + sparse) so the collection
can be queried through the existing hybrid retriever.

New collection: ``besiii_dsl_desc_v1``. Independent of any pre-existing index.

Usage (single-GPU):
    python token_conf/dsl_gate/build_index.py --drop-existing
Sanity query:
    python token_conf/dsl_gate/build_index.py --sanity-query "psi(2S) to pi+ pi- J/psi radiative"
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hepuke import HepUKE, load_config  # noqa: E402


DESC_DIR = ROOT / "data" / "corpus" / "descriptions"
RB_DIR = ROOT / "data" / "corpus" / "generated_dsl"


def _stems() -> list[str]:
    return sorted(p.stem for p in DESC_DIR.glob("*.txt"))


def _truncate_utf8(s: str, max_bytes: int) -> str:
    if not s:
        return ""
    b = s.encode("utf-8", "surrogatepass")
    return s if len(b) <= max_bytes else b[:max_bytes].decode("utf-8", "ignore")


def _row(stem: str, description: str) -> dict:
    rb_path = str(RB_DIR / f"{stem}.rb")
    meta = {"stem": stem, "rb_path": rb_path, "index_version": "besiii_dsl_desc_v1"}
    parent_key = f"{stem}::desc::0"
    return {
        "doc_id": _truncate_utf8(stem, 128),
        "chunk_type": "child",
        "parent_id": _truncate_utf8(parent_key, 128),
        "section_type": "description",
        "section_number": "0",
        "section_title": "",
        "title": _truncate_utf8(stem, 512),
        "authors": "",
        "date": "",
        "source_url": _truncate_utf8(f"arxiv://{stem}", 256),
        "content": _truncate_utf8(description, 65535),
        "assets": "",
        "metadata_json": _truncate_utf8(json.dumps(meta, ensure_ascii=False), 65535),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--collection", default="besiii_dsl_desc_v1")
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--drop-existing", action="store_true")
    ap.add_argument("--limit", type=int, default=None,
                    help="only ingest the first N stems (smoke)")
    ap.add_argument("--sanity-query", default=None,
                    help="hybrid_search top-3 after ingest")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s")
    log = logging.getLogger("dsl_gate.build_index")

    stems = _stems()
    if args.limit:
        stems = stems[: args.limit]
    if not stems:
        log.error("no descriptions found under %s", DESC_DIR)
        return 2
    log.info("found %d descriptions", len(stems))

    cfg = load_config()
    cfg.milvus.default_collection = args.collection
    rag = HepUKE(cfg)

    existing = rag.get_collections()
    if args.drop_existing and args.collection in existing:
        log.info("dropping existing collection %s", args.collection)
        rag.drop_collection(args.collection)
    if args.collection not in rag.get_collections():
        log.info("creating collection %s", args.collection)
        rag.create_collection(args.collection)
    rag.connect_collection(args.collection)
    retriever = rag.retriever

    t0 = time.time()
    pending: list[dict] = []
    n_missing_rb = 0

    def _drain(force: bool = False) -> None:
        nonlocal pending
        if not pending:
            return
        if not force and len(pending) < args.batch_size:
            return
        retriever.upsert_chunks(
            collection=args.collection,
            chunks=pending,
            embed_text_key="content",
            content_key="content",
        )
        pending = []

    for i, stem in enumerate(stems, 1):
        desc_path = DESC_DIR / f"{stem}.txt"
        rb_path = RB_DIR / f"{stem}.rb"
        if not rb_path.exists():
            n_missing_rb += 1
        description = desc_path.read_text(encoding="utf-8")
        pending.append(_row(stem, description))
        _drain(force=False)
        if i % 100 == 0 or i == len(stems):
            log.info("[%d/%d] rows staged | %.1fs", i, len(stems), time.time() - t0)

    _drain(force=True)

    log.info("=== done ===")
    log.info("  rows      : %d", len(stems))
    log.info("  missing rb: %d", n_missing_rb)
    log.info("  wall      : %.1fs", time.time() - t0)

    if args.sanity_query:
        hits = retriever.search(args.sanity_query, args.collection, top_k=3)
        log.info("=== sanity_query %r ===", args.sanity_query)
        for h in hits:
            log.info("  doc=%s content[:160]=%r",
                     h.get("doc_id"), (h.get("content") or "")[:160])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
