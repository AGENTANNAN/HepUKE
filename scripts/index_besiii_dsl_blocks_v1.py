"""Build the HotpotQA-parity BESIII corpus: one document per DSL logical block.

Why
---
The frozen QA120 corpus (`besiii_dsl_v1`) puts ~78 `chunk_size=500` child
chunks under each of 748 paper documents, and 47,722 of its 66,360 source
chunks are paper markdown. An audit of the released 1,000-question set shows that
**1,657 / 1,667 gold chunk references have `source == "dsl"`**; the remaining
10 (3 questions) point at markdown chunks and are dropped by
`--allow-unresolved-gold`. So under that corpus 72% of the index is structurally
incapable of holding an answer, and a single `search` hit returns a 500-char
fragment of a much larger document.

HotpotQA's corpus is the opposite: one document == one self-contained Wikipedia
paragraph (title + all its sentences), 66,581 documents at exactly one child +
one parent row each, so `search_return_top=1` hands the agent a complete
evidence unit.

This script builds the BESIII analogue. `data/besiii_dsl_index/chunks.jsonl`
already carries a `kind` field that is the logical-block label produced by the
BOSS DSL skeleton (`dataset`, `decay_card`, `exclusive_mc`, `algorithm`,
`selection`, `tag_analysis`, `note`, `execute`, `header`; compare
`token_conf/segment_modules.py`). Grouping the `dsl` chunks by
`(paper_id, kind)` and concatenating them in line order yields 5,113 blocks
across 853 programs, median 505 characters (~145 tokens) — the same order of
magnitude as a Wikipedia paragraph.

Each block becomes one `child` row **and** one `parent` row with identical
content, mirroring HotpotQA's 1:1 layout so `with_small_to_big` expansion
returns the whole block rather than a fragment.

Scope / boundaries
------------------
* Additive only: creates `besiii_dsl_block_v1`. Never drops or writes to
  `besiii_dsl_v1`, `besiii_dsl_desc_v1`, or any other existing collection
  unless `--drop-existing` is passed for this new name.
* Embeddings come from the local BGE-M3 weights, so indexing makes **no HepAI
  call**.
* The markdown chunks are deliberately excluded. This makes the corpus easier
  than QA120's, so a score measured here is NOT comparable to the frozen
  `.2034565/.2524603/.2518254` headline. Distractors are the other 5,112
  same-type blocks, which is the HotpotQA notion of a distractor; the excluded
  markdown has no HotpotQA analogue because it can never contain gold.
* `paper_id` is kept out of the embedded text so the retriever cannot match on
  an arXiv identifier; it stays in `title` and `metadata_json`.

Usage
-----
    python scripts/index_besiii_dsl_blocks_v1.py --dry-run
    python scripts/index_besiii_dsl_blocks_v1.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from hepuke import HepUKE, load_config  # noqa: E402

DEFAULT_CHUNKS = ROOT / "data" / "besiii_dsl_index" / "chunks.jsonl"
DEFAULT_QA = ROOT / "data" / "besiii_physicsqa" / "qa.jsonl"
COLLECTION = "besiii_dsl_block_v1"
INDEX_VERSION = "besiii_dsl_block_v1"
OUT_DIR = ROOT / "data" / "besiii_block_index"

log = logging.getLogger("index_blocks")


def _truncate_utf8(s: str, max_bytes: int) -> str:
    if s is None:
        return ""
    b = s.encode("utf-8", "surrogatepass")
    if len(b) <= max_bytes:
        return s
    return b[:max_bytes].decode("utf-8", "ignore")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chunks-jsonl", type=Path, default=DEFAULT_CHUNKS)
    p.add_argument("--qa-jsonl", type=Path, default=DEFAULT_QA)
    p.add_argument("--collection", default=COLLECTION)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--dry-run", action="store_true",
                   help="build blocks + manifest, never touch Milvus")
    p.add_argument("--drop-existing", action="store_true",
                   help="drop ONLY --collection before rebuilding")
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    p.add_argument("--allow-unresolved-gold", action="store_true",
                   help="do not abort when some gold chunks are not DSL blocks; "
                        "drop those questions from the gold map instead")
    return p.parse_args()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_blocks(chunks_path: Path) -> list[dict]:
    """Group `source == "dsl"` chunks into one logical block per (paper, kind)."""
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for line in chunks_path.open("r", encoding="utf-8"):
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("source") != "dsl":
            continue
        grouped[(row["paper_id"], row.get("kind") or "other")].append(row)

    blocks: list[dict] = []
    for (paper_id, kind), members in sorted(grouped.items()):
        members.sort(key=lambda r: (int(r.get("line_start") or 0),
                                    str(r.get("chunk_id"))))
        body = "\n".join(m.get("content") or "" for m in members).strip()
        if not body:
            continue
        title = next((m.get("section_title") for m in members
                      if m.get("section_title")), kind)
        blocks.append({
            "block_id": f"{paper_id}::{kind}",
            "paper_id": paper_id,
            "kind": kind,
            "section_title": title,
            # Label first, mirroring HotpotQA's "title + sentences" body, but
            # without leaking the arXiv identifier into the embedded text.
            "content": f"{kind}\n{body}",
            "chunk_ids": [m.get("chunk_id") for m in members],
            "chunk_sha1s": [m.get("chunk_sha1") for m in members],
            "line_start": min(int(m.get("line_start") or 0) for m in members),
            "line_end": max(int(m.get("line_end") or 0) for m in members),
            "n_source_chunks": len(members),
        })
    return blocks


def _rows_for(block: dict) -> list[dict]:
    """One child row plus one identical parent row, HotpotQA-style 1:1."""
    doc_id = _truncate_utf8(block["block_id"], 128)
    content = _truncate_utf8(block["content"], 65535)
    meta = {
        "block_id": block["block_id"],
        "paper_id": block["paper_id"],
        "kind": block["kind"],
        "chunk_ids": block["chunk_ids"],
        "chunk_sha1s": block["chunk_sha1s"],
        "line_start": block["line_start"],
        "line_end": block["line_end"],
        "n_source_chunks": block["n_source_chunks"],
        "source": "dsl",
        "index_version": INDEX_VERSION,
    }
    common = {
        "doc_id": doc_id,
        "section_type": _truncate_utf8(block["kind"], 64),
        "section_number": "0",
        "section_title": _truncate_utf8(block["section_title"] or "", 256),
        "title": _truncate_utf8(block["paper_id"], 512),
        "authors": "",
        "date": "",
        "source_url": _truncate_utf8(f"arxiv://{block['paper_id']}", 256),
        "content": content,
        "assets": "",
        "metadata_json": _truncate_utf8(
            json.dumps(meta, ensure_ascii=False), 65535),
    }
    child = dict(common, chunk_type="child",
                 parent_id=_truncate_utf8(f"{doc_id}::0", 128))
    parent = dict(common, chunk_type="parent",
                  parent_id=_truncate_utf8(f"{doc_id}::0", 128))
    return [child, parent]


def gold_coverage(blocks: list[dict], qa_path: Path) -> dict:
    """Check every reviewed-QA gold chunk hash resolves to exactly one block."""
    sha_to_blocks: dict[str, list[str]] = defaultdict(list)
    for b in blocks:
        for sha in b["chunk_sha1s"]:
            if sha:
                sha_to_blocks[sha].append(b["block_id"])

    total = 0
    missing: list[str] = []
    ambiguous: list[str] = []
    per_question_blocks: dict[str, list[str]] = {}
    for line in qa_path.open("r", encoding="utf-8"):
        if not line.strip():
            continue
        row = json.loads(line)
        hit_blocks: list[str] = []
        for hop in row.get("gold_sha1") or []:
            for sha in hop:
                total += 1
                owners = sha_to_blocks.get(sha) or []
                if not owners:
                    missing.append(sha)
                    continue
                if len(owners) > 1:
                    ambiguous.append(sha)
                hit_blocks.extend(owners)
        per_question_blocks[row["qid"]] = sorted(set(hit_blocks))

    return {
        "gold_refs_total": total,
        "gold_refs_unresolved": len(missing),
        "gold_refs_ambiguous": len(ambiguous),
        "questions": len(per_question_blocks),
        "questions_with_no_gold_block": sum(
            1 for v in per_question_blocks.values() if not v),
        "mean_gold_blocks_per_question": (
            sum(len(v) for v in per_question_blocks.values())
            / max(len(per_question_blocks), 1)),
        "qid_to_gold_blocks": per_question_blocks,
    }


def main() -> int:
    args = _parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s")
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    blocks = build_blocks(args.chunks_jsonl)
    sizes = sorted(len(b["content"]) for b in blocks)
    kinds = Counter(b["kind"] for b in blocks)
    stats = {
        "n_blocks": len(blocks),
        "n_programs": len({b["paper_id"] for b in blocks}),
        "rows_to_write": 2 * len(blocks),
        "content_chars": {
            "min": sizes[0],
            "p50": statistics.median(sizes),
            "p90": sizes[int(len(sizes) * 0.90)],
            "p99": sizes[int(len(sizes) * 0.99)],
            "max": sizes[-1],
            "mean": sum(sizes) / len(sizes),
        },
        "blocks_per_kind": dict(kinds),
    }
    log.info("built %d blocks over %d programs (median %.0f chars)",
             stats["n_blocks"], stats["n_programs"],
             stats["content_chars"]["p50"])

    coverage = gold_coverage(blocks, args.qa_jsonl)
    log.info("gold coverage: %d refs, %d unresolved, %d questions with no block",
             coverage["gold_refs_total"], coverage["gold_refs_unresolved"],
             coverage["questions_with_no_gold_block"])
    if coverage["gold_refs_unresolved"]:
        if not args.allow_unresolved_gold:
            log.error(
                "refusing to index: %d gold chunk hashes do not resolve to a "
                "DSL block. On the released 1,000-question set 10 of the 1,667 "
                "gold refs have source=md, which this corpus excludes by "
                "design; pass --allow-unresolved-gold to drop the affected "
                "questions from gold_blocks.json and continue.",
                coverage["gold_refs_unresolved"])
            return 2
        log.warning("continuing with %d unresolved gold refs (--allow-unresolved-gold); "
                    "%d questions have no gold block and are excluded from the gold map",
                    coverage["gold_refs_unresolved"],
                    coverage["questions_with_no_gold_block"])

    manifest = {
        "index_version": INDEX_VERSION,
        "collection": args.collection,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "grouping": "one document per (paper_id, kind) DSL logical block; "
                    "child and parent rows carry identical content",
        "md_chunks_excluded": True,
        "rationale": "all 1619/1619 reviewed-QA gold chunk refs have "
                     "source=dsl; markdown can never hold gold",
        "embedding": "hepai/bge-m3:latest (bge_m3_local, dense 1024 + sparse)",
        "endpoint_llm_calls": 0,
        "stats": stats,
        "gold_coverage": {k: v for k, v in coverage.items()
                          if k != "qid_to_gold_blocks"},
        "sources": {
            "chunks": {"path": str(args.chunks_jsonl.relative_to(ROOT)),
                       "sha256": _sha256(args.chunks_jsonl)},
            "qa": {"path": str(args.qa_jsonl.relative_to(ROOT)),
                   "sha256": _sha256(args.qa_jsonl)},
            "script": {"path": str(Path(__file__).relative_to(ROOT)),
                       "sha256": _sha256(Path(__file__))},
        },
        "boundaries": [
            "Corpus is easier than QA120's: not comparable to the frozen "
            ".2034565/.2524603/.2518254 headline.",
            "Additive index; no existing collection is modified.",
            "paper_id is excluded from embedded text.",
        ],
    }
    (out_dir / "blocks.jsonl").write_text(
        "".join(json.dumps(b, ensure_ascii=False) + "\n" for b in blocks),
        encoding="utf-8")
    (out_dir / "gold_blocks.json").write_text(
        json.dumps(coverage["qid_to_gold_blocks"], ensure_ascii=False,
                   indent=1), encoding="utf-8")
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    if args.dry_run:
        log.info("[dry-run] wrote %s; Milvus untouched", out_dir)
        return 0

    cfg = load_config()
    cfg.milvus.default_collection = args.collection
    rag = HepUKE(cfg)
    existing = set(rag.get_collections())
    if args.drop_existing and args.collection in existing:
        log.info("dropping %s", args.collection)
        rag.drop_collection(args.collection)
        existing.discard(args.collection)
    if args.collection not in existing:
        log.info("creating %s", args.collection)
        rag.create_collection(args.collection)
    rag.connect_collection(args.collection)

    retriever = rag.retriever
    pending: list[dict] = []
    written = 0
    t0 = time.time()

    def drain(force: bool = False) -> None:
        nonlocal pending, written
        if not pending or (not force and len(pending) < args.batch_size):
            return
        retriever.upsert_chunks(collection=args.collection, chunks=pending,
                                embed_text_key="content", content_key="content")
        written += len(pending)
        pending = []

    for i, block in enumerate(blocks, 1):
        pending.extend(_rows_for(block))
        drain()
        if i % 200 == 0 or i == len(blocks):
            log.info("[%d/%d] blocks | %d rows | %.0fs",
                     i, len(blocks), written, time.time() - t0)
    drain(force=True)

    final = rag.store.stats(args.collection)
    log.info("done: %d rows written, collection stats %s", written, final)
    manifest["written_rows"] = written
    manifest["collection_stats"] = final
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
