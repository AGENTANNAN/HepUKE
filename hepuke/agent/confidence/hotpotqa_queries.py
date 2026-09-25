"""HotpotQA loader → queries.jsonl for T7 replay collector.

HotpotQA (Yang et al. 2018) is served as a single JSON array of examples.
We don't depend on the `datasets` library — the file is small enough
(~600 MB fullwiki, ~50 MB distractor dev) to parse with the stdlib.

Each output row is::

    {"id": "<hotpot _id>", "query": "<question text>",
     "answer": "<gold answer>",
     "supporting_facts": [[title, sent_idx], ...],
     "type": "hotpotqa"}

Only ``id`` + ``query`` are consumed by :func:`collect_trajectories`; the
extra fields are kept so downstream evaluation (T8 metrics, T10 fusion
calibration) can pull them from the same JSONL without re-touching the
raw HotpotQA blob.
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


@dataclass
class HotpotExample:
    id: str
    question: str
    answer: str
    supporting_facts: list[list]  # [[title, sent_idx], ...]

    def as_query_row(self) -> dict:
        return {
            "id": self.id,
            "query": self.question,
            "answer": self.answer,
            "supporting_facts": self.supporting_facts,
            "type": "hotpotqa",
        }


def load_hotpotqa_json(path: str | Path) -> list[HotpotExample]:
    """Read the raw HotpotQA JSON (list of dicts) and coerce each example.

    Skips rows missing any of ``_id / question / answer`` — the released
    dev/train sets are well-formed, but we keep the tolerance in case a
    subset or private slice is passed.
    """
    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        raw = json.load(f)
    if not isinstance(raw, list):
        raise ValueError(f"expected a JSON array at {p}, got {type(raw).__name__}")

    out: list[HotpotExample] = []
    for row in raw:
        if not isinstance(row, dict):
            continue
        qid = row.get("_id") or row.get("id")
        q = row.get("question")
        a = row.get("answer", "")
        if not qid or not q:
            continue
        sup = row.get("supporting_facts") or []
        out.append(HotpotExample(
            id=str(qid), question=str(q), answer=str(a),
            supporting_facts=list(sup),
        ))
    return out


def sample_examples(
    examples: list[HotpotExample],
    *,
    n: int,
    seed: int = 0,
) -> list[HotpotExample]:
    """Deterministic random subsample without replacement."""
    if n >= len(examples):
        return list(examples)
    rng = random.Random(seed)
    return rng.sample(examples, n)


def dump_query_rows(examples: Iterable[HotpotExample], out_path: str | Path) -> int:
    """Write the ``as_query_row`` dicts as JSONL (one row per line)."""
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with p.open("w", encoding="utf-8") as f:
        for ex in examples:
            f.write(json.dumps(ex.as_query_row(), ensure_ascii=False))
            f.write("\n")
            n += 1
    return n


def iter_query_rows(path: str | Path) -> Iterator[dict]:
    """Stream the JSONL produced by :func:`dump_query_rows`. Handy for the
    collector's ``queries`` argument when the file is too big to load
    fully into memory."""
    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


__all__ = [
    "HotpotExample",
    "load_hotpotqa_json",
    "sample_examples",
    "dump_query_rows",
    "iter_query_rows",
]
