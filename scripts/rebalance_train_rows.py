"""Rebalance the predictor training set by up-sampling multi-round trajectories.

MLPRegressor.fit() does not accept sample_weight, so we express the weighting
as row duplication: rows belonging to a query whose max_t >= 2 are copied
`weight` times. Passing --weights 1 2 4 4 8 8 16 16 means rows in rounds=1
queries appear once, rounds=2 twice, rounds=3 four times, etc.

Writes a new JSONL. Feed that into scripts/train_predictor.py as usual.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="input", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--weights", type=int, nargs="+",
                   default=[1, 2, 4, 4, 8, 8, 16, 16],
                   help="One weight per bucket: index i -> rounds i+1 (last bucket catches all >= len).")
    args = p.parse_args()

    weights = list(args.weights)

    # First pass: per-query max_t.
    max_t: dict[str, int] = {}
    rows: list[dict] = []
    with args.input.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            rows.append(r)
            qid = r.get("query_id")
            t = int(r.get("t", 0) or 0)
            if qid is not None and t > max_t.get(qid, 0):
                max_t[qid] = t

    def _weight_for(mt: int) -> int:
        if mt <= 0:
            return 1
        idx = min(mt - 1, len(weights) - 1)
        return max(1, int(weights[idx]))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    per_bucket_qids = defaultdict(set)
    per_bucket_rows = Counter()
    with args.out.open("w", encoding="utf-8") as fout:
        for r in rows:
            qid = r.get("query_id")
            mt = max_t.get(qid, 0)
            w = _weight_for(mt)
            per_bucket_qids[mt].add(qid)
            per_bucket_rows[mt] += w
            for _ in range(w):
                fout.write(json.dumps(r, ensure_ascii=False) + "\n")
                n_written += 1

    print(f"input rows: {len(rows)}")
    print(f"output rows: {n_written}")
    print("per-bucket stats (max_t -> unique_qids, weighted_rows):")
    for mt in sorted(per_bucket_qids):
        print(f"  max_t={mt}: qids={len(per_bucket_qids[mt])}  rows_after_weight={per_bucket_rows[mt]}  weight={_weight_for(mt)}x")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
