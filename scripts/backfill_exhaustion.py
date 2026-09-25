"""Backfill the `exhaustion` φ feature onto legacy predictor-training rows.

Legacy rows (collected before commit adding `exhaustion` to
`extract_features`) already carry `doc_ids` (new IDs this round) but not the
new feature. Compute it in-place as `1 - min(1, len(doc_ids) / expected)`
so retrained predictors see the same feature distribution as newly-collected
rows. Writes to a new file; the original is left untouched.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="input", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--expected-docs-per-round", type=int, default=5)
    args = p.parse_args()

    denom = max(1, int(args.expected_docs_per_round))
    n_read = 0
    n_backfilled = 0
    n_skipped = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.input.open("r", encoding="utf-8") as fin, args.out.open("w", encoding="utf-8") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            n_read += 1
            feats = row.get("features")
            if not isinstance(feats, dict):
                n_skipped += 1
                fout.write(json.dumps(row, ensure_ascii=False) + "\n")
                continue
            if "exhaustion" in feats:
                n_skipped += 1
                fout.write(json.dumps(row, ensure_ascii=False) + "\n")
                continue
            doc_ids = row.get("doc_ids") or []
            n_new = len({d for d in doc_ids if d})
            exhaustion = 1.0 - min(1.0, n_new / denom)
            feats["exhaustion"] = float(exhaustion)
            n_backfilled += 1
            fout.write(json.dumps(row, ensure_ascii=False) + "\n")

    print(f"read={n_read} backfilled={n_backfilled} skipped_or_existing={n_skipped}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
