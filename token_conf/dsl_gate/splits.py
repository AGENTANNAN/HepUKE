"""Deterministic train/dev/test split over the released DSL corpus.

The split is built from whatever is in ``data/corpus/descriptions/`` by
shuffling under a fixed seed, so it is reproducible but corpus-size dependent.
Default proportions (0.69 / 0.115 / rest) match the ratio used in the paper,
which partitioned a larger internal corpus 600/100/168. **Running this over
the released corpus does not reproduce the paper's partition** — define and
report your own split.

Also emits two perturbation subsets, ``pert20`` and ``pert50``, drawn from the
dev split under the same seed; the perturbation ablation runs on these.

Usage:
    python token_conf/dsl_gate/splits.py                 # write if missing
    python token_conf/dsl_gate/splits.py --force         # rewrite
    python token_conf/dsl_gate/splits.py --train 484 --dev 81 --test 135
    from token_conf.dsl_gate.splits import load_splits
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DESC_DIR = ROOT / "data" / "corpus" / "descriptions"
SPLIT_FILE = ROOT / "data" / "corpus" / "splits.json"

SEED = 42
# Proportions, applied to however many descriptions are present.
FRAC_TRAIN, FRAC_DEV = 0.69, 0.115
N_PERT20, N_PERT50 = 20, 50


def build_splits(n_train: int | None = None,
                 n_dev: int | None = None,
                 n_test: int | None = None,
                 seed: int = SEED) -> dict:
    stems = sorted(p.stem for p in DESC_DIR.glob("*.txt"))
    if not stems:
        raise FileNotFoundError(
            f"no descriptions found in {DESC_DIR} — the released corpus is "
            f"missing, or DSL_CORPUS_DIR points elsewhere")
    rng = random.Random(seed)
    rng.shuffle(stems)

    total = len(stems)
    if n_train is None:
        n_train = int(round(total * FRAC_TRAIN))
    if n_dev is None:
        n_dev = int(round(total * FRAC_DEV))
    if n_test is None:
        n_test = total - n_train - n_dev
    if min(n_train, n_dev, n_test) < 1 or n_train + n_dev + n_test > total:
        raise ValueError(
            f"invalid split {n_train}/{n_dev}/{n_test} for {total} stems")

    train = stems[:n_train]
    dev = stems[n_train:n_train + n_dev]
    test = stems[n_train + n_dev:n_train + n_dev + n_test]

    # Perturbation subsets: deterministic prefixes of a reshuffled dev split,
    # so pert20 is a subset of pert50.
    pert_pool = list(dev)
    random.Random(seed).shuffle(pert_pool)
    pert50 = pert_pool[:min(N_PERT50, len(pert_pool))]
    pert20 = pert50[:min(N_PERT20, len(pert50))]

    return {"seed": seed, "n_corpus": total,
            "train": train, "dev": dev, "test": test,
            "pert20": pert20, "pert50": pert50,
            "note": "Local split over the released corpus; NOT the partition "
                    "used in the paper."}


def load_splits() -> dict:
    if not SPLIT_FILE.exists():
        raise FileNotFoundError(
            f"{SPLIT_FILE} is missing — run:\n"
            f"    python token_conf/dsl_gate/splits.py")
    return json.loads(SPLIT_FILE.read_text())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true", help="overwrite an existing split")
    ap.add_argument("--train", type=int, default=None)
    ap.add_argument("--dev", type=int, default=None)
    ap.add_argument("--test", type=int, default=None)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    if SPLIT_FILE.exists() and not args.force:
        d = load_splits()
        print(f"splits already exist: train={len(d['train'])} dev={len(d['dev'])} "
              f"test={len(d['test'])} (use --force to rewrite)")
        return 0

    splits = build_splits(args.train, args.dev, args.test, args.seed)
    SPLIT_FILE.write_text(json.dumps(splits, indent=2))
    print(f"wrote {SPLIT_FILE}: corpus={splits['n_corpus']} "
          f"train={len(splits['train'])} dev={len(splits['dev'])} "
          f"test={len(splits['test'])} pert20={len(splits['pert20'])} "
          f"pert50={len(splits['pert50'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
