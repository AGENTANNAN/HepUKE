"""Train the look-ahead KL-gain predictor from collected HotpotQA trajectories.

Reads ``storage/predictor_train/*.jsonl`` (produced by
``scripts/collect_hotpotqa_trajectories.py``) and fits three artifacts that
together form the ``predictor_ckpt`` consumed by the runtime controller
(``hepuke/agent/confidence/predictor.py``):

  * **KL head** — ``MLPRegressor`` trained on ``(φ, g_next)`` pairs, i.e.
    rows where the next round exists. Predicts the marginal KL gain of
    running one more retrieval.
  * **Isotonic calibrator** — monotone map applied to the MLP's raw output
    so downstream fusion sees a comparable [0,1] score.
  * **Early-stop head** — ``LogisticRegression`` trained on **all** rows
    with a label ``y = 1 iff this row is the last round of an
    already-terminated query and ``t == 1```` (i.e. the one-shot success
    signal). Uses the same 9-dim φ so no separate feature engineering.

Rationale for the auxiliary head: 60–70 % of HotpotQA queries terminate at
t=1 in our collector, and those rows have ``g_next = None`` so they cannot
enter the KL regression. Discarding them wastes majority of the samples;
the classifier reuses them as evidence for "what does an easy query's φ
look like?" without polluting the KL head's label distribution.

Usage::

    python scripts/train_predictor.py \\
        --data storage/predictor_train/hotpotqa_train.gpt4omini.jsonl \\
        --out  storage/predictor.pkl \\
        --seed 0

Diagnostic mode (no checkpoint written)::

    python scripts/train_predictor.py --data ... --dry-run
"""
from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any


FEATURE_KEYS = (
    "t",
    "entropy",
    "top1_prob",
    "posterior_var",
    "n_unique_docs",
    "novelty",
    "gain_mean",
    "gain_var",
    "gain_slope",
    "exhaustion",
)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True, type=Path,
                   help="JSONL of per-round rows from collect_hotpotqa_trajectories.py.")
    p.add_argument("--out", type=Path, default=Path("storage/predictor.pkl"),
                   help="Where to pickle the (mlp, iso, aux, meta) bundle.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--val-frac", type=float, default=0.15,
                   help="Held-out fraction for MLP early-stopping / calibrator.")
    p.add_argument("--hidden", type=int, nargs="+", default=[32, 16],
                   help="MLPRegressor hidden layer sizes.")
    p.add_argument("--max-iter", type=int, default=500)
    p.add_argument("--filter-onehot", action="store_true", default=True,
                   help="Drop rows whose posterior is essentially one-hot "
                        "(entropy<--onehot-entropy AND top1>--onehot-top1) "
                        "from the KL head training set. Default: on. These "
                        "rows have degenerate KL/JSD labels that inject "
                        "orders-of-magnitude noise into g_next.")
    p.add_argument("--no-filter-onehot", dest="filter_onehot", action="store_false")
    p.add_argument("--onehot-entropy", type=float, default=0.01)
    p.add_argument("--onehot-top1", type=float, default=0.99)
    p.add_argument("--dry-run", action="store_true",
                   help="Fit but do not write the pickle — for smoke checks.")
    return p.parse_args()


def _load_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    return rows


def _is_onehot(row: dict, entropy_thresh: float, top1_thresh: float) -> bool:
    """Row's posterior is essentially a delta at the top-1 token.

    Detected via the stored feature summary (avoids re-parsing the raw
    posterior). Falls back to False when the required features are missing.
    """
    feats = row.get("features") or {}
    try:
        ent = float(feats.get("entropy", 0.0))
        top1 = float(feats.get("top1_prob", 0.0))
    except (TypeError, ValueError):
        return False
    return ent < entropy_thresh and top1 > top1_thresh


def _extract_kl_pool(
    rows: list[dict],
    *,
    filter_onehot: bool = True,
    entropy_thresh: float = 0.01,
    top1_thresh: float = 0.99,
) -> tuple[list[list[float]], list[float], int]:
    """Rows where g_next is defined → training set for the KL head.

    Optionally drops rows whose posterior is essentially one-hot; those have
    degenerate JSD labels that add noise without signal (see
    ``docs/predictor_data_hygiene.md``). Returns the third element as the
    count of dropped one-hot rows so the caller can log the effective size.
    """
    X: list[list[float]] = []
    y: list[float] = []
    dropped_onehot = 0
    for r in rows:
        g = r.get("g_next")
        if g is None:
            continue
        if filter_onehot and _is_onehot(r, entropy_thresh, top1_thresh):
            dropped_onehot += 1
            continue
        feats = r.get("features") or {}
        try:
            vec = [float(feats[k]) for k in FEATURE_KEYS]
        except (KeyError, TypeError, ValueError):
            continue
        X.append(vec)
        y.append(float(g))
    return X, y, dropped_onehot


def _extract_earlystop_pool(rows: list[dict]) -> tuple[list[list[float]], list[int]]:
    """Auxiliary classifier: y=1 iff this row is a one-shot terminal round.

    We label at the *query* level: mark all rows whose query had ``max(t)==1``
    (agent stopped after the first retrieval) as positive; all others as
    negative. Every row participates (bridge to the KL head, which drops
    one-shot rows entirely).
    """
    max_t_by_qid: dict[str, int] = defaultdict(int)
    for r in rows:
        qid = r.get("query_id")
        t = int(r.get("t", 0) or 0)
        if qid is not None and t > max_t_by_qid[qid]:
            max_t_by_qid[qid] = t

    X: list[list[float]] = []
    y: list[int] = []
    for r in rows:
        feats = r.get("features") or {}
        try:
            vec = [float(feats[k]) for k in FEATURE_KEYS]
        except (KeyError, TypeError, ValueError):
            continue
        qid = r.get("query_id")
        label = 1 if max_t_by_qid.get(qid, 0) == 1 else 0
        X.append(vec)
        y.append(label)
    return X, y


def _split(X: list, y: list, val_frac: float, seed: int) -> tuple[list, list, list, list]:
    import random
    rng = random.Random(seed)
    idx = list(range(len(X)))
    rng.shuffle(idx)
    n_val = max(1, int(len(idx) * val_frac))
    val_idx = set(idx[:n_val])
    Xt, yt, Xv, yv = [], [], [], []
    for i, (x, yi) in enumerate(zip(X, y)):
        if i in val_idx:
            Xv.append(x); yv.append(yi)
        else:
            Xt.append(x); yt.append(yi)
    return Xt, yt, Xv, yv


def _fit_kl_head(
    X_train: list[list[float]],
    y_train: list[float],
    X_val: list[list[float]],
    y_val: list[float],
    *,
    hidden: tuple[int, ...],
    max_iter: int,
    seed: int,
    log: logging.Logger,
):
    from sklearn.neural_network import MLPRegressor
    from sklearn.isotonic import IsotonicRegression

    if not X_train:
        log.warning("KL pool empty — skipping MLP fit; returning identity stub.")
        return None, None

    mlp = MLPRegressor(
        hidden_layer_sizes=tuple(hidden),
        activation="relu",
        solver="adam",
        max_iter=max_iter,
        early_stopping=False,  # tiny dataset — full-batch is fine
        random_state=seed,
    )
    t0 = time.time()
    mlp.fit(X_train, y_train)
    train_loss = float(mlp.loss_)
    log.info("MLP trained in %.1fs, final train loss=%.4f", time.time() - t0, train_loss)

    # Calibrate on the held-out set (raw MLP output → isotonic to true g_next).
    iso: IsotonicRegression | None = None
    if X_val:
        raw_val = mlp.predict(X_val)
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=None)
        iso.fit(raw_val, y_val)
        cal_val = iso.transform(raw_val)
        val_mse = float(sum((c - y) ** 2 for c, y in zip(cal_val, y_val)) / len(y_val))
        log.info("Isotonic calibration val MSE=%.4f", val_mse)
    else:
        log.warning("No validation data — isotonic calibrator skipped.")
    return mlp, iso


def _fit_earlystop_head(
    X: list[list[float]],
    y: list[int],
    *,
    seed: int,
    log: logging.Logger,
):
    from sklearn.linear_model import LogisticRegression

    if not X or len(set(y)) < 2:
        log.warning("Early-stop pool degenerate (n=%d, classes=%s) — skipping.",
                    len(X), set(y))
        return None

    clf = LogisticRegression(
        random_state=seed,
        max_iter=1000,
        class_weight="balanced",  # keep the minority class visible
    )
    t0 = time.time()
    clf.fit(X, y)
    acc = float(clf.score(X, y))
    log.info(
        "Early-stop head trained in %.1fs, in-sample acc=%.3f (n=%d, pos=%d)",
        time.time() - t0, acc, len(X), sum(y),
    )
    return clf


def main() -> int:
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    log = logging.getLogger("train_predictor")

    rows = _load_rows(args.data)
    log.info("loaded %d rows from %s", len(rows), args.data)
    if not rows:
        log.error("input file empty — nothing to train")
        return 2

    # -------- KL head (MLP + isotonic)
    X_kl, y_kl, dropped = _extract_kl_pool(
        rows,
        filter_onehot=args.filter_onehot,
        entropy_thresh=args.onehot_entropy,
        top1_thresh=args.onehot_top1,
    )
    log.info(
        "KL pool: %d labeled rows (of %d total, %d one-shot dropped, "
        "%d one-hot dropped [filter=%s])",
        len(X_kl), len(rows),
        len(rows) - len(X_kl) - dropped,  # one-shot rows (g_next=None)
        dropped, args.filter_onehot,
    )
    Xt, yt, Xv, yv = _split(X_kl, y_kl, args.val_frac, args.seed)
    log.info("KL split: train=%d val=%d", len(Xt), len(Xv))
    mlp, iso = _fit_kl_head(
        Xt, yt, Xv, yv,
        hidden=tuple(args.hidden),
        max_iter=args.max_iter,
        seed=args.seed,
        log=log,
    )

    # -------- Early-stop auxiliary head
    X_es, y_es = _extract_earlystop_pool(rows)
    log.info("Early-stop pool: %d rows, %d positive (one-shot)",
             len(X_es), sum(y_es))
    aux = _fit_earlystop_head(X_es, y_es, seed=args.seed, log=log)

    # -------- Bundle
    bundle: dict[str, Any] = {
        "feature_keys": list(FEATURE_KEYS),
        "kl_mlp": mlp,
        "kl_isotonic": iso,
        "early_stop_clf": aux,
        "meta": {
            "n_rows_total": len(rows),
            "n_kl_labeled": len(X_kl),
            "n_onehot_dropped": dropped,
            "n_earlystop_positive": int(sum(y_es)),
            "seed": args.seed,
            "hidden": list(args.hidden),
            "val_frac": args.val_frac,
            "filter_onehot": bool(args.filter_onehot),
            "onehot_entropy": float(args.onehot_entropy),
            "onehot_top1": float(args.onehot_top1),
            "label_type": "jsd",  # switched from KL to JSD in T8.1
            "data_path": str(args.data),
            "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
    }

    if args.dry_run:
        log.info("dry-run: skipping pickle write")
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("wb") as f:
        pickle.dump(bundle, f)
    log.info("wrote %s (%.1f KB)", args.out, args.out.stat().st_size / 1024)
    return 0


if __name__ == "__main__":
    sys.exit(main())
