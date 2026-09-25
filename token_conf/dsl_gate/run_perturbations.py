"""Run the perturbation ablation: generate candidates under each variant.

Two families, one setting each (see variants.py):

    k10     retrieval depth doubled from the unperturbed top-5
    para    LLM paraphrase of the task description, retrieved at top-5

Each variant is one `generate_candidates.py` subprocess over the same split,
run sequentially, writing `<prefix>_<variant>.jsonl` into `runs_gate/`. The
`para` variant needs `build_paraphrases.py` to have run first.

Usage:
    python token_conf/dsl_gate/run_perturbations.py                  # both, pert20
    python token_conf/dsl_gate/run_perturbations.py --split pert50
    python token_conf/dsl_gate/run_perturbations.py --variants k10
    python token_conf/dsl_gate/run_perturbations.py --workers 8 --retries 1

Generation is resumable: `generate_candidates.py` skips targets already
present in the output, so re-running continues where the last run stopped.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

from variants import VARIANT_NAMES, desc_dir, tag, top_k   # noqa: E402

GATE = ROOT / "token_conf" / "dsl_gate" / "runs_gate"

# Each worker loads BGE-M3, which spawns one OMP thread pool per process. On a
# host with a low per-user process limit that exhausts it quickly, so pin the
# numeric libraries to a single thread and keep the worker count modest.
CHILD_ENV = {
    **os.environ,
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "TOKENIZERS_PARALLELISM": "false",
}


def run_one(variant: str, split: str, prefix: str, workers: int,
            temperature: float, log_dir: Path) -> int:
    k, dd = top_k(variant), desc_dir(variant)
    out_tag = tag(variant, prefix)
    cmd = [
        sys.executable, "-u",
        str(ROOT / "token_conf" / "dsl_gate" / "generate_candidates.py"),
        "--split", split,
        "--workers", str(workers),
        "--out-tag", out_tag,
        "--top-k", str(k),
        "--temperature", str(temperature),
    ]
    if dd is not None:
        if not dd.is_dir():
            print(f"[error] {variant} needs paraphrased descriptions in {dd}\n"
                  f"        run: python token_conf/dsl_gate/build_paraphrases.py "
                  f"--split {split}", flush=True)
            return 2
        cmd += ["--desc-dir", str(dd)]

    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{out_tag}.log"
    print(f"[run] {out_tag}  k={k}  desc={dd or '(unperturbed)'}  log={log}", flush=True)
    t0 = time.time()
    with log.open("w") as fp:
        rc = subprocess.call(cmd, stdout=fp, stderr=subprocess.STDOUT, env=CHILD_ENV)
    print(f"[done] {out_tag} rc={rc} wall={(time.time() - t0) / 60:.1f}min", flush=True)
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--split", default="pert20",
                    help="split name from data/corpus/splits.json (default: pert20)")
    ap.add_argument("--variants", nargs="+", default=VARIANT_NAMES,
                    choices=VARIANT_NAMES,
                    help="which perturbation settings to run (default: both)")
    ap.add_argument("--prefix", default=None,
                    help="output tag prefix (default: the split name)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--temperature", type=float, default=0.1)
    ap.add_argument("--retries", type=int, default=0,
                    help="retry a failed variant this many times before giving "
                         "up; generation returns non-zero if ANY target failed, "
                         "which is usually a transient API error")
    ap.add_argument("--keep-going", action="store_true",
                    help="continue with the next variant instead of aborting")
    ap.add_argument("--log-dir", type=Path, default=ROOT / "token_conf" / "dsl_gate" / "logs")
    args = ap.parse_args()

    prefix = args.prefix or args.split
    failed = []
    for variant in args.variants:
        rc = run_one(variant, args.split, prefix, args.workers,
                     args.temperature, args.log_dir)
        for attempt in range(args.retries):
            if rc == 0:
                break
            print(f"[retry {attempt + 1}/{args.retries}] {variant}", flush=True)
            rc = run_one(variant, args.split, prefix, args.workers,
                         args.temperature, args.log_dir)
        if rc != 0:
            failed.append(variant)
            if not args.keep_going:
                print(f"[abort] variant {variant} failed rc={rc} "
                      f"(use --keep-going to continue)", flush=True)
                return rc

    if failed:
        print(f"[done] {len(args.variants) - len(failed)}/{len(args.variants)} "
              f"variants ok; failed: {', '.join(failed)}", flush=True)
        return 1
    print(f"[all-done] {len(args.variants)} variants on split={args.split}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
