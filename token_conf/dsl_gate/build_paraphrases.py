"""Build the paraphrased task descriptions used by the `para` perturbation.

Each description in the chosen split is rewritten once by an LLM: different
sentence order and different connective phrases, with every physical fact
preserved verbatim — energy points, particle names, decay chains, dataset ids,
cut values, sample sizes, PID methods and kinematic-fit specs must survive
unchanged. Only the wording moves, so a score change downstream is attributable
to surface form rather than content.

Output goes to `token_conf/batch/desc_para/<target>.txt`, which
`run_perturbations.py --variants para` then passes to generation as
`--desc-dir`.

Usage:
    python token_conf/dsl_gate/build_paraphrases.py                 # pert20
    python token_conf/dsl_gate/build_paraphrases.py --split pert50
    python token_conf/dsl_gate/build_paraphrases.py --split dev --force

Resumable: a target whose paraphrase already exists and is non-empty is
skipped unless `--force` is given.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "token_conf"))
sys.path.insert(0, str(ROOT / "token_conf" / "dsl_gate"))

import pipeline_common as pc                            # noqa: E402
from splits import load_splits                          # noqa: E402
from variants import PARA_DESC_DIR                      # noqa: E402

DESC_ROOT = ROOT / "data" / "corpus" / "descriptions"

PARA_MODEL = pc.model_id("deepseek-ai/deepseek-v4.1-flash")
PARA_TEMP = 0.5
PARA_MAX_TOK = 4000

PARA_SYS = ("You are a paraphrasing assistant for BESIII high-energy physics "
            "analysis descriptions. Rewrite the given description in different "
            "wording but preserve every physical fact verbatim: energy points, "
            "particle names, decay chains, dataset ids, cut values, sample sizes, "
            "PID methods, and kinematic-fit specs must all appear unchanged. Do "
            "not add commentary; output only the rewritten description.")

PARA_HINT = "Rewrite with a different sentence order and different connective phrases."


def paraphrase(desc: str, base: str, key: str) -> str:
    user = f"{PARA_HINT}\n\nOriginal description:\n{desc}\n\nRewritten description:"
    r = pc.chat(system=PARA_SYS, user=user, model=PARA_MODEL,
                base_url=base, api_key=key,
                max_tokens=PARA_MAX_TOK, temperature=PARA_TEMP,
                logprobs=False, stream=False, timeout=300.0,
                extra_payload={"reasoning_effort": "none"})
    return (r.get("answer") or "").strip()


def build(targets: list[str], force: bool = False) -> Path:
    base, key = pc.load_cfg()
    PARA_DESC_DIR.mkdir(parents=True, exist_ok=True)
    made = skipped = failed = 0
    for q in targets:
        out = PARA_DESC_DIR / f"{q}.txt"
        if out.exists() and out.stat().st_size > 0 and not force:
            skipped += 1
            continue
        src = DESC_ROOT / f"{q}.txt"
        if not src.is_file():
            print(f"  [skip] no description for {q}", flush=True)
            failed += 1
            continue
        t0 = time.time()
        text = paraphrase(src.read_text().strip(), base, key)
        if not text:
            print(f"  [fail] empty paraphrase for {q}", flush=True)
            failed += 1
            continue
        out.write_text(text, encoding="utf-8")
        made += 1
        print(f"  {q}  {time.time() - t0:.1f}s  {len(text)} chars", flush=True)
    print(f"[done] {PARA_DESC_DIR}: {made} written, {skipped} already present, "
          f"{failed} failed", flush=True)
    return PARA_DESC_DIR


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--split", default="pert20",
                    help="split name from data/corpus/splits.json (default: pert20)")
    ap.add_argument("--force", action="store_true",
                    help="rewrite paraphrases that already exist")
    args = ap.parse_args()

    splits = load_splits()
    if args.split not in splits:
        raise SystemExit(f"no split {args.split!r} in splits.json "
                         f"(have: {', '.join(k for k, v in splits.items() if isinstance(v, list))})")
    targets = splits[args.split]
    print(f"[paraphrase] split={args.split} targets={len(targets)} model={PARA_MODEL}",
          flush=True)
    build(targets, force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
