# BESIII-PhysicsQA

1,000 question–answer pairs over BESIII analysis programs and their source
publications. Each item names the scorer it is to be graded with and points at
the gold evidence chunks that support the answer.

`manifest.json` carries the frozen counts, the SHA-256 of `qa.jsonl`, the list
of corrections applied, and the known issues repeated below.

## Schema (`qa.jsonl`, one JSON object per line)

| field | meaning |
|---|---|
| `qid` | stable item id, e.g. `besiii-A-0055` |
| `question` | the question string |
| `answer` | gold answer, graded under `scorer` |
| `scorer` | `exact` (352) / `numeric` (517) / `set_f1` (128) / `bool_plus_numeric` (3) |
| `n_hops` | number of gold facts that must be combined (1: 871, ≥2: 129) |
| `hop_chain` | ordered list of gold chunk references, one group per hop |
| `gold_sha1` | SHA-1 of each gold chunk; resolves against `../besiii_dsl_index/chunks.jsonl` (`chunk_sha1`) |
| `task_class` | A (869) / B (75) / C (27) / D (3) / E (21) / F (5) |
| `scenario` | one-line description of the analysis task the question comes from |
| `cluster_key`, `paper_ids` | source publication(s) |
| `provenance` | fact ids, locators, index version, prompt hash, build timestamp |

## Scoring

Reference implementation: `scripts/besiii_set_scorer.py` in this repository.

- `exact` — normalised string equality.
- `numeric` — the gold's numeric tokens must appear in the prediction. It is
  deliberately loose because predictions embed numbers in prose.
- `set_f1` — set overlap over normalised items.

Grade by `scorer`; do not pool a single metric across all four families
without saying so.

## How it was built

Every factual decision — which chunks form the hop chain, which value is the
answer, which scorer applies — comes from a deterministic query over a
pre-built fact index. A language model was used only to phrase an
already-fixed chain into a natural question, and the resulting items were
reviewed by an author.

## Corrections applied to this snapshot

Seven items carried the gold value as a Ruby integer literal with digit
separators (`300_000`). The numeric scorer splits that into `300` and `000`,
so a correct prediction scored zero. The separators were removed; nothing else
about those items changed. The exact list is in `manifest.json`
(`corrections_applied`).

## Known issues — read before reporting numbers

1. **127 numeric-typed dataset labels.** Items whose answer is a dataset
   identifier such as `data_3773` carry `scorer="numeric"`. The numeric scorer
   extracts the digit run, so they still discriminate correctly against other
   datasets, and they are left as-is here. A stricter exact-match reading is
   defensible and would make these items harder; if you re-type them, say so.
2. **Review provenance.** Items were reviewed by an author. This release
   carries no per-item review record, no independent multi-rater annotation
   and no inter-annotator agreement statistic. Do not describe the set as
   multiply-annotated or as carrying expert correctness labels for anything
   beyond the answer fields.
3. **Class imbalance.** 869 of 1,000 items are class A and 871 are
   single-hop. The multi-hop portion is small; report it separately rather
   than averaging it away.
4. **Gold resolution needs the index.** `gold_sha1` is only meaningful against
   the accompanying `besiii_dsl_index`. 1,657 of the 1,667 gold references
   point at DSL chunks and resolve completely. The remaining 10 (3 questions)
   point at source-paper chunks whose text is withheld for licensing reasons —
   the rows and their hashes are present, the text is not. See
   `../LICENSE-DATA.md`.
