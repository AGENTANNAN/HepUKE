# Terms for the data in this directory

Three data products with three different provenances, so three different
answers.

## 1. `besiii_physicsqa/` — 1,000 question–answer pairs

**Authored for this work.** Question strings, answer values, scorer
assignments, hop chains and provenance records were produced by our pipeline
and reviewed by an author. The only third-party material is the short factual
values the answers quote (a dataset label, a cut value, an EvtGen model name)
— individual facts, not expression.

**Licence: CC BY 4.0.** Cite the paper when you use it.

## 2. `besiii_dsl_index/` — chunked index (no third-party text)

`chunks.jsonl` has 66,360 rows in two families:

| rows | source | content |
|---|---|---|
| 18,638 | `dsl` | full text of the DSL programs — see §3 |
| 47,722 | `md` | **text withheld**: `chunk_sha1`, `paper_id`, `kind`, section title and line range only, flagged `"content_withheld": "arxiv_license"` |

Why the text is withheld: the underlying 862 papers are on arXiv under
per-paper licences. A 36-paper random sample (seed 0, arXiv OAI `<license>`)
found **72% under `arxiv.org/licenses/nonexclusive-distrib/1.0`**, which grants
arXiv a distribution licence and grants us nothing; only ~22% were CC (BY 4.0,
BY 3.0, CC0, one BY-NC-SA 3.0). Auditing all 862 would leave roughly a fifth
redistributable, and nothing in this repository reads that text anyway —
`index_besiii_dsl_blocks_v1.py`, the only consumer, excludes `md` chunks by
design.

What this costs you: 1,657 of the 1,667 gold references point at `dsl` chunks
and are fully verifiable from this release. The remaining 10 (3 questions)
point at paper text you must fetch yourself. Every arXiv id is listed in
`besiii_dsl_index/arxiv_ids.txt`; note that chunk hashes were computed over a
particular PDF→markdown conversion, so refetching will not reproduce them
byte-for-byte.

`facts.jsonl` is filtered the same way: 48,865 `dsl`-derived facts kept,
170,424 `md`-derived facts (which quoted paper text in their `text` field)
removed. `annotations/` and `index.json` carry only section titles, table
captions, line ranges and extracted numeric values — structural metadata, no
prose.

**Licence: CC BY 4.0** for the structural metadata and the DSL-derived rows,
subject to §3 for the DSL text itself.

## 3. `corpus/` — 700 analyses

- `dsl/` — reference BOSS DSL programs, reviewed and confirmed by the authors
- `generated_dsl/` — model regenerations; these are what the knowledge base stores
- `descriptions/` — natural-language task descriptions derived from the reference programs

The DSL language, its manual and the generation prompt package originate in a
separate BOSS/DSL tooling project, not in this work; the programs in `dsl/`
are derived from it. The prompt package itself is **not** redistributed here —
set `DSL_SKILL_DIR` to your own copy if you want to regenerate.

**Licence: CC BY 4.0**, released with the permission of the upstream BOSS/DSL
tooling project.

> **Attribution — to be completed before de-anonymised release.** CC BY 4.0
> requires naming the upstream project. This package is anonymised for review,
> so the credit line is deliberately generic here. At camera-ready, replace
> this note with the upstream project's name and its preferred citation.

BESIII collaboration policy on analysis-specification material may apply
independently of the DSL project's permission. Check it separately.

## Not included anywhere in this release

Internal collaboration memos, the HyperNews mailing-list mirror, any
collaboration-internal corpus, API credentials, and service endpoints.
