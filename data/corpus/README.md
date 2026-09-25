# BESIII DSL corpus (700 analyses)

Three parallel directories, keyed by arXiv stem (e.g. `1001.5328v1`), plus
`MANIFEST.csv` listing every stem with the SHA-256 of all three files.

| directory | n | what it is |
|---|---|---|
| `dsl/` | 700 | reference BOSS DSL program for the analysis — reviewed and confirmed by the authors |
| `descriptions/` | 700 | natural-language task description derived from the reference program; this is the retrieval query and the generator's input |
| `generated_dsl/` | 700 | model regeneration of the same analysis from its description |

`generated_dsl/` is not a duplicate of `dsl/`. The knowledge base stores and
retrieves the **generated** programs together with their token-level entropy;
`dsl/` is the reference the evaluator scores against. `build_index.py`
indexes `descriptions/` and attaches the matching `generated_dsl/` program to
each row.

## Which 700

The full internal corpus holds 868 analyses. This release contains:

- the **594** programs that actually participated in the reported runs — every
  evaluation target plus every program ever returned as a retrieved precedent;
- **106** further programs, taken in the frozen seed-42 split order
  (`train` then `dev` then `test`, first 106 not already among the 594), so
  that the selection rule is reproducible rather than ad hoc.

## Two consequences for reproduction

1. **Retrieval ranking will not match the paper exactly.** The reported index
   was built over all 868 descriptions; rebuilding over these 700 changes
   which candidates enter the top-5 for some queries. Scoring, entropy
   computation and the gating logic reproduce exactly; the candidate pool does
   not.
2. **The evaluation split is not shipped.** `splits.py` regenerates a
   600/100/168 partition by shuffling whatever is in `descriptions/` under
   seed 42. Run over these 700 files it produces a *different* partition than
   the one used in the paper. Define your own split and report it.

## Regenerating instead of using the snapshot

`token_conf/build_descriptions.py` and `token_conf/run_skill_batch.py`
regenerate `descriptions/` and `generated_dsl/` from reference programs. They
need the DSL generation prompt package, which is a third-party artifact and is
**not** redistributed here; point `DSL_SKILL_DIR` at your own copy, and
`DSL_CORPUS_DIR` at your own reference programs.

Licence: CC BY 4.0, released with the permission of the upstream BOSS/DSL
tooling project. See `../LICENSE-DATA.md` — the attribution line is generic
in this anonymised package and must be completed at de-anonymisation.
