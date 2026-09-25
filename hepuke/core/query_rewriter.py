"""LLM-driven multi-candidate query rewriter.

The old `Preprocessor` (CN→EN dictionary + regex) is a static fallback.
This module is the ROOT solution to the "语料库全英，用户问中文" problem:
a small LLM turns any-language physics question into a set of English
retrieval candidates, each with its own dense (full sentence) and sparse
(keyword bag with technical expansions) text.

Downstream (`HybridRetriever.search`) fans out over the candidates and
RRF-merges the ranked lists. That covers the "Zc(3900) → charmonium /
charmonium-like / exotic hadron / tetraquark" recall problem WITHOUT a
hand-maintained词表 — the LLM produces the expansions, so novel entities
work on day one.

Design notes
------------
* Deterministic backward compat: this class is a drop-in `Callable[[str],
  tuple[str, str]]` — calling it directly yields ONE pair (the top
  candidate). Multi-candidate retrieval calls `.rewrite_multi(...)`
  explicitly, so old paths keep working unchanged.
* LRU cache keyed on the raw question. Same question in a session or
  across sessions never re-pays the LLM cost.
* Hard three-tier fallback: LLM fine → LLM broken/timeout → old
  `Preprocessor`. The rewriter must NEVER be a source of retrieval
  failure — worst case is that we're no better than the old path.
"""
from __future__ import annotations

import functools
import json
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Optional

from hepuke.core.query_preprocessor import Preprocessor, normalize_unicode_physics


__all__ = ["Candidate", "LLMQueryRewriter"]


_JSON_RE = re.compile(r"\{.*\}", re.DOTALL)


@dataclass
class Candidate:
    """One retrieval hypothesis. `dense_text` feeds the semantic vector
    channel; `sparse_text` feeds the lexical (BM25-like) channel."""
    dense_text: str
    sparse_text: str
    label: str = ""  # human-readable tag e.g. "translated", "expansion:charmonium"
    expansions: list[str] = field(default_factory=list)

    def as_pair(self) -> tuple[str, str]:
        return self.dense_text, self.sparse_text


REWRITER_SYSTEM_PROMPT = """\
You rewrite a physics research question (any language) into multiple
English retrieval queries for a hybrid dense+sparse search over an
ENGLISH corpus.

Output STRICT JSON, no prose:
{
  "en_question": "<question translated + normalised to English>",
  "candidates": [
    {
      "label": "<short tag>",
      "query": "<a full English query sentence>",
      "expansions": ["<related English technical term>", ...]
    },
    ...
  ]
}

Rules:
- Produce 3-5 candidates covering complementary angles:
    (1) a faithful translation of the question,
    (2) a broader concept query (parent field / physics category),
    (3) one or two entity-anchored queries expanding named states
        into the term families a paper would actually use.
- Preserve every named entity verbatim in `query`: `Zc(3900)`,
  `BAM-00097`, `J/psi`, `pi+pi-J/psi`, `arXiv:2401.01234`, etc.
- `expansions`: 3-8 English technical terms that MIGHT appear in an
  English paper on the topic but that the user did NOT say. Examples:
    Zc(3900) → charmonium-like, exotic hadron, tetraquark, hidden charm,
              charged charmonium state, e+e- annihilation
    X(3872)  → exotic charmonium, DDbar molecule, charmonium-like state
    Higgs    → electroweak symmetry breaking, Yukawa coupling, gauge boson
    lattice  → lattice QCD, quenched, Wilson action
- Do NOT invent doc_ids or citations. Do NOT answer the question.
- If the question is already fully English, still emit multiple angles.
"""

_FEW_SHOT_USER = "Zc(3900) 的观测工作有哪些"
_FEW_SHOT_ASSISTANT = json.dumps({
    "en_question": "Observation of Zc(3900)",
    "candidates": [
        {
            "label": "faithful",
            "query": "Observation of Zc(3900)",
            "expansions": ["e+e- annihilation", "BESIII", "pi+pi-J/psi"],
        },
        {
            "label": "concept",
            "query": "charged charmonium-like exotic hadron states",
            "expansions": ["charmonium-like", "exotic hadron", "tetraquark", "hidden charm"],
        },
        {
            "label": "entity-expansion",
            "query": "Zc(3900) tetraquark interpretation and decay channels",
            "expansions": ["tetraquark", "DDbar* molecule", "charged charmonium state"],
        },
    ],
}, ensure_ascii=False)


# --------------------------------------------------------------------------- #
# Code-register variant — used when the corpus is BOSS Ruby DSL source, not
# prose. Both REWRITER_SYSTEM_PROMPT ("terms a paper would use") and
# GENERIC_REWRITER_SYSTEM_PROMPT ("encyclopedia-style corpus") expand into
# publication vocabulary, which is out of distribution against Ruby code: with
# prompt_style="generic" the rewriter was measured emitting queries like
# "BESIII dataset 2.396 GeV luminosity sample" against a corpus whose text is
# `DatasetManager.real_data.find("713_Rscan_2396")`.
# Purely additive — neither existing style is altered.
# --------------------------------------------------------------------------- #
DSL_REWRITER_SYSTEM_PROMPT = """\
You rewrite a question about a BESIII analysis into English retrieval queries
for a hybrid dense+sparse search over a corpus of **BOSS Ruby DSL source
code**, not papers. Every indexed document is one logical block of code.
Write queries in the register of the code itself.

Each document carries a `kind` drawn from exactly this list:
  dataset       DatasetManager.real_data / .inclusive_mc .find("...")
  decay_card    EvtGen heredoc: Decay / Enddecay / PHSP / ETA_DALITZ / VSS
  exclusive_mc  create_exclusive_mc, config.sample_name, config.events
  algorithm     Algorithm.new, set_header, set_constant, set_alias
  selection     Selection.new -> select_track / select_photon / pid /
                kinematic_fit, with cos_theta / Vz / Vr / nChrp / nChrn /
                chi2_cut / mass_window / nC
  tag_analysis  TagAnalysis, D-tag / Lambda_c-tag, mBC / deltaE, other side
  note          .note(...) free-text constraints
  execute       with_decay_card / apply / execute_on / root_files /
                save_to_config
  header        leading comments only — rarely holds a setting

Output STRICT JSON, no prose:
{
  "en_question": "<question normalised to English>",
  "candidates": [
    {"label": "<one kind from the list above>",
     "query": "<retrieval query in code register>",
     "expansions": ["<DSL identifier or method name>", ...]},
    ...
  ]
}

Rules:
- Produce 3-4 candidates. Put in `label` the block `kind` you believe holds
  the answer, chosen only from the list above. The caller filters retrieval by
  that kind, so a wrong label costs recall.
- Never label a candidate `header` unless the question is about a comment.
- `query` must read like the code, not like a paper. Prefer DSL method and
  parameter names over physics prose: write
  `select_photon energy threshold barrel endcap`, not
  "minimum photon energy required in the event selection".
- Preserve every named entity verbatim: `psi(3686)`, `J/psi`, `chi_c0`,
  `pi+pi-J/psi`, `Zc(3900)`, and energy points such as `3.686`, `2.396`.
- `expansions`: 3-8 DSL identifiers, method names, or parameter keys that
  would literally appear in the target block. NOT paper vocabulary — never
  emit terms like "charmonium-like", "tetraquark", "cross section
  measurement", "luminosity".
- Dataset and decay-model names are code identifiers with underscore
  spellings (`psip_data`, `data_3686`, `713_Rscan_2396`, `PHSP`,
  `ETA_DALITZ`). Put plausible spelling variants in `expansions`, because the
  same physical sample is named differently across programs.
- Do NOT answer the question. Do NOT explain.
"""

_DSL_FEW_SHOT_USER = (
    "For the chi_c0 -> pi0 pi0 pi0 pi0 analysis, what is the minimum number "
    "of pi0 required in the event selection?"
)
_DSL_FEW_SHOT_ASSISTANT = json.dumps({
    "en_question": (
        "For the chi_c0 -> pi0 pi0 pi0 pi0 analysis, what is the minimum "
        "number of pi0 required in the event selection?"
    ),
    "candidates": [
        {
            "label": "selection",
            "query": ("Selection.new select_photon pi0 count nPi0 minimum "
                      "chi_c0 four pi0 final state"),
            "expansions": ["select_photon", "nPi0", "select_track",
                           "cos_theta", "Vz", "kinematic_fit", "mass_window"],
        },
        {
            "label": "selection",
            "query": "kinematic_fit nC chi2_cut mass_window pi0 gamma gamma",
            "expansions": ["kinematic_fit", "nC", "chi2_cut", "mass_window",
                           "pi0"],
        },
        {
            "label": "exclusive_mc",
            "query": ("create_exclusive_mc config.sample_name chi_c0 4pi0 "
                      "config.events"),
            "expansions": ["create_exclusive_mc", "config.sample_name",
                           "config.events", "config.related_dataset"],
        },
        {
            "label": "decay_card",
            "query": "Decay chi_c0 pi0 pi0 pi0 pi0 PHSP Enddecay",
            "expansions": ["Decay", "Enddecay", "PHSP", "chi_c0", "pi0"],
        },
    ],
}, ensure_ascii=False)


# --------------------------------------------------------------------------- #
# Domain-agnostic variant — used for open-domain multi-hop QA (HotpotQA etc.)
# where the physics-specific system prompt causes small LLMs to "helpfully"
# answer the question in prose instead of emitting JSON candidates.
# --------------------------------------------------------------------------- #
GENERIC_REWRITER_SYSTEM_PROMPT = """\
You rewrite a user question into multiple English retrieval queries for
a hybrid dense+sparse search over an English encyclopedia-style corpus.

DO NOT answer the question. DO NOT explain your reasoning.
Output STRICT JSON matching this exact shape, and NOTHING else:

{
  "en_question": "<question, translated to English if needed>",
  "candidates": [
    {"label": "<short tag>", "query": "<retrieval query sentence>",
     "expansions": ["<related term>", ...]},
    ...
  ]
}

Rules:
- Produce 3 candidates covering complementary angles:
    (1) a faithful phrasing of the question,
    (2) a phrasing that isolates one of the two entities in the question
        so the retriever finds THAT entity's page,
    (3) a phrasing that isolates the OTHER entity.
- Preserve every named entity verbatim (people, films, places, dates).
- `expansions`: 3-6 related English terms that MIGHT appear on the
  target page but are NOT in the question itself (occupation, era,
  medium, genre, etc.).
- If you are tempted to write prose, STOP and emit the JSON instead.
"""

_GENERIC_FEW_SHOT_USER = (
    "Which film starring both Humphrey Bogart and Ingrid Bergman won "
    "the Academy Award for Best Picture?"
)
_GENERIC_FEW_SHOT_ASSISTANT = json.dumps({
    "en_question": (
        "Which film starring both Humphrey Bogart and Ingrid Bergman "
        "won the Academy Award for Best Picture?"
    ),
    "candidates": [
        {
            "label": "faithful",
            "query": (
                "Film starring Humphrey Bogart and Ingrid Bergman that "
                "won the Academy Award for Best Picture"
            ),
            "expansions": ["Casablanca", "1943 film", "Warner Bros"],
        },
        {
            "label": "entity-A",
            "query": "Humphrey Bogart filmography Academy Award",
            "expansions": ["American actor", "1940s Hollywood", "film noir"],
        },
        {
            "label": "entity-B",
            "query": "Ingrid Bergman filmography Academy Award",
            "expansions": ["Swedish actress", "1940s Hollywood", "Oscar-winning film"],
        },
    ],
}, ensure_ascii=False)


class LLMQueryRewriter:
    """Callable rewriter with `.rewrite_multi(...)` for candidate lists.

    Parameters
    ----------
    llm : object with `.chat_with_tools(messages, tools, tool_choice)`
        (the same LLM protocol the RagAgent uses). We only ever call it
        with `tools=[]` / `tool_choice="none"` — this is a plain
        completion request, not tool use.
    fallback : Preprocessor used when the LLM fails, times out, or emits
        broken JSON. Defaults to a plain `Preprocessor()`.
    cache_size : LRU cap on `.rewrite_multi()` calls keyed by the raw
        question. 0 disables the cache (tests use this).
    max_candidates : upper bound on how many rewrites we keep from the
        LLM's output — protects against a runaway LLM emitting dozens.
    """

    def __init__(
        self,
        *,
        llm: Any,
        fallback: Optional[Preprocessor] = None,
        cache_size: int = 512,
        max_candidates: int = 5,
        keep_original: bool = False,
        prompt_style: str = "physics",
    ) -> None:
        self.llm = llm
        self.fallback = fallback or Preprocessor()
        self.max_candidates = int(max_candidates)
        self.keep_original = bool(keep_original)
        # `prompt_style` picks the system prompt: "physics" (default, keeps
        # BESIII paper-corpus behaviour untouched), "generic" for open-domain
        # multi-hop QA like HotpotQA, or "dsl" when the corpus is BOSS Ruby
        # source rather than prose. See 5ade9c9c probe: the physics prompt made
        # flash answer the question in prose instead of emitting JSON
        # candidates. "dsl" exists because both other styles expand into
        # publication vocabulary, which is out of distribution against code.
        style = (prompt_style or "physics").strip().lower()
        if style == "dsl":
            self._system_prompt = DSL_REWRITER_SYSTEM_PROMPT
            self._few_shot_user = _DSL_FEW_SHOT_USER
            self._few_shot_assistant = _DSL_FEW_SHOT_ASSISTANT
        elif style == "generic":
            self._system_prompt = GENERIC_REWRITER_SYSTEM_PROMPT
            self._few_shot_user = _GENERIC_FEW_SHOT_USER
            self._few_shot_assistant = _GENERIC_FEW_SHOT_ASSISTANT
        else:
            self._system_prompt = REWRITER_SYSTEM_PROMPT
            self._few_shot_user = _FEW_SHOT_USER
            self._few_shot_assistant = _FEW_SHOT_ASSISTANT
        if cache_size and cache_size > 0:
            self._cached = functools.lru_cache(maxsize=cache_size)(
                self._rewrite_multi_uncached,
            )
        else:
            self._cached = self._rewrite_multi_uncached

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def rewrite_multi(self, question: str) -> list[Candidate]:
        """Return 1..max_candidates candidates. Never empty (fallback
        guarantees at least one)."""
        q = (question or "").strip()
        if not q:
            # Empty question: return a single trivial candidate rather
            # than raise, so downstream code has a stable contract.
            return [Candidate(dense_text="", sparse_text="", label="empty")]
        # lru_cache returns a tuple so it stays hashable; unpack.
        return list(self._cached(q))

    def __call__(self, question: str) -> tuple[str, str]:
        """Backward-compat: yield the top candidate as (dense, sparse)."""
        cands = self.rewrite_multi(question)
        return cands[0].as_pair()

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------

    def _rewrite_multi_uncached(self, question: str) -> tuple[Candidate, ...]:
        try:
            payload = self._call_llm(question)
        except Exception:
            return self._with_original(question, self._fallback_candidates(question))
        parsed = self._parse_payload(payload)
        if parsed is None:
            return self._with_original(question, self._fallback_candidates(question))
        cands = self._candidates_from_parsed(parsed, original=question)
        if not cands:
            return self._with_original(question, self._fallback_candidates(question))
        return self._with_original(question, tuple(cands[: self.max_candidates]))

    def _with_original(
        self, question: str, cands: tuple[Candidate, ...],
    ) -> tuple[Candidate, ...]:
        """Ensure the raw user question is the first candidate. Guards
        against rewrites that drop critical keywords (see 5ade9c9c
        Carfagno probe: `LLM dropped 'screenplay'` → gold recall crashed
        from top-2 to MISS)."""
        if not self.keep_original:
            return cands
        original = (question or "").strip()
        if not original:
            return cands
        # Dedup: skip if any existing cand already uses this exact dense_text
        for c in cands:
            if c.dense_text == original:
                return cands
        head = Candidate(
            dense_text=original,
            sparse_text=original,
            label="original",
        )
        return (head,) + cands[: max(0, self.max_candidates - 1)]

    def _call_llm(self, question: str) -> str:
        msg = self.llm.chat_with_tools(
            [
                {"role": "system", "content": self._system_prompt},
                {"role": "user", "content": self._few_shot_user},
                {"role": "assistant", "content": self._few_shot_assistant},
                {"role": "user", "content": question},
            ],
            tools=[],
            tool_choice="none",
        )
        return (getattr(msg, "content", "") or "").strip()

    @staticmethod
    def _parse_payload(text: str) -> Optional[dict]:
        if not text:
            return None
        m = _JSON_RE.search(text)
        if not m:
            return None
        try:
            obj = json.loads(m.group(0))
        except json.JSONDecodeError:
            return None
        if not isinstance(obj, dict):
            return None
        return obj

    def _candidates_from_parsed(
        self, parsed: dict, *, original: str,
    ) -> list[Candidate]:
        out: list[Candidate] = []
        seen_dense: set[str] = set()
        raw_candidates = parsed.get("candidates")
        if not isinstance(raw_candidates, list):
            raw_candidates = []
        for entry in raw_candidates:
            if not isinstance(entry, dict):
                continue
            query = (entry.get("query") or "").strip()
            if not query:
                continue
            expansions_raw = entry.get("expansions") or []
            if not isinstance(expansions_raw, list):
                expansions_raw = []
            expansions = [str(x).strip() for x in expansions_raw if str(x).strip()]
            label = str(entry.get("label") or "").strip() or "rewrite"
            dense_text = normalize_unicode_physics(query)
            # Sparse text: query + expansions + verbatim named entities
            # from the ORIGINAL question so BM25 keeps entity spelling
            # even if the LLM paraphrased it. Order preserved.
            sparse_bag = _dedup_tokens(
                _tokenise(dense_text)
                + expansions
                + _tokenise(normalize_unicode_physics(original))
            )
            sparse_text = " ".join(sparse_bag).strip() or dense_text
            if dense_text in seen_dense:
                continue
            seen_dense.add(dense_text)
            out.append(Candidate(
                dense_text=dense_text,
                sparse_text=sparse_text,
                label=label,
                expansions=expansions,
            ))
        return out

    def _fallback_candidates(self, question: str) -> tuple[Candidate, ...]:
        dense, sparse = self.fallback(question)
        return (Candidate(
            dense_text=dense, sparse_text=sparse, label="fallback",
        ),)


# ---------------------------------------------------------------------------
# tokenisation helpers — deliberately dumb, we're building BM25 input not
# doing NLP. Split on whitespace, keep punctuation that BM25 tokenisers
# usually retain in physics contexts (parens, +, -, /, .).
# ---------------------------------------------------------------------------

_TOKEN_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9\+\-\(\)/\.]*")


def _tokenise(text: str) -> list[str]:
    return _TOKEN_RE.findall(text or "")


def _dedup_tokens(tokens: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for t in tokens:
        key = t.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(t)
    return out
