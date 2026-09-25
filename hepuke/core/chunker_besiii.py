"""BESIII BAM-memo section chunker.

BAM memos use a mix of Roman-numeral top-level headings ("III. EVENT
SELECTION", "VII. SYSTEMATIC UNCERTAINTIES") and idiosyncratic subsection
labels ("NUMBER OF ψ' EVENTS", "BOSS SIMULATION") that fall through the
generic ICLR-flavoured classifier in ``chunker.PaperSectionChunker``. This
subclass:

* strips a leading Roman-numeral prefix (``I.``, ``II.``, ``VII.``) before
  matching, so BAM headings share the same lookup path as Arabic-numbered
  ones,
* keeps ``systematics`` as its own ``section_type`` (per user decision) —
  the systematic-error questions in the benchmark need distinct SP-level
  targets from regular ``results``,
* adds BESIII-specific aliases for event selection, fitting/branching
  fraction, detector-and-MC, and BOSS/simulation headings,
* falls back to the base ``PaperSectionChunker._classify`` for anything
  that doesn't match the BESIII table (so summary/conclusion/etc. still
  work).
"""
from __future__ import annotations

import re
from typing import Any

from hepuke.core.chunker import PaperSectionChunker


_ROMAN_RE = re.compile(r"^\s*([IVX]+)\.\s+(.*)$")


# Order matters: earlier entries win. `systematics` must appear before
# `results` because "SYSTEMATIC UNCERTAINTIES" also contains no other
# keyword — but if we later add a broader results pattern it must not
# accidentally swallow "SYSTEMATIC ..." rows.
_BESIII_SECTION_PATTERNS: list[tuple[str, str, str]] = [
    (r"SYSTEMATIC\s+UNCERTAINT(Y|IES)?", "systematics", "Systematic Uncertainties"),
    (r"EVENT\s+SELECTION|SELECTION\s+CRITERIA", "method", "Event Selection"),
    (r"BRANCHING\s+FRACTION|CROSS\s+SECTION|FITTING\s+RESULTS?|"
     r"NUMBER\s+OF\s+.*\s+EVENTS", "results", "Results"),
    (r"DETECTOR\s+AND\s+(DATA|MC|MONTE\s*CARLO)|MONTE\s+CARLO|BOSS\s+SIMULATION",
     "background", "Detector and MC"),
    (r"BESIII\s+AND\s+BEPCII|THE\s+BESIII\s+EXPERIMENT", "background", "BESIII/BEPCII"),
    (r"BACKGROUND\s+ANALYSIS|BACKGROUND\s+STUDY", "background", "Background Analysis"),
]


class BesiiiMemoSectionParser(PaperSectionChunker):
    """PaperSectionChunker specialised for BESIII BAM memo headings."""

    def _classify(self, heading: str) -> tuple[str, str, str] | None:
        text = self._TRAILING_PAGE_NUM.sub("", heading).strip()

        # Strip leading Roman numeral prefix if any. The Roman "number" is
        # captured as the section_number so downstream SP-picker can still
        # emit (bam_id, "III") if the memo uses Roman numerals.
        m = _ROMAN_RE.match(text)
        if m is not None:
            sec_num = m.group(1)
            tail = m.group(2).strip()
        else:
            sec_num = ""
            tail = text

        for pat, sec_type, canonical in _BESIII_SECTION_PATTERNS:
            if re.search(pat, tail, re.IGNORECASE):
                return sec_type, sec_num, tail or canonical

        # Fall back to Arabic-numbered generic classifier. If we found a
        # Roman prefix we already stripped it; pass the tail down so the
        # generic classifier sees a clean heading.
        if m is not None:
            base = super()._classify(tail)
            if base is not None:
                sec_type, _base_num, sec_title = base
                return sec_type, sec_num, sec_title
            return None
        return super()._classify(text)
