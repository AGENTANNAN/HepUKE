"""14-slot extractor for BOSS DSL scripts.

Each slot corresponds to a single methodological decision a physicist would
list as an independent line in a systematic-uncertainty table. Values are
bucketed to make L-way JSD meaningful under L=5 replicates.

Special values:
    MISSING  — no match at all (model did not touch this decision)
    NONE     — model explicitly wrote "no cut" / "no fit" / equivalent
    N/A      — the slot does not apply for this analysis class

    python token_conf/extract_slots.py path/to/file.rb
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

# ---------------------------------------------------------------------------
# helpers

_NUM = r"[-+]?\d+(?:\.\d+)?"


def _strip_comments(src: str) -> str:
    """Remove Ruby `# ...` line comments (outside strings) to avoid regex hits
    inside them. STRING BODIES are preserved — they carry semantic content
    like dataset ids ("712_3773") and heredoc decay cards.
    """
    out = []
    for line in src.splitlines():
        in_str = False
        quote = ""
        buf = []
        i = 0
        while i < len(line):
            ch = line[i]
            if in_str:
                buf.append(ch)
                if ch == "\\" and i + 1 < len(line):
                    buf.append(line[i + 1]); i += 2; continue
                if ch == quote:
                    in_str = False
                i += 1; continue
            if ch in ('"', "'"):
                in_str = True; quote = ch; buf.append(ch); i += 1; continue
            if ch == "#":
                break  # rest of line is comment
            buf.append(ch); i += 1
        out.append("".join(buf))
    return "\n".join(out)


def _bucket(x: float, edges: list[float]) -> str:
    """Right-open buckets: return "<e0" / "[e0,e1)" / ... / ">=eN"."""
    for i, e in enumerate(edges):
        if x < e:
            return "<" + str(e) if i == 0 else f"[{edges[i-1]},{e})"
    return f">={edges[-1]}"


# ---------------------------------------------------------------------------
# per-slot extractors  (each returns a normalized string token)

def slot_dataset_id(src: str) -> str:
    m = re.findall(r'\.find\(\s*"(\d{3}_\d{3,4})"', src)
    if not m:
        return "MISSING"
    # if multiple datasets, sort+join so order does not fake divergence
    return "+".join(sorted(set(m)))


def slot_exclusive_mc_events(src: str) -> str:
    m = re.search(r"config\.events\s*=\s*([\d_]+)", src)
    if not m:
        return "NONE"
    n = int(m.group(1).replace("_", ""))
    for cutoff, label in [(200_000, "<200k"), (300_000, "~200k"),
                          (700_000, "~500k"), (2_000_000, ">=1M")]:
        if n < cutoff:
            return label
    return ">=2M"


def slot_decay_card_topology(src: str) -> str:
    """Return normalized 'mother→sorted daughters' of the FIRST Decay block."""
    heredoc = re.search(r"<<[~-]?DECAYCARD\s*(.*?)\bDECAYCARD\b", src, re.S)
    body = heredoc.group(1) if heredoc else src
    m = re.search(r"Decay\s+([^\s]+)\s*\n\s*[\d.]+\s+([^;]+)", body)
    if not m:
        return "MISSING"
    mother = m.group(1).strip()
    prods = [p for p in re.split(r"\s+", m.group(2).strip()) if p and p != "PHSP"]
    prods = [p for p in prods if not re.fullmatch(r"[A-Z_]+", p)]  # drop models
    return f"{mother}→{'+'.join(sorted(prods))}"


def slot_track_costheta_cut(src: str) -> str:
    m = re.search(rf"cos.?theta[^<>=]*[<=]\s*({_NUM})", src, re.I)
    if not m:
        return "MISSING"
    return f"{float(m.group(1)):.2f}"


def slot_track_Vz_cut(src: str) -> str:
    # common notations: Vz < 100, |Vz|<10 cm, config.vz_cut = 10
    m = re.search(rf"\bV[_z]?z?\b[^<>=\n]{{0,20}}[<=]\s*({_NUM})", src)
    if not m:
        return "MISSING"
    x = float(m.group(1))
    for cutoff, label in [(5, "<5"), (12, "10cm"), (25, "20cm"), (150, "100cm")]:
        if x < cutoff:
            return label
    return ">=150"


def slot_photon_energy_barrel_MeV(src: str) -> str:
    # look for barrel-energy min: number near "barrel" and "MeV" / "energy"
    m = re.search(
        rf"barrel[^\n]{{0,60}}?({_NUM})\s*MeV|energy[^\n]{{0,60}}?({_NUM})\s*MeV",
        src, re.I,
    )
    if not m:
        return "MISSING"
    val = float(m.group(1) or m.group(2))
    for cutoff, label in [(30, "25"), (45, "40"), (60, "50")]:
        if val < cutoff:
            return label
    return ">=60"


def slot_photon_track_isolation_deg(src: str) -> str:
    m = re.search(rf"(?:iso|angle|isolat)[^\n]{{0,60}}?({_NUM})\s*(?:deg|°)?", src, re.I)
    if not m:
        return "MISSING"
    x = float(m.group(1))
    for cutoff, label in [(12, "10"), (25, "20"), (35, "30")]:
        if x < cutoff:
            return label
    return ">=35"


def slot_pid_method(src: str) -> str:
    s = src.lower()
    if "probability" in s and "pid" in s:
        return "probability"
    if re.search(r"chi.?2.*(dedx|tof|combined)", s):
        return "chi2_combined"
    if "tof" in s and "dedx" in s:
        return "tof_dedx_chi2"
    if "pid" in s or "particle_id" in s:
        return "pid_other"
    return "NONE"


def slot_KS_mass_window_MeV(src: str) -> str:
    m = re.search(rf"M[^\n]{{0,10}}K.?S[^\n]{{0,30}}?({_NUM})\s*MeV", src, re.I)
    if not m:
        return "MISSING"
    x = float(m.group(1))
    for cutoff, label in [(12, "10"), (18, "15"), (25, "20"), (60, "50")]:
        if x < cutoff:
            return label
    return ">=60"


def slot_KS_Lsigma_cut(src: str) -> str:
    m = re.search(rf"L\s*/\s*(?:sigma|σ)[_L]?\s*[><=]\s*({_NUM})", src, re.I)
    if m:
        return f"L/σ>{float(m.group(1)):.1f}"
    if re.search(r"K.?S.?0?.*(?:vertex_fit|VertexFit|VertexAlg)", src, re.I):
        return "vertex_fit_no_cut"
    return "NONE"


def slot_kinematic_fit_type(src: str) -> str:
    kinds = set()
    if re.search(r"constrain_four_momentum|4C|four.?momentum", src, re.I):
        kinds.add("4C")
    if re.search(r"\b5C\b|kalman.*pi0|pi0.*kalman", src, re.I):
        kinds.add("5C")
    if re.search(r"kalman", src, re.I) and "5C" not in kinds:
        kinds.add("Kalman")
    if not kinds:
        return "NONE"
    return "+".join(sorted(kinds))


def slot_kfit_chi2_cut(src: str) -> str:
    m = re.findall(rf"chi.?2[^<>=\n]{{0,20}}[<=]\s*({_NUM})", src, re.I)
    if not m:
        m = re.findall(rf"chi2_cut\s+({_NUM})", src)
    if not m:
        return "MISSING"
    x = min(float(v) for v in m)   # tightest cut is the binding one
    for cutoff, label in [(30, "<30"), (55, "50"), (110, "100"), (210, "200")]:
        if x < cutoff:
            return label
    return ">=210"


def slot_analysis_class(src: str) -> str:
    if re.search(r"\bTagAnalysis\.new\b", src):
        return "TagAnalysis"
    if re.search(r"\bPartialRec\.new\b|\bPartialReconstruction\b", src, re.I):
        return "PartialRec"
    if re.search(r"\bSelection\.new\b", src):
        return "Selection"
    if re.search(r"\bConExc\.new\b", src):
        return "ConExc"
    return "other"


def slot_tag_side_count(src: str) -> str:
    return str(len(re.findall(r"\.tag_side\b", src)))


SLOTS = {
    "dataset_id": slot_dataset_id,
    "exclusive_mc_events": slot_exclusive_mc_events,
    "decay_card_topology": slot_decay_card_topology,
    "track_costheta_cut": slot_track_costheta_cut,
    "track_Vz_cut": slot_track_Vz_cut,
    "photon_energy_barrel_MeV": slot_photon_energy_barrel_MeV,
    "photon_track_isolation_deg": slot_photon_track_isolation_deg,
    "pid_method": slot_pid_method,
    "KS_mass_window_MeV": slot_KS_mass_window_MeV,
    "KS_Lsigma_cut": slot_KS_Lsigma_cut,
    "kinematic_fit_type": slot_kinematic_fit_type,
    "kfit_chi2_cut": slot_kfit_chi2_cut,
    "analysis_class": slot_analysis_class,
    "tag_side_count": slot_tag_side_count,
}


def extract(src: str) -> dict[str, str]:
    stripped = _strip_comments(src)
    return {name: fn(stripped) for name, fn in SLOTS.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="+")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    for p in args.path:
        src = Path(p).read_text()
        slots = extract(src)
        if args.json:
            print(json.dumps({"path": p, "slots": slots}, ensure_ascii=False))
        else:
            print(f"# {p}")
            for k, v in slots.items():
                print(f"  {k:<30} {v}")
            print()


if __name__ == "__main__":
    main()
