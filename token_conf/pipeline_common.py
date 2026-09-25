"""Shared helpers for the description->DSL->entropy batch pipeline.

- config: read HepAI base_url/api_key from ../config.yaml
- assemble_skill_context(): build the FULL description2DSL-BOSS skill prompt
  (SKILL.md + guidelines + full manual + all examples + dataset tables), so a
  batch over diverse papers sees the same complete skill a real invocation would
  — otherwise a missing manual section inflates the generated DSL's entropy.
- chat(): one OpenAI-compatible chat/completions call, optional logprobs.
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from pathlib import Path

import os
import sys
import yaml

REPO = Path(__file__).resolve().parents[1]
SKILL = Path(os.environ.get("DSL_SKILL_DIR", REPO / "third_party" / "description2DSL-BOSS"))
RB_DIR = Path(os.environ.get("DSL_CORPUS_DIR", REPO / "data" / "corpus" / "dsl"))

# batch output layout under token_conf/batch/
BATCH = REPO / "token_conf" / "batch"
DESC_DIR = BATCH / "descriptions"   # <stem>.txt   (reverse-engineered NL description)
RUN_DIR = BATCH / "runs"            # <stem>.json  (skill output + content_logprobs)
GOLD_DIR = BATCH / "gold"           # <stem>.json  (per-token entropy gold)


def _raw_cfg() -> dict:
    path = REPO / "config.yaml"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} not found — copy the template first:\n"
            f"    cp config.example.yaml config.yaml")
    return yaml.safe_load(path.read_text()) or {}


def load_cfg() -> tuple[str, str]:
    """(base_url, api_key) for the OpenAI-compatible chat endpoint.

    Resolved through ``hepuke.load_config`` so that ``${VAR}`` placeholders in
    config.yaml are expanded from the environment or .env, with a plain YAML
    read as a fallback.
    """
    base = key = ""
    try:
        sys.path.insert(0, str(REPO))
        from hepuke import load_config          # noqa: PLC0415
        c = load_config()
        base, key = str(c.llm.base_url or ""), str(c.llm.api_key or "")
    except Exception:
        llm = _raw_cfg().get("llm", {})
        base, key = str(llm.get("base_url") or ""), str(llm.get("api_key") or "")
    base = base.rstrip("/")
    if not base or base.startswith("${") or "your-openai-compatible-endpoint" in base:
        raise RuntimeError(
            "LLM endpoint is not configured. Set LLM_BASE_URL and LLM_API_KEY in "
            ".env (or edit llm.base_url / llm.api_key in config.yaml).\n"
            "    cp .env.example .env && $EDITOR .env\n"
            "Check with: python scripts/check_setup.py --live")
    return base, key


def model_id(frozen: str) -> str:
    """Chat model to use.

    Defaults to the identifier frozen in the paper, so a rerun on the original
    provider reproduces the reported setting. Override with the ``DSL_MODEL``
    environment variable, or by setting ``llm.model`` in config.yaml, when
    running against a different provider — the frozen identifier is not
    guaranteed to resolve anywhere else.
    """
    env = os.environ.get("DSL_MODEL")
    if env:
        return env
    m = ""
    try:
        sys.path.insert(0, str(REPO))
        from hepuke import load_config          # noqa: PLC0415
        m = str(load_config().llm.model or "")
    except Exception:
        try:
            m = str((_raw_cfg().get("llm") or {}).get("model") or "")
        except Exception:
            m = ""
    if m and not m.startswith("${"):
        return m
    return frozen


_SKILL_CTX_CACHE: dict[str, str] = {}

TEMPLATE_SKELETON = REPO / "token_conf" / "skill_templates" / "template_skeleton.md"

SKILL_MODES = ("full", "template", "no-example", "no-example-clean-manual")


_MANUAL_DEFAULT_REDACTIONS: list[tuple[str, str]] = [
    # ---- CppTemplate class docstring examples (lines ~1349-1363) ----
    ("cos_theta: 0.93,",                    "cos_theta: <float>,"),
    ("Vz: 15.0,",                            "Vz: <float cm>,"),
    ("Vr: 1.5",                              "Vr: <float cm>"),
    (":energyThreshold, 0.035",              ":energyThreshold, <float GeV>"),

    # ---- CppTemplate::TEMPLATES defaults, track/photon (lines ~1444-1460) ----
    (":cos_theta => 0.93,   # |cosθ| < 0.93",
     ":cos_theta => <float>,        # |cosθ| upper bound"),
    (":Vz => 10.0,          # Vz cut in cm",
     ":Vz        => <float cm>,     # Vz cut"),
    (":Vr => 1.0,           # Vr cut in cm",
     ":Vr        => <float cm>,     # Vr cut"),
    (":tdc_emc_start => 0,          # EMC timing start (ns)",
     ":tdc_emc_start     => <int ns>,   # EMC timing start"),
    (":tdc_emc_end => 14,           # EMC timing end (ns)",
     ":tdc_emc_end       => <int ns>,   # EMC timing end"),
    (":angle_to_track => 10.0,      # Minimum angle to tracks (degrees)",
     ":angle_to_track    => <float deg>, # Min angle to tracks"),
    (":energyThreshold_b => 0.025,  # Barrel energy threshold (GeV)",
     ":energyThreshold_b => <float GeV>, # Barrel energy threshold"),
    (":energyThreshold_e => 0.050   # Endcap energy threshold (GeV)",
     ":energyThreshold_e => <float GeV>  # Endcap energy threshold"),

    # ---- CppTemplate::TEMPLATES defaults, isolated photon / PID / lepton (lines ~1468-1500) ----
    (":angle_to_prpprm => 20.0,     # Angle to (anti)proton (degrees)",
     ":angle_to_prpprm => <float deg>, # Angle to (anti)proton"),
    (":prob_cut => 0.001,           # Minimum probability cut",
     ":prob_cut => <float>,         # Minimum probability cut"),
    (":momentum_cut => 1.0,         # = treat_as_lepton_if_momentum_above (GeV)",
     ":momentum_cut => <float GeV>, # = treat_as_lepton_if_momentum_above"),
    (":energy_cut => 0.6,           # = treat_as_electron_if_energy_above (GeV)",
     ":energy_cut => <float GeV>,   # = treat_as_electron_if_energy_above"),
    (":prob_cut => 0.001            # Probability cut",
     ":prob_cut => <float>          # Probability cut"),

    # ---- Selection class @example docstrings (lines ~1548-1614) ----
    ("    cos_theta     0.93   # |cosθ| < 0.93",
     "    cos_theta     <float>   # |cosθ| upper bound"),
    ("    Vz            10.0   # |Vz| < 10 cm (along beam axis)",
     "    Vz            <float cm>   # |Vz| cut (along beam axis)"),
    ("    Vr            1.0    # Vr < 1 cm (transverse plane)",
     "    Vr            <float cm>    # Vr cut (transverse plane)"),
    ("               tdc_emc_start   0 # the flight time in the TOF start at 0 (unit: 700 ns).",
     "               tdc_emc_start   <int ns> # EMC TDC start."),
    ("               tdc_emc_end   14  # end at 14 (unit: 700 ns).",
     "               tdc_emc_end   <int ns>  # EMC TDC end."),
    ("               energyThreshold_b 0.025 # The deposited energy of each shower must be more than 25 MeV in the barrel region (| cos theta| < 0.80)",
     "               energyThreshold_b <float GeV> # Barrel shower energy threshold."),
    ("               energyThreshold_e  0.050 # and more than 50 MeV in the end cap region (0.86 < | cos theta| < 0.92).",
     "               energyThreshold_e <float GeV> # Endcap shower energy threshold."),
    ("               angle_to_track   10.0   # angle to nearest charged track shower > 10 degrees",
     "               angle_to_track   <float deg>   # angle to nearest charged track shower."),
    ("               angle_to_prm_track   20.0  ",
     "               angle_to_prm_track   <float deg>"),
    ("               angle_to_prp_track   20.0  ",
     "               angle_to_prp_track   <float deg>"),
    ("    prob_cut 0.001",
     "    prob_cut <float>"),
    ("    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,\n  #                                    treat_as_electron_if_energy_above: 0.6",
     "    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: <float GeV>,\n  #                                    treat_as_electron_if_energy_above: <float GeV>"),
]


def _redact_manual_defaults(manual_text: str) -> str:
    """Strip the concrete BESIII default cut values from the dsl_manual.md
    CppTemplate::TEMPLATES block (lines ~1349–1460), preserving field names and
    keyword arguments so the DSL schema stays intact.

    Raises if any target literal is missing — the manual has drifted and this
    redaction must be re-audited before shipping.
    """
    out = manual_text
    missing: list[str] = []
    for src, dst in _MANUAL_DEFAULT_REDACTIONS:
        if src not in out:
            missing.append(src)
            continue
        out = out.replace(src, dst)
    if missing:
        raise RuntimeError(
            "dsl_manual.md redaction targets not found — audit the CppTemplate "
            "block before re-running:\n  " + "\n  ".join(repr(m) for m in missing)
        )
    return out


def assemble_skill_context(mode: str = "full") -> str:
    """Full skill package as one authoritative system-prompt body (cached per mode).

    mode:
      full                     — SKILL.md + guidelines + full manual + example1/2/3 + dataset tables
      template                 — same, but example1/2/3 replaced by a single structural skeleton
                                 (API signatures + field names, no concrete values / modes / cuts)
      no-example               — same, but example1/2/3 entirely removed (guidelines/manual only)
      no-example-clean-manual  — same as no-example, AND the CppTemplate defaults block
                                 in dsl_manual.md (lines ~1349–1460) has its concrete cut
                                 values (0.93, 10.0, 1.0, 0.025, 0.050, 15.0, 1.5, 0.035)
                                 replaced by <float>/<int ns>/... placeholders. Field
                                 names and API signatures are preserved.
    """
    if mode not in SKILL_MODES:
        raise ValueError(f"unknown skill mode: {mode!r} (want one of {SKILL_MODES})")
    if mode in _SKILL_CTX_CACHE:
        return _SKILL_CTX_CACHE[mode]
    if not (SKILL / "SKILL.md").is_file():
        raise FileNotFoundError(
            f"DSL generation prompt package not found at {SKILL}.\n"
            f"  It is a third-party artifact and is not redistributed with this\n"
            f"  release. Point DSL_SKILL_DIR at your own copy:\n"
            f"    export DSL_SKILL_DIR=/path/to/description2DSL-BOSS\n"
            f"  Only the DSL *generation* stages need it; scoring, entropy,\n"
            f"  gating and the controller do not.")

    def _read(name: str, path: Path) -> str:
        text = path.read_text()
        if mode == "no-example-clean-manual" and name == "references/dsl_manual.md (full)":
            text = _redact_manual_defaults(text)
        return text

    base_parts: list[tuple[str, Path]] = [
        ("SKILL.md", SKILL / "SKILL.md"),
        ("references/translation_guidelines.md", SKILL / "references/translation_guidelines.md"),
        ("references/dsl_manual/index.md", SKILL / "references/dsl_manual/index.md"),
        ("references/dsl_manual.md (full)", SKILL / "references/dsl_manual.md"),
    ]
    example_parts: list[tuple[str, Path]] = [
        ("references/example1.md", SKILL / "references/example1.md"),
        ("references/example2.md", SKILL / "references/example2.md"),
        ("references/example3.md", SKILL / "references/example3.md"),
    ]
    tail_parts: list[tuple[str, Path]] = [
        ("assets/BES3_dataset.md", SKILL / "assets/BES3_dataset.md"),
        ("assets/BES3_incMC.md", SKILL / "assets/BES3_incMC.md"),
    ]

    if mode == "full":
        parts = base_parts + example_parts + tail_parts
        body = "\n\n---\n\n".join(f"# {name}\n\n{_read(name, p)}" for name, p in parts)
    elif mode == "template":
        parts = base_parts + tail_parts
        body = "\n\n---\n\n".join(f"# {name}\n\n{_read(name, p)}" for name, p in parts)
        skel = TEMPLATE_SKELETON.read_text()
        body += (
            "\n\n---\n\n# references/example_skeleton.md\n\n"
            "The three worked examples (example1/2/3.md) are ABSENT from this run.\n"
            "In their place, this structural skeleton shows the required Ruby DSL\n"
            "shape (API signatures + field names only). Fill in physics-appropriate\n"
            "values yourself; do NOT copy the placeholders literally.\n\n" + skel
        )
    else:  # no-example, no-example-clean-manual
        parts = base_parts + tail_parts
        body = "\n\n---\n\n".join(f"# {name}\n\n{_read(name, p)}" for name, p in parts)
        tail = (
            "\n\n---\n\n# references/examples\n\n"
            "The three worked examples (example1/2/3.md) are ABSENT from this run.\n"
            "Rely solely on SKILL.md, the translation guidelines, and the DSL manual.\n"
        )
        if mode == "no-example-clean-manual":
            tail += (
                "\nAdditionally, the CppTemplate::TEMPLATES default-values block in\n"
                "dsl_manual.md has had its concrete BESIII cut values (|cosθ|, Vz,\n"
                "Vr, EMC timing, energy thresholds) replaced by typed placeholders.\n"
                "Infer physics-appropriate values from the target description.\n"
            )
        body += tail

    _SKILL_CTX_CACHE[mode] = body
    return body


def chat(
    *,
    system: str,
    user: str,
    model: str,
    base_url: str,
    api_key: str,
    max_tokens: int,
    temperature: float = 0.0,
    logprobs: bool = False,
    top_logprobs: int = 20,
    stream: bool | None = None,
    timeout: float = 600.0,
    extra_payload: dict | None = None,
) -> dict:
    """One chat/completions call.

    NOTE on logprobs + HepAI: the HepAI gateway silently drops logprobs on
    non-streaming responses (choices[0].logprobs == null). It returns the real
    top-K distribution ONLY when stream=true. So whenever logprobs are requested
    we force streaming and reassemble content / reasoning / logprobs from the SSE
    chunks. temp=0 still yields a genuine (non-degenerate) distribution here —
    unlike api.deepseek.com, which returns -9999 sentinels at temp=0.

    ``extra_payload`` is merged verbatim into the request body — use it to pass
    provider-specific knobs like ``{"reasoning_effort": "none"}`` (accepted by
    deepseek-v4.1-flash to skip the reasoning stream; measured 0 reasoning
    chunks in probing). Unknown keys are silently ignored by most gateways;
    verify on your target model before shipping.
    """
    if stream is None:
        stream = logprobs  # must stream to get logprobs on HepAI
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        "max_tokens": max_tokens,
        "temperature": temperature,
        "stream": stream,
    }
    if logprobs:
        payload["logprobs"] = True
        payload["top_logprobs"] = min(max(top_logprobs, 0), 20)
    if extra_payload:
        payload.update(extra_payload)
    req = urllib.request.Request(
        f"{base_url}/chat/completions",
        data=json.dumps(payload).encode(),
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            if not stream:
                data = json.loads(r.read())
                ch = data["choices"][0]
                return {
                    "answer": ch["message"].get("content") or "",
                    "reasoning": ch["message"].get("reasoning_content") or "",
                    "finish_reason": ch.get("finish_reason"),
                    "content_logprobs": (ch.get("logprobs") or {}).get("content") or [],
                    "reasoning_logprobs": (ch.get("logprobs") or {}).get("reasoning_content") or [],
                    "usage": data.get("usage"),
                }
            # ---- streaming: reassemble from SSE ----
            answer_parts: list[str] = []
            reasoning_parts: list[str] = []
            content_lp: list[dict] = []
            reasoning_lp: list[dict] = []
            finish_reason = None
            usage = None
            for raw in r:
                line = raw.decode("utf-8", "replace").strip()
                if not line.startswith("data:"):
                    continue
                data = line[len("data:"):].strip()
                if data == "[DONE]":
                    break
                try:
                    obj = json.loads(data)
                except json.JSONDecodeError:
                    continue
                if obj.get("usage"):
                    usage = obj["usage"]
                choices = obj.get("choices") or []
                if not choices:
                    continue
                ch = choices[0]
                delta = ch.get("delta") or {}
                if delta.get("content"):
                    answer_parts.append(delta["content"])
                if delta.get("reasoning_content"):
                    reasoning_parts.append(delta["reasoning_content"])
                lp = ch.get("logprobs") or {}
                if lp.get("content"):
                    content_lp.extend(lp["content"])
                if lp.get("reasoning_content"):
                    reasoning_lp.extend(lp["reasoning_content"])
                if ch.get("finish_reason"):
                    finish_reason = ch["finish_reason"]
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"HTTP {e.code}: {e.read()[:600].decode(errors='replace')}") from e
    return {
        "answer": "".join(answer_parts),
        "reasoning": "".join(reasoning_parts),
        "finish_reason": finish_reason,
        "content_logprobs": content_lp,
        "reasoning_logprobs": reasoning_lp,
        "usage": usage,
    }


def stems(limit: int | None = None, only: list[str] | None = None) -> list[str]:
    if only:
        return only
    xs = sorted(p.stem for p in RB_DIR.glob("*.rb"))
    return xs[:limit] if limit else xs


def strip_code_fence(text: str) -> str:
    """If the model wrapped output in ```ruby ... ```, return the inner code."""
    t = text.strip()
    if t.startswith("```"):
        lines = t.splitlines()
        if lines and lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].strip().startswith("```"):
            lines = lines[:-1]
        return "\n".join(lines)
    return text


# ---------------------------------------------------------------------------
# Description tiers
# ---------------------------------------------------------------------------
# The full ~200-word descriptions in batch/descriptions/<qid>.txt were reverse-
# engineered from the gold DSL and inevitably restate concrete cut values, PID
# thresholds and mass windows in natural language. That leaks selection knobs
# into the No-RAG baseline. bare_desc() collapses a full description to two
# facts — energy point + process arrow — so a No-RAG run has only the task's
# physics identity and must invent everything else. RAG / oracle prompts also
# read from bare_desc() so the only variable across policies is the reference.
# ---------------------------------------------------------------------------

_BARE_ENERGY_PAT = re.compile(r"(\d\.\d{2,4})\s*GeV")
_BARE_ARROW = "→"


def bare_desc(desc: str) -> str:
    """Strip a full description down to energy point + reaction arrow.

    Uses the FIRST "A → B" arrow in the description as the physics process,
    and the FIRST "x.xx GeV" as the energy point. Falls back gracefully when
    those markers are missing.
    """
    s = desc.strip()
    en = _BARE_ENERGY_PAT.search(s)
    en_s = en.group(0) if en else "the relevant energy point"
    if _BARE_ARROW in s:
        idx = s.index(_BARE_ARROW)
        left = s[:idx].rstrip()
        mother_m = re.search(
            r"([A-Za-zψΨηΦφρωπΞΛΣ'′*_0-9\+\-\(\)^{}\\/]+)\s*$", left
        )
        mother = mother_m.group(1) if mother_m else left.split()[-1]
        right = s[idx + 1 :].lstrip()
        stop_m = re.search(
            rf"({_BARE_ARROW}|,\s|;\s|\bwith\b|\busing\b|\bthrough\b|\bvia\b|\bin\b|\.$|\.\s)",
            right,
        )
        prod = (right[: stop_m.start()] if stop_m else right[:60]).strip().rstrip(".,;")
        process = f"{mother} → {prod}"
    else:
        first = re.split(r"\.(\s|$)", s, maxsplit=1)[0]
        process = re.split(
            r",|\bwith\b|\busing\b|\bthrough\b|\bvia\b", first, maxsplit=1
        )[0].strip()
    return (
        f"Analyze BESIII data collected at {en_s}.\n"
        f"Target process: {process}.\n"
        f"Produce a BOSS DSL script for this measurement."
    )
