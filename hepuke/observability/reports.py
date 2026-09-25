"""Report dataclasses with JSON + Markdown three-line-table serialisation."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class RetrievalReport:
    suite: str
    config: str
    n_cases: int
    ks: list[int]
    recall: dict[int, float]
    mrr: float
    ndcg: dict[int, float]
    trace_file: str = ""
    per_case: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": "retrieval",
            "suite": self.suite,
            "config": self.config,
            "n_cases": self.n_cases,
            "ks": self.ks,
            "recall": {str(k): v for k, v in self.recall.items()},
            "mrr": self.mrr,
            "ndcg": {str(k): v for k, v in self.ndcg.items()},
            "trace_file": self.trace_file,
            "per_case": self.per_case,
        }

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)

    def to_markdown(self) -> str:
        lines: list[str] = []
        lines.append(f"### Retrieval — {self.suite} / {self.config} (n={self.n_cases})")
        header = ["k"] + [f"Recall@{k}" for k in self.ks] + [f"nDCG@{k}" for k in self.ks] + ["MRR"]
        row = ["-"] \
            + [f"{self.recall[k]:.4f}" for k in self.ks] \
            + [f"{self.ndcg[k]:.4f}" for k in self.ks] \
            + [f"{self.mrr:.4f}"]
        lines.append("| " + " | ".join(header) + " |")
        lines.append("|" + "|".join(["---"] * len(header)) + "|")
        lines.append("| " + " | ".join(row) + " |")
        return "\n".join(lines)


@dataclass
class AnswerReport:
    suite: str
    config: str
    n_cases: int
    em: float
    f1: float
    per_case: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "kind": "answer",
            "suite": self.suite,
            "config": self.config,
            "n_cases": self.n_cases,
            "em": self.em,
            "f1": self.f1,
            "per_case": self.per_case,
        }

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)

    def to_markdown(self) -> str:
        lines: list[str] = []
        lines.append(f"### Answer — {self.suite} / {self.config} (n={self.n_cases})")
        lines.append("| EM | F1 |")
        lines.append("|---|---|")
        lines.append(f"| {self.em:.4f} | {self.f1:.4f} |")
        return "\n".join(lines)


def write_report(
    out_dir: str | Path,
    suite: str,
    config: str,
    stamp: str,
    retrieval: RetrievalReport | None,
    answer: AnswerReport | None,
    trace_file: str = "",
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "suite": suite,
        "config": config,
        "stamp": stamp,
        "trace_file": trace_file,
    }
    if retrieval is not None:
        retrieval.trace_file = trace_file
        payload["retrieval"] = retrieval.to_dict()
    if answer is not None:
        payload["answer"] = answer.to_dict()

    json_path = out / f"{suite}_{config}_{stamp}.json"
    md_path = out / f"{suite}_{config}_{stamp}.md"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    md_parts: list[str] = [f"# {suite} / {config} — {stamp}"]
    if trace_file:
        md_parts.append(f"trace: `{trace_file}`")
    if retrieval is not None:
        md_parts.append(retrieval.to_markdown())
    if answer is not None:
        md_parts.append(answer.to_markdown())
    md_path.write_text("\n\n".join(md_parts) + "\n", encoding="utf-8")

    return {
        "json": str(json_path),
        "markdown": str(md_path),
        "payload": payload,
    }
