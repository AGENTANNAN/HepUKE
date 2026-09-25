"""Native Document dataclass (replaces `llama_index.core.Document`)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class Document:
    text: str
    metadata: dict[str, Any] = field(default_factory=dict)
    doc_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {"text": self.text, "metadata": dict(self.metadata), "doc_id": self.doc_id}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Document":
        if not data.get("text"):
            raise ValueError("Document dict must contain non-empty 'text' field.")
        return cls(
            text=data["text"],
            metadata=dict(data.get("metadata") or {}),
            doc_id=data.get("doc_id"),
        )
