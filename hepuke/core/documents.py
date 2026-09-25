"""Document helpers — JSON loader + generic dict-to-Document coercion.

Uses `hepuke.document.Document` (native dataclass) instead of the
llama-index `Document`. The old `_depth_first_yield` helper is inlined here
to avoid the private llama-index dependency.
"""
from __future__ import annotations

import json
from glob import glob
from typing import Any, Iterator, List

from hepuke.document import Document


def _depth_first_yield(
    payload: Any,
    *,
    levels_back: int = 1,
    path: list[str] | None = None,
    ensure_ascii: bool = False,
) -> Iterator[str]:
    """Yield leaf lines of the form 'k1 k2 ... kn: value'.

    Replacement for llama-index `readers.json._depth_first_yield` (private API).
    """
    path = list(path or [])
    if isinstance(payload, dict):
        for k, v in payload.items():
            yield from _depth_first_yield(
                v, levels_back=levels_back, path=path + [str(k)], ensure_ascii=ensure_ascii
            )
    elif isinstance(payload, list):
        for v in payload:
            yield from _depth_first_yield(
                v, levels_back=levels_back, path=path, ensure_ascii=ensure_ascii
            )
    else:
        key_path = " ".join(path[-levels_back:]) if levels_back else ""
        val = payload if isinstance(payload, str) else json.dumps(
            payload, ensure_ascii=ensure_ascii
        )
        yield f"{key_path}: {val}" if key_path else str(val)


def json_to_documents(path: str) -> List[Document]:
    with open(path, encoding="utf-8") as f:
        payload = json.load(f)
    docs: List[Document] = []
    items = payload if isinstance(payload, list) else [payload]
    for item in items:
        lines = list(_depth_first_yield(item, levels_back=1))
        text = "\n".join(lines)
        docs.append(Document(text=text, metadata={"file": path}))
    return docs


def load_json_dir(directory: str, pattern: str = "*.json") -> List[Document]:
    documents: List[Document] = []
    for f in glob(f"{directory}/{pattern}", recursive=True):
        documents.extend(json_to_documents(f))
    return documents


def dict_to_document(data: dict[str, Any]) -> Document:
    """Coerce a REST payload {text, metadata, ...} into a Document."""
    return Document.from_dict(data)
