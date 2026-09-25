"""Tool registry — OpenAI function-calling shape.

Every tool has a name, a description, a JSON Schema for its arguments, and a
Python callable. `as_openai_tools()` emits the `tools=` payload the OpenAI
`chat.completions.create` API expects. `dispatch()` decodes JSON arguments
and invokes the callable with keyword arguments (never positional — models
sometimes reorder JSON keys).
"""
from __future__ import annotations

import inspect
import json
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable


@dataclass
class Tool:
    name: str
    description: str
    parameters: dict  # JSON Schema object
    call: Callable[..., Any]

    def as_openai(self) -> dict:
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description,
                "parameters": self.parameters,
            },
        }


class ToolRegistry:
    def __init__(self) -> None:
        self._tools: dict[str, Tool] = {}

    def register(self, tool: Tool) -> None:
        self._tools[tool.name] = tool

    def register_all(self, tools: Iterable[Tool]) -> None:
        for t in tools:
            self.register(t)

    def has(self, name: str) -> bool:
        return name in self._tools

    def get(self, name: str) -> Tool:
        if name not in self._tools:
            raise KeyError(f"tool not registered: {name}")
        return self._tools[name]

    def names(self) -> list[str]:
        return list(self._tools.keys())

    def as_openai_tools(self) -> list[dict]:
        return [t.as_openai() for t in self._tools.values()]

    def dispatch(self, name: str, arguments: dict | str) -> Any:
        tool = self.get(name)
        args = json.loads(arguments) if isinstance(arguments, str) else (arguments or {})
        if not isinstance(args, dict):
            raise TypeError(f"tool arguments must be a JSON object, got {type(args)!r}")
        # Strip kwargs the callable doesn't accept (models occasionally add
        # decorative keys). Keep call kwargs positional-safe.
        sig = inspect.signature(tool.call)
        accepted = set(sig.parameters.keys())
        has_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())
        if not has_kwargs:
            args = {k: v for k, v in args.items() if k in accepted}
        return tool.call(**args)
