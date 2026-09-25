"""Specialist subagents.

Each subagent is a *nested* `RagAgent` with a narrow toolset and a
corpus-specific system prompt, wrapped as a single `Tool` for the main agent
to call. Rationale: the main prompt no longer has to memorize every
corpus's conventions — those live inside the subagent that owns them.

Contract every subagent respects (input → output):

    input : {question: str, max_steps?: int, doc_id_hint?: str}
    output: {
        "answer":   str | None,   # None only when nothing salvageable
        "evidence": [ {source, doc_id?, path?, snippet?}, ... ],
        "dry":      bool,          # True ⇒ corpus had nothing relevant
        "steps":    int,           # inner tool-call rounds actually used
    }

The main agent uses `dry` to decide whether to try a sibling subagent
instead of concluding "no data".
"""
from hepuke.agent.subagents.base import make_subagent_tool  # noqa: F401
from hepuke.agent.subagents.memo_search import (  # noqa: F401
    build_memo_search_subagent,
)
from hepuke.agent.subagents.milvus_search import (  # noqa: F401
    build_milvus_search_subagent,
)
from hepuke.agent.subagents.hypernews_db import (  # noqa: F401
    build_hypernews_db_subagent,
)
