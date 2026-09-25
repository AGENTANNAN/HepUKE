"""Observability package (trace + eval)."""
from hepuke.observability.trace import (  # noqa: F401
    clear_trace_context,
    configure_trace,
    emit,
    set_trace_context,
)
from hepuke.observability.sinks import (  # noqa: F401
    JsonlSink,
    SqliteSink,
    TraceSink,
    is_sqlite_path,
    iter_sqlite_events,
    known_run_ids,
)
from hepuke.observability.eval import (  # noqa: F401
    Case,
    eval_answer,
    eval_retrieval,
    load_cases,
)
from hepuke.observability.reports import (  # noqa: F401
    AnswerReport,
    RetrievalReport,
    write_report,
)

