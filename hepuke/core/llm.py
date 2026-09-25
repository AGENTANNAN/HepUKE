"""LLM chat client (OpenAI-compatible /chat/completions).

Adds `chat_with_tools()` for OpenAI-style function calling (agent loop backbone)
and `ask_json()` for JSON-mode responses. Plain `chat()` preserved for the
existing FastAPI RAG facade.
"""
from __future__ import annotations

import json
import math
from typing import Any, List

from openai import BadRequestError, OpenAI

_DEFAULT_SYSTEM_PROMPT = (
    "你是严谨的助手。请仅根据给定的上下文回答用户问题；"
    "如果上下文不足以回答，请明确说不知道，不要编造。"
)


class LLMClient:
    # Class-level default so `__new__`-constructed clients (tests) still
    # resolve `self.extra_body` without going through __init__.
    extra_body: dict | None = None

    def __init__(
        self,
        *,
        base_url: str,
        api_key: str,
        model: str,
        temperature: float = 0.1,
        max_tokens: int = 1024,
        extra_body: dict | None = None,
    ) -> None:
        self.client = OpenAI(base_url=base_url, api_key=api_key)
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.extra_body = extra_body
        self._logprobs_unsupported = False

    # ---- plain chat --------------------------------------------------------
    def chat(
        self,
        question: str,
        contexts: List[str] | None = None,
        history: List[dict] | None = None,
        system_prompt: str = _DEFAULT_SYSTEM_PROMPT,
    ) -> str:
        messages: list[dict] = [{"role": "system", "content": system_prompt}]
        if history:
            messages.extend(history)
        if contexts:
            joined = "\n\n---\n\n".join(contexts)
            user_content = f"上下文:\n{joined}\n\n问题: {question}"
        else:
            user_content = question
        messages.append({"role": "user", "content": user_content})

        resp = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            temperature=self.temperature,
            max_tokens=self.max_tokens,
            extra_body=self.extra_body,
        )
        return resp.choices[0].message.content or ""

    # ---- OpenAI-style function calling (agent loop) ------------------------
    def chat_with_tools(
        self,
        messages: list[dict],
        tools: list[dict],
        tool_choice: str | dict = "auto",
        *,
        temperature: float | None = None,
        max_tokens: int | None = None,
        logprobs: bool = False,
        top_logprobs: int | None = None,
    ) -> Any:
        """Return the raw `ChatCompletionMessage` (with .content and .tool_calls).

        Caller inspects `.tool_calls` to dispatch tools, then appends
        `{'role':'tool', 'tool_call_id': ..., 'content': ...}` and re-enters
        the loop.

        Emits an `llm_call` trace record with usage/cache fields so the
        outer agent loop can see prefix-cache hit rates without an extra
        wrapper. HepAI/deepseek surfaces `prompt_tokens_details.cached_tokens`
        and (as bonus) `std_amounts.input_cache_read` — both are captured
        when present. Cache misses on the first turn are expected; steady
        state should read ~2048 cached tokens per subsequent turn.

        When ``logprobs=True`` and ``top_logprobs>0``, the message's first
        content-token top-B distribution is stashed on
        ``message._hepuke_posterior`` as a probability vector (softmax of
        the returned log-probs). This feeds the confidence controller's
        answer posterior p_t without adding a second round-trip.
        """
        from hepuke.observability.trace import emit

        extra: dict[str, Any] = {}
        want_logprobs = bool(
            logprobs and top_logprobs and top_logprobs > 0
            and not self._logprobs_unsupported
        )
        if want_logprobs:
            extra["logprobs"] = True
            extra["top_logprobs"] = int(top_logprobs)

        # Some OpenAI-compatible upstreams (HepAI's hosted deepseek family)
        # only surface `logprobs` when `stream=True` — non-stream silently
        # returns `choices[0].logprobs = null`. Force stream in that case
        # and aggregate the chunks into a synthetic message with the same
        # shape as a non-stream response.
        stream_mode = want_logprobs

        try:
            if stream_mode:
                resp_stream = self.client.chat.completions.create(
                    model=self.model,
                    messages=messages,
                    tools=tools,
                    tool_choice=tool_choice,
                    temperature=self.temperature if temperature is None else temperature,
                    max_tokens=self.max_tokens if max_tokens is None else max_tokens,
                    stream=True,
                    extra_body=self.extra_body,
                    **extra,
                )
                resp = _aggregate_stream(resp_stream)
            else:
                resp = self.client.chat.completions.create(
                    model=self.model,
                    messages=messages,
                    tools=tools,
                    tool_choice=tool_choice,
                    temperature=self.temperature if temperature is None else temperature,
                    max_tokens=self.max_tokens if max_tokens is None else max_tokens,
                    extra_body=self.extra_body,
                    **extra,
                )
        except BadRequestError as exc:
            # Some upstream providers (e.g. HepAI's hosted deepseek) return 400
            # for `logprobs=True`. Retry once without logprobs and remember so
            # subsequent calls skip the parameter entirely — the confidence
            # controller degrades gracefully (empty posterior → entropy signal
            # unavailable that round).
            if "logprobs" in extra and _is_logprobs_unsupported(exc):
                self._logprobs_unsupported = True
                extra.pop("logprobs", None)
                extra.pop("top_logprobs", None)
                resp = self.client.chat.completions.create(
                    model=self.model,
                    messages=messages,
                    tools=tools,
                    tool_choice=tool_choice,
                    temperature=self.temperature if temperature is None else temperature,
                    max_tokens=self.max_tokens if max_tokens is None else max_tokens,
                    extra_body=self.extra_body,
                    **extra,
                )
            else:
                raise
        usage = getattr(resp, "usage", None)
        message = resp.choices[0].message
        if usage is not None:
            details = getattr(usage, "prompt_tokens_details", None)
            cached = getattr(details, "cached_tokens", None) if details else None
            completion = getattr(usage, "completion_tokens", None)
            # Stash on the message so RagAgent's token-budget accounting can
            # read it inline without a second round-trip / thread-local
            # side channel. Attribute name is private on purpose (leading
            # underscore, hepuke_ prefix) so it never collides with real
            # OpenAI SDK fields.
            try:
                setattr(message, "_hepuke_completion_tokens", completion)
            except Exception:
                pass
            emit(
                "llm_call",
                model=self.model,
                prompt_tokens=getattr(usage, "prompt_tokens", None),
                completion_tokens=completion,
                cached_tokens=cached,
            )

        # Posterior extraction: first content-token top-B logprobs → softmax.
        # Multi-token length normalisation is a no-op at B=1 token, so this
        # matches the paper's "length-normalized softmax over top-B" up to a
        # single-token approximation (see method §4.2, App.~B-sensitivity).
        if extra:
            posterior = _extract_top_b_posterior(resp)
            if posterior is not None:
                try:
                    setattr(message, "_hepuke_posterior", posterior)
                except Exception:
                    pass
        return message

    # ---- JSON mode --------------------------------------------------------
    def ask_json(
        self,
        prompt: str,
        schema: dict | None = None,
        *,
        system_prompt: str | None = None,
        temperature: float | None = None,
    ) -> dict:
        """Ask for a JSON object. `schema` is prepended as a hint (models
        that don't support structured outputs still tend to comply)."""
        sys_msg = system_prompt or (
            "You return ONLY valid JSON. No prose, no fences."
            + (
                f"\nExpected JSON schema:\n{json.dumps(schema, ensure_ascii=False)}"
                if schema
                else ""
            )
        )
        resp = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": sys_msg},
                {"role": "user", "content": prompt},
            ],
            temperature=self.temperature if temperature is None else temperature,
            max_tokens=self.max_tokens,
            response_format={"type": "json_object"},
            extra_body=self.extra_body,
        )
        content = resp.choices[0].message.content or "{}"
        try:
            return json.loads(content)
        except json.JSONDecodeError as e:
            raise ValueError(f"ask_json: model returned non-JSON: {content!r}") from e


def _extract_top_b_posterior(resp: Any) -> list[float] | None:
    """Read the first content-token's top-B logprobs → softmax → prob vector.

    OpenAI-compatible response shape (SDK object or dict):
        choices[0].logprobs.content[0].top_logprobs = [{"token": ..., "logprob": ...}, ...]

    Returns ``None`` when the field is absent (server didn't send logprobs)
    or malformed — the caller treats missing posterior as "no observation
    this round" and the controller skips the entropy signal for that step.
    """
    try:
        choice = resp.choices[0]
        logprobs_obj = getattr(choice, "logprobs", None) or (
            choice.get("logprobs") if isinstance(choice, dict) else None
        )
        if not logprobs_obj:
            return None
        content = getattr(logprobs_obj, "content", None) or (
            logprobs_obj.get("content") if isinstance(logprobs_obj, dict) else None
        )
        if not content:
            return None
        first = content[0]
        top = getattr(first, "top_logprobs", None) or (
            first.get("top_logprobs") if isinstance(first, dict) else None
        )
        if not top:
            return None
        logps: list[float] = []
        for row in top:
            lp = getattr(row, "logprob", None) or (
                row.get("logprob") if isinstance(row, dict) else None
            )
            if lp is None:
                continue
            logps.append(float(lp))
        if not logps:
            return None
        m = max(logps)
        exps = [math.exp(lp - m) for lp in logps]
        z = sum(exps)
        if z <= 0.0:
            return None
        return [e / z for e in exps]
    except Exception:
        return None


def _aggregate_stream(chunks) -> Any:
    """Fold an OpenAI streaming response into a non-stream response shape.

    We only need the fields downstream code touches: `choices[0].message`
    (with `.content` and optional `.tool_calls`), `choices[0].logprobs`
    (with `.content = [<first token's TopLogprob list>]`), and `usage`.

    HepAI's deepseek surfaces `logprobs` chunk-by-chunk; we accumulate the
    per-token `top_logprobs` list and keep it in the same order as
    OpenAI's non-stream response would (first content token first).
    """
    from types import SimpleNamespace

    content_parts: list[str] = []
    tool_calls_by_idx: dict[int, dict] = {}
    logprob_tokens: list[Any] = []
    usage_obj = None

    for ch in chunks:
        if not getattr(ch, "choices", None):
            u = getattr(ch, "usage", None)
            if u is not None:
                usage_obj = u
            continue
        choice = ch.choices[0]
        delta = getattr(choice, "delta", None)
        if delta is None:
            continue
        # content
        piece = getattr(delta, "content", None)
        if piece:
            content_parts.append(piece)
        # tool calls (streamed as deltas keyed by index)
        tc_deltas = getattr(delta, "tool_calls", None) or []
        for tc in tc_deltas:
            idx = getattr(tc, "index", 0) or 0
            slot = tool_calls_by_idx.setdefault(
                idx, {"id": None, "type": "function",
                       "function": {"name": None, "arguments": ""}},
            )
            if getattr(tc, "id", None):
                slot["id"] = tc.id
            fn = getattr(tc, "function", None)
            if fn is not None:
                if getattr(fn, "name", None):
                    slot["function"]["name"] = fn.name
                if getattr(fn, "arguments", None):
                    slot["function"]["arguments"] += fn.arguments
        # logprobs (per-token entries)
        lp = getattr(choice, "logprobs", None)
        if lp is not None:
            lp_content = getattr(lp, "content", None) or []
            for tok in lp_content:
                logprob_tokens.append(tok)
        # capture usage if it lands on a choice-carrying chunk
        u = getattr(ch, "usage", None)
        if u is not None:
            usage_obj = u

    # Build a message that mimics ChatCompletionMessage well enough for the
    # rest of the pipeline (which reads .content and .tool_calls).
    tool_calls_list: list[Any] = []
    for idx in sorted(tool_calls_by_idx):
        slot = tool_calls_by_idx[idx]
        tool_calls_list.append(SimpleNamespace(
            id=slot["id"],
            type=slot["type"],
            function=SimpleNamespace(
                name=slot["function"]["name"],
                arguments=slot["function"]["arguments"],
            ),
        ))
    message = SimpleNamespace(
        role="assistant",
        content="".join(content_parts) or None,
        tool_calls=tool_calls_list or None,
    )

    logprobs_obj = SimpleNamespace(content=logprob_tokens) if logprob_tokens else None
    choice = SimpleNamespace(index=0, message=message, logprobs=logprobs_obj,
                              finish_reason=None)
    return SimpleNamespace(choices=[choice], usage=usage_obj)


def _is_logprobs_unsupported(exc: BadRequestError) -> bool:
    """Heuristic: does this 400 mean 'logprobs not supported'?

    HepAI's upstream deepseek returns:
        {"detail": {"error_code": "UPSTREAM_INVALID_REQUEST",
                    "message": "<400> InternalError.Algo.InvalidParameter:
                                The parameters `logprobs` is not supported."}}
    Match on the string 'logprobs' in the error body — good enough for the
    fallback, and any other 400 still bubbles up.
    """
    try:
        body = str(getattr(exc, "response", None) and exc.response.text) or str(exc)
    except Exception:
        body = str(exc)
    return "logprobs" in body


def build_llm(cfg) -> LLMClient:
    return LLMClient(
        base_url=cfg.base_url,
        api_key=cfg.api_key,
        model=cfg.model,
        temperature=cfg.temperature,
        max_tokens=cfg.max_tokens,
    )
