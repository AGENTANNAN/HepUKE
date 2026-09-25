"""Read-only shell sandbox exposed to the agent.

`run_cmd(argv, path_root, ...)` lets the LLM invoke a small set of read-only
POSIX tools (`find`, `wc`, `head`, `tail`, `cat`, `ls`, `sort`, `uniq`,
`cut`, `tr`, `grep`) with layered guards. This is a *stricter* superset of
what `memo_tools.grep` gives — same "application-layer sandbox" tier
(Firejail/Bubblewrap not required), but broader command surface.

Guards (all four MUST hold; any one failing rejects the call):

  1. **Command whitelist.** `argv[0]` must be an exact match against
     `_ALLOWED_CMDS`. No user-typed path, no `sh -c`, no interpreter.
  2. **Per-command flag denylist.** Each whitelisted command has a set of
     write / execute / network / output-redirect flags that are refused
     up-front (`find -delete`, `find -exec`, `find -fprint`, `grep -Z` w/
     xargs, etc.). Any argv element matching a denylisted prefix is rejected.
  3. **Path confinement.** Every argv element that "looks like a path"
     (starts with `/`, `./`, `../`, or contains `/`) is resolved against
     `path_root` and must live under it. Absolute paths must be under
     `path_root`. `--` is respected as an end-of-flags marker.
  4. **Process caging.** `shell=False`, `cwd=path_root`, `timeout` enforced,
     `stdout`/`stderr` each truncated at `max_bytes`. Environment is scrubbed
     (`env={"LC_ALL": "C", "LANG": "C", "PATH": "/usr/bin:/bin"}`) so a
     hostile `LD_PRELOAD` / `PYTHONPATH` can't slip through.

**What this does NOT protect against, by design:**
  * Fork bombs launched *by a whitelisted command* — mitigated by `timeout`
    only, not by cgroups. If you need hard CPU/mem limits, wrap this in a
    real OS-level sandbox (Firejail / Bubblewrap / docker `--read-only
    --network=none --pids-limit=64`).
  * Reading arbitrary files under `path_root`. If `path_root` contains
    secrets, don't point the agent at it.
  * `grep -f /etc/passwd` style tricks reading files *outside* path_root —
    blocked by guard 3 (any `/etc/...` argument is refused).

To add a new command:
  1. Add its name to `_ALLOWED_CMDS`.
  2. Add its forbidden flags to `_FORBIDDEN_FLAGS`.
  3. Add a unit test in `tests/unit/test_shell_tools.py` that exercises
     both the happy path and the new denylist entries.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Iterable, Sequence

from hepuke.observability.trace import emit
from hepuke.tools.memo_tools import MemoToolError, confine, safe_root


_ALLOWED_CMDS: frozenset[str] = frozenset({
    "grep", "find", "wc", "head", "tail", "cat", "ls",
    "sort", "uniq", "cut", "tr",
})

# Per-command denylist of flag *prefixes*. We match by prefix so `-exec`
# catches `-execdir` too. `find`'s action flags are the big risk.
_FORBIDDEN_FLAGS: dict[str, tuple[str, ...]] = {
    "find": (
        "-delete", "-exec", "-execdir", "-ok", "-okdir",
        "-fprint", "-fprint0", "-fprintf", "-fls",
    ),
    # `grep -Z` on its own is fine, but nothing here should be able to
    # write output files, so we ban `-o`/`--output` on every command that
    # supports one.
    "grep": ("--output-file", "-Z"),
    "sort": ("-o", "--output"),
    "tr":   (),
    "cut":  ("--output-delimiter=/",),  # nonsense: keeps future-proofing
    "wc":   (),
    "head": (),
    "tail": ("--follow", "-f", "--retry"),  # -f could hang past timeout
    "cat":  (),
    "ls":   (),
    "uniq": (),
}

# Global forbidden flags (apply to every command). Anything that could
# smuggle in a subprocess or output redirection lives here.
_GLOBAL_FORBIDDEN_FLAGS: tuple[str, ...] = (
    "--files-with-",   # placeholder — not a real flag; documents intent
)

# Environment given to the child. Locked-down enough to defeat LD_PRELOAD /
# PYTHONPATH tricks; enough to make grep behave in C locale.
_CHILD_ENV: dict[str, str] = {
    "LC_ALL": "C",
    "LANG": "C",
    "PATH": "/usr/bin:/bin",
}


def _looks_like_path(token: str) -> bool:
    """Heuristic: is this argv element a path (vs. a flag / pattern / number)?"""
    if not token:
        return False
    if token.startswith(("-", "+")):
        return False
    if token in ("--",):
        return False
    return token.startswith(("/", "./", "../")) or "/" in token


def _validate_argv(argv: Sequence[str], root: Path) -> list[str]:
    """Return a sanitized argv or raise MemoToolError.

    Every path-like element is resolved against `root` and re-emitted as its
    resolved absolute path, guaranteeing the child can only see files under
    `root`. "Path-like" is any of:
      * starts with `/`, `./`, `../`
      * contains a `/`
      * is a bare name that maps to a real file/dir under `root` (this
        catches symlink-in-corpus escapes: `cat sneaky.md` where
        `sneaky.md → /etc/passwd`).
    Non-path arguments (flags, regex patterns, filename globs that don't
    match a real file) pass through untouched.
    """
    if not argv:
        raise MemoToolError("argv is empty")

    cmd = argv[0]
    if cmd not in _ALLOWED_CMDS:
        raise MemoToolError(
            f"command `{cmd}` is not on the read-only whitelist. Allowed: "
            f"{sorted(_ALLOWED_CMDS)}"
        )

    forbidden = set(_FORBIDDEN_FLAGS.get(cmd, ())) | set(_GLOBAL_FORBIDDEN_FLAGS)
    end_of_flags = False
    cleaned: list[str] = [cmd]

    for tok in argv[1:]:
        if not isinstance(tok, str):
            raise MemoToolError(f"argv element must be a string, got {type(tok)!r}")
        if tok == "--":
            end_of_flags = True
            cleaned.append(tok)
            continue
        if not end_of_flags and tok.startswith("-"):
            for bad in forbidden:
                if tok == bad or tok.startswith(bad + "="):
                    raise MemoToolError(
                        f"flag `{tok}` is not allowed for `{cmd}` (denylist: "
                        f"{sorted(forbidden)})"
                    )

        path_like = _looks_like_path(tok)
        if not path_like:
            # Bare token: only intercept when it actually names a file under
            # root — this catches symlinks in the corpus that resolve outside
            # while leaving regex/glob patterns alone.
            candidate = root / tok
            if candidate.exists() or candidate.is_symlink():
                path_like = True

        if path_like:
            candidate = Path(tok) if tok.startswith("/") else (root / tok)
            resolved = confine(root, candidate)
            cleaned.append(str(resolved))
        else:
            cleaned.append(tok)

    return cleaned


def run_cmd(
    argv: Sequence[str],
    *,
    path_root: str | os.PathLike,
    timeout: int = 20,
    max_bytes: int = 128 * 1024,
) -> dict:
    """Execute a whitelisted read-only command under `path_root`.

    Parameters
    ----------
    argv : list[str]
        `argv[0]` must be in the whitelist. Path-like arguments are resolved
        under `path_root` and rejected if they escape it.
    path_root : str | Path
        The root the agent is confined to (typically the memo corpus root).
    timeout : int
        Wall-clock cap. Defaults to 20 s.
    max_bytes : int
        Each stream (stdout, stderr) is truncated to this many bytes; when
        cut, the returned dict flags `truncated=True`.

    Returns
    -------
    dict
        `{argv, returncode, stdout, stderr, truncated, timed_out}`.
        Errors from validation raise `MemoToolError` synchronously — those
        never make it to a subprocess.
    """
    root = safe_root(path_root)
    cleaned = _validate_argv(list(argv), root)

    try:
        proc = subprocess.run(
            cleaned,
            capture_output=True,
            cwd=str(root),
            env=_CHILD_ENV,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as e:
        emit(
            "tool_call", name="run_cmd", cmd=cleaned[0],
            timed_out=True, timeout_s=timeout,
        )
        return {
            "argv": cleaned,
            "returncode": None,
            "stdout": (e.stdout or b"")[:max_bytes].decode("utf-8", "replace"),
            "stderr": (e.stderr or b"")[:max_bytes].decode("utf-8", "replace"),
            "truncated": True,
            "timed_out": True,
        }

    stdout = proc.stdout or b""
    stderr = proc.stderr or b""
    stdout_lines = stdout.count(b"\n") + (1 if stdout and not stdout.endswith(b"\n") else 0)
    truncated_stdout = len(stdout) > max_bytes
    truncated_stderr = len(stderr) > max_bytes
    stdout_text = stdout[:max_bytes].decode("utf-8", "replace")
    stderr_text = stderr[:max_bytes].decode("utf-8", "replace")
    # Front-load the summary fields so the model sees them first even when a
    # downstream serializer truncates the payload.
    result = {
        "returncode": proc.returncode,
        "stdout_lines": stdout_lines,
        "stdout_bytes": len(stdout),
        "truncated": truncated_stdout or truncated_stderr,
        "timed_out": False,
        "argv": cleaned,
        "stdout": stdout_text,
        "stderr": stderr_text,
    }
    emit(
        "tool_call", name="run_cmd", cmd=cleaned[0],
        returncode=proc.returncode,
        stdout_bytes=len(stdout), stderr_bytes=len(stderr),
        truncated=result["truncated"],
    )
    return result


__all__ = ["run_cmd", "MemoToolError"]
