"""gwsignal4gwsurr package initialization."""

import os
import sys
import subprocess
import platform


def _get_git_info():
    """Read git info directly from the repo at import time."""
    pkg_dir = os.path.dirname(os.path.abspath(__file__))

    def run(cmd):
        try:
            return subprocess.check_output(
                cmd, cwd=pkg_dir, text=True, stderr=subprocess.DEVNULL
            ).strip()
        except Exception:
            return None

    git_hash = run(["git", "rev-parse", "--short=9", "HEAD"])
    if git_hash is None:
        return None

    describe = run(["git", "describe", "--tags", "--always", "--dirty"])
    status = run(["git", "status", "--porcelain"])
    dirty_files = status.splitlines() if status else []

    return {
        "git_hash": git_hash,
        "describe": describe,
        "dirty_files": dirty_files,
        "dirty": len(dirty_files) > 0,
    }

_git_info = _get_git_info()

if _git_info:
    __version__ = _git_info["describe"]
else:
    try:
        from importlib.metadata import version
        __version__ = version("gwsignal4gwsurr")
    except Exception:
        __version__ = "unknown"

def _print_diagnostics():
    line = f"gwsignal4gwsurr {__version__} | python {sys.version.split()[0]} | {platform.platform()}"
    print(line)
    if _git_info and _git_info["dirty"]:
        for f in _git_info["dirty_files"]:
            print(f"    {f.lstrip()}")

_print_diagnostics()
