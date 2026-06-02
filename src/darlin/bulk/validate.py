from __future__ import annotations

import os
import shutil
from pathlib import Path


class BulkInputError(Exception):
    """Raised when bulk CLI input validation fails."""


def require_readable_files(labels_and_paths: list[tuple[str, str | Path]]) -> None:
    """Require each path to exist and be a regular file (not a directory).

    Raises BulkInputError with a multi-line message listing all failures.
    """
    failures: list[str] = []
    for label, p in labels_and_paths:
        path = Path(p).expanduser()
        if path.is_file():
            continue
        if path.exists() and not path.is_file():
            reason = "not a regular file"
        else:
            reason = "not found"
        failures.append(f"  {label}: {p} ({reason})")
    if failures:
        msg = "Invalid or missing input files:\n" + "\n".join(failures)
        raise BulkInputError(msg)


def require_pear_executable(pear_path: str) -> None:
    """Ensure PEAR can be invoked: explicit path must exist as a file; bare name must be on PATH."""
    expanded = os.path.expanduser(pear_path)
    p = Path(expanded)
    explicit = p.is_absolute() or len(p.parts) > 1
    if explicit:
        if not p.is_file():
            raise BulkInputError(f"PEAR executable not found: {pear_path}")
    elif shutil.which(pear_path) is None:
        raise BulkInputError(f"PEAR executable not found in PATH: {pear_path}")


def require_valid_locus(locus: str) -> None:
    try:
        from darlin_core.config.amplicon_configs import load_carlin_config_by_locus  # type: ignore

        load_carlin_config_by_locus(locus=locus)
    except ValueError as e:
        raise BulkInputError(str(e)) from e
    except FileNotFoundError as e:
        raise BulkInputError(str(e)) from e
