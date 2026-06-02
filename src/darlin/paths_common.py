from __future__ import annotations

import re

_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def validate_sample_id(sample_id: str) -> None:
    if not sample_id or not sample_id.strip():
        raise ValueError("sample_id must be a non-empty identifier")
    if sample_id in (".", ".."):
        raise ValueError("sample_id cannot be '.' or '..'")
    if not _SAMPLE_ID_RE.fullmatch(sample_id):
        raise ValueError(
            "sample_id must contain only letters, digits, '.', '_', or '-' "
            f"(got {sample_id!r})"
        )
