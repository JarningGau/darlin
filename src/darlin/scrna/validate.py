from __future__ import annotations

from pathlib import Path

from darlin.scrna.protocols import ScrnaProtocol, get_protocol


class ScrnaInputError(ValueError):
    pass


def require_readable_files(items: list[tuple[str, str | Path]]) -> None:
    missing: list[str] = []
    for label, raw_path in items:
        path = Path(raw_path)
        if not path.exists():
            missing.append(f"{label}: file not found: {path}")
        elif not path.is_file():
            missing.append(f"{label}: not a file: {path}")
    if missing:
        raise ScrnaInputError("\n".join(missing))


def require_protocol(name: str) -> ScrnaProtocol:
    try:
        return get_protocol(name)
    except ValueError as exc:
        raise ScrnaInputError(str(exc)) from exc


def resolve_whitelist_path(protocol: ScrnaProtocol, override: str | None) -> Path:
    path = Path(override) if override else protocol.default_whitelist_path
    require_readable_files([("Whitelist (--whitelist or protocol default)", path)])
    return path

