from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ScrnaProtocol:
    name: str
    cb_len: int
    umi_len: int
    default_whitelist_path: Path


_PROTOCOLS: dict[str, ScrnaProtocol] = {
    "10xv3": ScrnaProtocol(
        name="10xv3",
        cb_len=16,
        umi_len=12,
        default_whitelist_path=Path("reference/whitelist/10xv3.txt.gz"),
    ),
}


def get_protocol(name: str) -> ScrnaProtocol:
    key = str(name).strip()
    try:
        return _PROTOCOLS[key]
    except KeyError as exc:
        choices = ", ".join(sorted(_PROTOCOLS))
        raise ValueError(f"Unknown protocol: {name!r}. Available protocols: {choices}") from exc

