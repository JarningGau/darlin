from __future__ import annotations

import argparse


def positive_int(value: str) -> int:
    """Argparse type: integer strictly greater than zero."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"must be an integer, got {value!r}") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"must be a positive integer, got {parsed}")
    return parsed


def positive_float(value: str) -> float:
    """Argparse type: float strictly greater than zero."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"must be a number, got {value!r}") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"must be a positive number, got {parsed}")
    return parsed
