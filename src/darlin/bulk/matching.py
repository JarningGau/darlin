from __future__ import annotations

from dataclasses import dataclass


def get_mm_dist(seq: str, rate: float = 0.05, n: int = 2) -> int:
    return max(int(len(seq) * rate), n)


def find_exact_matches(seq_str: str, patterns: dict[str, str]) -> dict[str, list[int]]:
    matches: dict[str, list[int]] = {}
    for pattern_name, pattern_seq in patterns.items():
        positions: list[int] = []
        start = 0
        while True:
            pos = seq_str.find(pattern_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1  # allow overlapping
        matches[pattern_name] = positions
    return matches


def find_fuzzy_matches(seq_str: str, pattern: str, max_errors: int):
    # Lazy import so `--help` works without optional deps.
    from fuzzysearch import find_near_matches  # type: ignore

    return find_near_matches(pattern, seq_str, max_l_dist=max_errors)


@dataclass
class MatchResult:
    p3_matches: list
    p5_matches: list
    match_method: str


def find_all_matches(seq_str: str, p3_seq: str, p5_seq: str, p3_mm: int, p5_mm: int) -> MatchResult:
    patterns = {"p3": p3_seq, "p5": p5_seq}
    exact = find_exact_matches(seq_str, patterns)

    exact_success = len(exact["p3"]) == 1 and len(exact["p5"]) == 1
    if exact_success:
        p3 = [type("Match", (), {"start": exact["p3"][0], "end": exact["p3"][0] + len(p3_seq)})()]
        p5 = [type("Match", (), {"start": exact["p5"][0], "end": exact["p5"][0] + len(p5_seq)})()]
        return MatchResult(p3_matches=p3, p5_matches=p5, match_method="exact")

    mixed_success = True
    methods_used: list[str] = []

    # p3
    if len(exact["p3"]) == 1:
        p3_matches = [type("Match", (), {"start": exact["p3"][0], "end": exact["p3"][0] + len(p3_seq)})()]
        methods_used.append("exact")
    else:
        p3_matches = find_fuzzy_matches(seq_str, p3_seq, p3_mm)
        methods_used.append("fuzzy")
        if len(p3_matches) != 1:
            mixed_success = False

    # p5
    if mixed_success:
        if len(exact["p5"]) == 1:
            p5_matches = [type("Match", (), {"start": exact["p5"][0], "end": exact["p5"][0] + len(p5_seq)})()]
            methods_used.append("exact")
        else:
            p5_matches = find_fuzzy_matches(seq_str, p5_seq, p5_mm)
            methods_used.append("fuzzy")
            if len(p5_matches) != 1:
                mixed_success = False
    else:
        p5_matches = []

    if mixed_success:
        method = "mixed" if "fuzzy" in methods_used else "exact"
        return MatchResult(p3_matches=p3_matches, p5_matches=p5_matches, match_method=method)

    return MatchResult(p3_matches=[], p5_matches=[], match_method="failed")

