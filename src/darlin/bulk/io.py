from __future__ import annotations

import gzip
from typing import Iterator, TextIO


def open_fastq_text(file_path: str) -> TextIO:
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def iter_fastq_raw(handle: TextIO) -> Iterator[tuple[str, str, str]]:
    while True:
        id_line = handle.readline()
        if not id_line:
            break
        seq_line = handle.readline()
        plus_line = handle.readline()
        qual_line = handle.readline()

        if not (seq_line and plus_line and qual_line):
            raise ValueError("Incomplete FASTQ record encountered.")
        if not id_line.startswith("@") or not plus_line.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line).")

        read_id = id_line[1:].strip()
        seq = seq_line.strip()
        qual = qual_line.strip()
        if len(seq) != len(qual):
            raise ValueError(
                f"Length mismatch (seq {len(seq)} vs qual {len(qual)}) at read {read_id}"
            )
        yield read_id, seq, qual

