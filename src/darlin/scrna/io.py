from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator, TextIO


def open_text_maybe_gzip(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def iter_fastq_paired(handle1: TextIO, handle2: TextIO) -> Iterator[tuple[str, str, str, str, str, str]]:
    while True:
        id_line1 = handle1.readline()
        if not id_line1:
            return
        seq_line1 = handle1.readline()
        plus_line1 = handle1.readline()
        qual_line1 = handle1.readline()

        id_line2 = handle2.readline()
        seq_line2 = handle2.readline()
        plus_line2 = handle2.readline()
        qual_line2 = handle2.readline()

        if not all([seq_line1, plus_line1, qual_line1, id_line2, seq_line2, plus_line2, qual_line2]):
            raise ValueError("Incomplete paired FASTQ record encountered.")
        if not id_line1.startswith("@") or not plus_line1.startswith("+"):
            raise ValueError("Invalid FASTQ structure in R1.")
        if not id_line2.startswith("@") or not plus_line2.startswith("+"):
            raise ValueError("Invalid FASTQ structure in R2.")

        seq1 = seq_line1.strip()
        qual1 = qual_line1.strip()
        seq2 = seq_line2.strip()
        qual2 = qual_line2.strip()
        if len(seq1) != len(qual1):
            raise ValueError("R1 sequence and quality lengths do not match.")
        if len(seq2) != len(qual2):
            raise ValueError("R2 sequence and quality lengths do not match.")

        yield id_line1[1:].strip(), seq1, qual1, id_line2[1:].strip(), seq2, qual2


def load_whitelist(path: str | Path) -> set[str]:
    with open_text_maybe_gzip(path) as handle:
        return {line.strip() for line in handle if line.strip()}
