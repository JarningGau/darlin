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


def iter_fastq_paired(
    handle1: TextIO, handle2: TextIO
) -> Iterator[tuple[str, str, str, str, str, str]]:
    """
    Paired FASTQ iterator: read synchronously from two FASTQ streams.

    Yields (read_id1, seq1, qual1, read_id2, seq2, qual2) per pair.
    Stops when file 1 ends; raises if file 2 ends first or records are malformed.
    """
    while True:
        id_line1 = handle1.readline()
        if not id_line1:
            break
        seq_line1 = handle1.readline()
        plus_line1 = handle1.readline()
        qual_line1 = handle1.readline()

        if not (seq_line1 and plus_line1 and qual_line1):
            raise ValueError("Incomplete FASTQ record encountered in file 1.")

        if not id_line1.startswith("@") or not plus_line1.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line) in file 1.")

        id_line2 = handle2.readline()
        if not id_line2:
            raise ValueError("File 2 ended before file 1.")
        seq_line2 = handle2.readline()
        plus_line2 = handle2.readline()
        qual_line2 = handle2.readline()

        if not (seq_line2 and plus_line2 and qual_line2):
            raise ValueError("Incomplete FASTQ record encountered in file 2.")

        if not id_line2.startswith("@") or not plus_line2.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line) in file 2.")

        read_id1 = id_line1[1:].strip()
        seq1 = seq_line1.strip()
        qual1 = qual_line1.strip()
        if len(seq1) != len(qual1):
            raise ValueError(
                f"Length mismatch (seq {len(seq1)} vs qual {len(qual1)}) at read {read_id1} in file 1"
            )

        read_id2 = id_line2[1:].strip()
        seq2 = seq_line2.strip()
        qual2 = qual_line2.strip()
        if len(seq2) != len(qual2):
            raise ValueError(
                f"Length mismatch (seq {len(seq2)} vs qual {len(qual2)}) at read {read_id2} in file 2"
            )

        yield read_id1, seq1, qual1, read_id2, seq2, qual2

