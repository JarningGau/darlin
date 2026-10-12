from __future__ import annotations

import logging
import subprocess
from itertools import islice, zip_longest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterator, TextIO

from darlin.bulk.pear import assemble_pe_reads
from darlin.scrna.io import open_text_maybe_gzip
from darlin.scrna.validate import ScrnaInputError


def normalize_read_id(header: str) -> str:
    """Use the first header token, without its optional terminal mate suffix."""
    tokens = header.lstrip("@").split()
    if not tokens:
        raise ScrnaInputError("Empty FASTQ read ID.")
    name = tokens[0]
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    if not name:
        raise ScrnaInputError("Empty FASTQ read ID.")
    return name


def iter_fastq_records(handle: TextIO) -> Iterator[tuple[str, str, str]]:
    """Read complete FASTQ records, retaining their original header."""
    while header := handle.readline():
        sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
        if not all((sequence, plus, quality)):
            raise ScrnaInputError("Incomplete FASTQ record encountered.")
        if not header.startswith("@") or not plus.startswith("+"):
            raise ScrnaInputError("Invalid FASTQ structure.")
        sequence, quality = sequence.strip(), quality.strip()
        if len(sequence) != len(quality):
            raise ScrnaInputError("FASTQ sequence and quality lengths do not match.")
        yield header.strip(), sequence, quality


def iter_assembled_camellia_reads(
    *,
    fq1: Path,
    fq2: Path,
    pear_dir: Path,
    barcode_len: int,
    max_reads: int | None,
    pear_path: str,
    threads: int,
    logger: logging.Logger,
) -> Iterator[tuple[str, str]]:
    """Yield assembled sequence and original R2 barcode prefix, joined by ID."""
    pear_dir.mkdir(parents=True, exist_ok=True)
    prefix = pear_dir / "pear"
    assembled = pear_dir / "pear.assembled.fastq"
    pear_log = pear_dir / "pear.log"
    barcodes: dict[str, str] = {}

    # Normalize both headers for PEAR and sample input pairs before assembly.
    # Temporary input copies are deleted; all PEAR products remain available.
    with TemporaryDirectory(prefix="input-", dir=pear_dir) as temporary:
        input1, input2 = Path(temporary) / "r1.fastq", Path(temporary) / "r2.fastq"
        with (
            open_text_maybe_gzip(fq1) as h1,
            open_text_maybe_gzip(fq2) as h2,
            input1.open("w") as out1,
            input2.open("w") as out2,
        ):
            pairs = zip_longest(iter_fastq_records(h1), iter_fastq_records(h2))
            for record1, record2 in islice(pairs, max_reads):
                if record1 is None or record2 is None:
                    raise ScrnaInputError("R1 and R2 contain different numbers of reads.")
                id1, seq1, qual1 = record1
                id2, seq2, qual2 = record2
                name = normalize_read_id(id1)
                if name != normalize_read_id(id2):
                    raise ScrnaInputError(f"Paired FASTQ ID mismatch: {id1!r} versus {id2!r}.")
                if name in barcodes:
                    raise ScrnaInputError(f"Duplicate paired FASTQ ID: {name!r}.")
                if len(seq2) < barcode_len:
                    raise ScrnaInputError(f"R2 is too short for CB+UMI: {name!r}.")
                barcodes[name] = seq2[:barcode_len]
                out1.write(f"@{name}\n{seq1}\n+\n{qual1}\n")
                out2.write(f"@{name}\n{seq2}\n+\n{qual2}\n")

        logger.info("Input read pairs: %d", len(barcodes))
        # Never consume an assembly left over from an earlier invocation.
        for suffix in ("assembled.fastq", "unassembled.forward.fastq", "unassembled.reverse.fastq", "discarded.fastq"):
            (pear_dir / f"pear.{suffix}").unlink(missing_ok=True)
        if barcodes:
            try:
                assemble_pe_reads(
                    in_fq1=str(input1), in_fq2=str(input2), out_prefix=prefix,
                    log_file=pear_log, pear_path=pear_path, threads=threads, logger=logger,
                )
            except (subprocess.CalledProcessError, OSError) as exc:
                raise ScrnaInputError(f"PEAR assembly failed; see {pear_log}: {exc}") from exc
            if not assembled.is_file():
                raise ScrnaInputError(f"PEAR assembly failed: assembled FASTQ missing; see {pear_log}.")
        else:
            assembled.write_text("")
            pear_log.write_text("No input pairs; PEAR was not invoked.\n")

    seen: set[str] = set()
    with assembled.open() as handle:
        for header, sequence, _quality in iter_fastq_records(handle):
            name = normalize_read_id(header)
            if name not in barcodes:
                raise ScrnaInputError(f"Assembled read has no original pair: {name!r}.")
            if name in seen:
                raise ScrnaInputError(f"Duplicate assembled read ID: {name!r}.")
            seen.add(name)
            yield sequence, barcodes[name]
    logger.info("Assembled read pairs: %d", len(seen))
    logger.info("Unassembled read pairs: %d", len(barcodes) - len(seen))
