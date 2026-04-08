from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path


def _display_path(path: str, bases: tuple[Path, ...]) -> str:
    """Prefer paths relative to bases (first match wins) for shorter log lines."""
    try:
        p = Path(path).expanduser().resolve()
    except OSError:
        return path
    for base in bases:
        try:
            b = base.expanduser().resolve()
            return str(p.relative_to(b))
        except ValueError:
            continue
    return str(p)


def assemble_pe_reads(
    in_fq1: str,
    in_fq2: str,
    out_prefix: str | Path,
    log_file: str | Path,
    pear_path: str = "pear",
    threads: int = 8,
    logger: logging.Logger | None = None,
    log_path_bases: tuple[Path, ...] = (),
) -> None:
    if logger is None:
        logger = logging.getLogger(__name__)

    pear_exec = os.path.expanduser(pear_path)
    out_prefix = str(out_prefix)
    cmd = [pear_exec, "-f", in_fq1, "-r", in_fq2, "-o", out_prefix, "-j", str(threads)]
    if log_path_bases:
        f1 = _display_path(in_fq1, log_path_bases)
        f2 = _display_path(in_fq2, log_path_bases)
        op = _display_path(out_prefix, log_path_bases)
        log_cmd = " ".join([pear_exec, "-f", f1, "-r", f2, "-o", op, "-j", str(threads)])
    else:
        log_cmd = " ".join(cmd)

    logger.info(f"Running PEAR assembly: {log_cmd}")
    with open(log_file, "w") as log_handle:
        subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=True)
    log_disp = _display_path(str(log_file), log_path_bases) if log_path_bases else str(log_file)
    logger.info(f"PEAR assembly completed. Log saved to: {log_disp}")

