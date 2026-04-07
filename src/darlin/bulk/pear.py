from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path


def assemble_pe_reads(
    in_fq1: str,
    in_fq2: str,
    out_prefix: str | Path,
    log_file: str | Path,
    pear_path: str = "pear",
    threads: int = 8,
    logger: logging.Logger | None = None,
) -> None:
    if logger is None:
        logger = logging.getLogger(__name__)

    pear_exec = os.path.expanduser(pear_path)
    out_prefix = str(out_prefix)
    cmd = [pear_exec, "-f", in_fq1, "-r", in_fq2, "-o", out_prefix, "-j", str(threads)]

    logger.info(f"Running PEAR assembly: {' '.join(cmd)}")
    with open(log_file, "w") as log_handle:
        subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=True)
    logger.info(f"PEAR assembly completed. Log saved to: {log_file}")

