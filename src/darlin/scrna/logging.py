from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd  # type: ignore


def log_step_title(logger: logging.Logger, title: str) -> None:
    logger.info('-' * 40)
    logger.info(title)
    logger.info('-' * 40)


def log_count(logger: logging.Logger, label: str, value: int) -> None:
    logger.info("%s: %s", label, f"{int(value):,}")


def log_count_with_pct(
    logger: logging.Logger,
    label: str,
    value: int,
    total: int,
) -> None:
    count_str = f"{int(value):,}"
    if total > 0:
        logger.info("%s: %s (%s)", label, count_str, f"{value / total:.2%}")
    else:
        logger.info("%s: %s", label, count_str)


def _barcode_nunique(df: pd.DataFrame, preferred: str, fallback: str) -> int:
    col = preferred if preferred in df.columns else fallback if fallback in df.columns else None
    if col is None or len(df) == 0:
        return 0
    return int(df[col].nunique())


def log_molecule_summary(
    logger: logging.Logger,
    df: pd.DataFrame,
    *,
    read_col: str = "reads",
) -> None:
    if len(df) == 0:
        n_reads = 0
        n_cr = n_ur = n_lr = 0
    else:
        n_reads = int(df[read_col].sum()) if read_col in df.columns else len(df)
        n_cr = _barcode_nunique(df, "CR", "CB")
        n_ur = _barcode_nunique(df, "UR", "UB")
        n_lr = _barcode_nunique(df, "LR", "LB")
    log_count(logger, "Number of reads", n_reads)
    log_count(logger, "Number of molecules", len(df))
    log_count(logger, "Number of UMIs", n_ur)
    log_count(logger, "Number of cell barcodes", n_cr)
    log_count(logger, "Number of lineage barcodes", n_lr)


def log_output(logger: logging.Logger, path: Path) -> None:
    logger.info("-> %s", path)


def setup_logging(log_file: Path, level: int) -> logging.Logger:
    logger = logging.getLogger(f"darlin.scrna.{log_file}")
    logger.setLevel(level)
    logger.propagate = False

    if logger.handlers:
        return logger

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger

