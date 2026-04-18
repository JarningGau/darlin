"""Inline diagnostic figures for the bulk pipeline (matplotlib Agg)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd


def _configure_matplotlib() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401

    _ = plt


def write_pre_denoise_plots(
    *,
    filtered_tsv: Path,
    out_dir: Path,
    reads_cutoff: int,
    unedited_bc_len: int,
    logger: logging.Logger | None = None,
) -> None:
    """Plots from aggregated filtered pairs (before denoise); matches legacy bulk script."""
    _configure_matplotlib()
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(filtered_tsv, sep="\t")
    if df.empty:
        if logger:
            logger.warning("Pre-denoise plots skipped: empty filtered table %s", filtered_tsv)
        return

    if "LB_len" not in df.columns and "LB" in df.columns:
        df = df.copy()
        df["LB_len"] = df["LB"].astype(str).str.len()

    _plot_barcode_length_distribution_by_pairs(df, out_dir, unedited_bc_len=unedited_bc_len)
    _plot_reads_cutoff_distribution(df, reads_cutoff, out_dir)
    _plot_cutoff_vs_num_umis(df, reads_cutoff, out_dir)
    _plot_cutoff_vs_fraction_reads_retained(df, reads_cutoff, out_dir)
    plt.close("all")
    if logger:
        logger.info("Pre-denoise diagnostic plots -> %s", out_dir)


def write_post_denoise_plots(
    *,
    denoised_agg_tsv: Path,
    out_dir: Path,
    unedited_bc_len: int,
    logger: logging.Logger | None = None,
) -> None:
    """Barcode length by UMI after denoising (from denoised_agg.tsv)."""
    _configure_matplotlib()
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(denoised_agg_tsv, sep="\t")
    if df.empty or "LR" not in df.columns:
        if logger:
            logger.warning("Post-denoise length plot skipped: missing data in %s", denoised_agg_tsv)
        plt.close("all")
        return

    agg = df.copy()
    agg["bc_len"] = agg["LR"].astype(str).str.len()

    plt.figure(figsize=(4, 2))
    plt.hist(agg["bc_len"], bins=range(1, 300, 1), edgecolor="black")
    plt.axvline(unedited_bc_len, color="red", linestyle="--", linewidth=0.6)
    plt.xlabel("Sequence Length")
    plt.ylabel("Number of UMIs")
    plt.title("Distribution of DARLIN Array Sequence\nLengths By UMI (After denoising)")
    plt.tight_layout()
    plt.savefig(out_dir / "barcode_length_distribution_by_UMI_after_denoising.png", dpi=150)
    plt.close()
    if logger:
        logger.info("Post-denoise diagnostic plot -> %s", out_dir)


def write_post_finalize_plots(
    *,
    alleles_tsv: Path,
    out_dir: Path,
    unedited_bc_len: int,
    logger: logging.Logger | None = None,
) -> None:
    """Clone size and length-by-editing-event style plots from alleles_by_umis.tsv."""
    _configure_matplotlib()
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(alleles_tsv, sep="\t")
    if df.empty:
        if logger:
            logger.warning("Post-finalize plots skipped: empty %s", alleles_tsv)
        plt.close("all")
        return

    if "UMIs" in df.columns:
        sub = df[df["UMIs"] < 500]
        plt.figure(figsize=(5, 2))
        plt.hist(sub["UMIs"], bins=range(0, 501, 1), color="steelblue", edgecolor=None)
        plt.yscale("log")
        plt.xlabel("UMIs (clone size)")
        plt.ylabel("Alleles (Clones)")
        plt.title("Clone Size Distribution")
        plt.tight_layout()
        plt.savefig(out_dir / "clone_size_distribution.png", dpi=150)
        plt.close()

    if "aligned_query" in df.columns:
        bc_lens = df["aligned_query"].astype(str).str.replace("-", "", regex=False).str.len()
        plt.figure(figsize=(4, 2))
        plt.hist(bc_lens, bins=range(1, 300, 1), edgecolor="black")
        plt.axvline(unedited_bc_len, color="red", linestyle="--", linewidth=0.6)
        plt.xlabel("Sequence Length")
        plt.ylabel("Number of Alleles")
        plt.title("Distribution of DARLIN Array Sequence\nLengths By Editing Events")
        plt.tight_layout()
        plt.savefig(out_dir / "barcode_length_distribution_by_editing_events.png", dpi=150)
        plt.close()

    plt.close("all")
    if logger:
        logger.info("Post-finalize diagnostic plots -> %s", out_dir)


def _plot_barcode_length_distribution_by_pairs(
    results_df: pd.DataFrame,
    output_dir: Path,
    *,
    unedited_bc_len: int,
) -> None:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(4, 2))
    plt.hist(results_df["LB_len"], bins=range(1, 300, 1), edgecolor="black")
    plt.axvline(unedited_bc_len, color="red", linestyle="--", linewidth=0.6)
    plt.xlabel("Sequence Length")
    plt.ylabel("Number of reads")
    plt.title("Distribution of DARLIN Array Sequence\nLengths By Reads")
    plt.tight_layout()
    plt.savefig(output_dir / "barcode_length_distribution_by_reads.png", dpi=150)
    plt.close()


def _plot_reads_cutoff_distribution(results_df: pd.DataFrame, reads_cutoff: int, output_dir: Path) -> None:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(4, 2))
    plt.hist(results_df["reads"], bins=20, edgecolor=None, color="skyblue")
    plt.xlabel("Number of Reads")
    plt.ylabel("Frequency")
    plt.axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    plt.title("Distribution of Read Counts")
    plt.grid(True, alpha=0.3)
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(output_dir / "reads_counts_distribution.png", dpi=150)
    plt.close()


def _cutoff_grid(max_reads: int) -> np.ndarray:
    max_reads = int(max_reads)
    if max_reads <= 0:
        return np.array([1], dtype=int)
    if max_reads <= 10000:
        return np.arange(1, max_reads + 1, dtype=int)
    return np.unique(np.linspace(1, max_reads, num=min(2000, max_reads), dtype=int))


def _plot_cutoff_vs_num_umis(results_df: pd.DataFrame, reads_cutoff: int, output_dir: Path) -> None:
    import matplotlib.pyplot as plt

    max_r = int(results_df["reads"].max())
    cutoff_values = _cutoff_grid(max_r)
    num_umis = [(results_df[results_df["reads"] >= int(c)]["UB"].nunique()) for c in cutoff_values]

    plt.figure(figsize=(4, 2))
    plt.plot(cutoff_values, num_umis, marker="o", markersize=2, linewidth=1)
    plt.axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    plt.xlabel("Reads Cutoff")
    plt.ylabel("Number of UMIs")
    plt.yscale("log")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "retained_UMIs_by_reads_cutoff.png", dpi=150)
    plt.close()


def _plot_cutoff_vs_fraction_reads_retained(results_df: pd.DataFrame, reads_cutoff: int, output_dir: Path) -> None:
    import matplotlib.pyplot as plt

    total = float(results_df["reads"].sum())
    if total <= 0:
        return
    max_r = int(results_df["reads"].max())
    cutoff_values = _cutoff_grid(max_r)
    fraction_reads_retained = [
        float(results_df[results_df["reads"] >= int(c)]["reads"].sum()) / total for c in cutoff_values
    ]

    plt.figure(figsize=(4, 2))
    plt.plot(cutoff_values, fraction_reads_retained, marker="o", markersize=2, linewidth=1)
    plt.axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    plt.xlabel("Reads Cutoff")
    plt.ylabel("Frac. of Reads\nRetained")
    plt.ylim(0, 1.05)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "fraction_retained_reads_by_reads_cutoff.png", dpi=150)
    plt.close()
