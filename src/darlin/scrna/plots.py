from __future__ import annotations

from pathlib import Path


def write_extract_plots(df, diagnostics_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore

    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4, 2.5))
    if len(df) > 0:
        lengths = df["LB_len"].astype(int)
        bins = range(1, max(int(lengths.max()) + 2, 3))
        ax.hist(lengths, bins=bins, edgecolor="white")
    else:
        ax.text(0.5, 0.5, "No matched reads", ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel("LB length")
    ax.set_ylabel("Reads")
    ax.set_title("Extracted lineage barcode lengths")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "extract_lb_length.png", dpi=150)
    plt.close(fig)


def write_qc_plots(
    df_all,
    cell_summary,
    df_final,
    diagnostics_dir: Path,
    *,
    major_fraction_threshold_molecule: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.ticker import ScalarFormatter  # type: ignore

    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(4, 2.5))
    ax.hist(df_all["reads_fraction"], bins=50, edgecolor="white")
    ax.axvline(major_fraction_threshold_molecule, color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Reads fraction")
    ax.set_ylabel("Count")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "qc_reads_fraction_hist.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4, 2.5))
    ax.scatter(df_all["reads_fraction"], df_all["reads"], s=4, alpha=0.2)
    ax.axvline(major_fraction_threshold_molecule, color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Reads fraction")
    ax.set_ylabel("Reads")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "qc_reads_fraction_scatter.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4, 3))
    if len(cell_summary) > 0:
        ax.scatter(cell_summary["n_reads"], cell_summary["n_UR"], s=6, alpha=0.4)
        ax.set_xscale("log")
        ax.set_yscale("log")
    ax.set_xlabel("Reads")
    ax.set_ylabel("UMIs")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "qc_reads_vs_umis.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4, 2.5))
    if len(cell_summary) > 0:
        k_cutoffs = list(range(1, 20))
        n_cr_above_k = [(cell_summary["k"] >= cutoff).sum() for cutoff in k_cutoffs]
        ax.plot(k_cutoffs, n_cr_above_k, marker="o")
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_xlabel("k cutoff")
    ax.set_ylabel("Cells >= cutoff")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "qc_k_cutoff_curve.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4, 2.5))
    if len(df_final) > 0:
        df_plot = df_final[["CR", "n_LR"]].drop_duplicates()
        bins = range(1, max(int(df_plot["n_LR"].max()) + 2, 3))
        ax.hist(df_plot["n_LR"], bins=bins, edgecolor="white")
        ax.set_yscale("log")
    ax.set_xlabel("LRs per cell")
    ax.set_ylabel("Cells")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "qc_n_lr_per_cr_hist.png", dpi=150)
    plt.close(fig)
