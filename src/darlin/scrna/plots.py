from __future__ import annotations

from pathlib import Path


def _k_category_reads_per_umi(k: float) -> str:
    """Bin k = reads / UMIs for scatter coloring (matches notebook darlin-scrna)."""
    if k <= 1:
        return "≤1"
    if k <= 5:
        return "≤5"
    if k <= 10:
        return "≤10"
    return ">10"


def write_extract_plots(df, diagnostics_dir: Path, *, unedited_bc_len: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore

    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4, 2))
    if len(df) > 0:
        ax.hist(df["LB_len"], bins=range(1, 300, 1), edgecolor="black")
    else:
        ax.text(0.5, 0.5, "No matched reads", ha="center", va="center", transform=ax.transAxes)
    ax.axvline(unedited_bc_len, color="red", linestyle="--", linewidth=0.6)
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Number of reads")
    # ax.set_title("Distribution of DARLIN Array Sequence\nLengths By Reads")
    fig.tight_layout()
    fig.savefig(diagnostics_dir / "fragment_length_distribution.png", dpi=150)
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
    import matplotlib.patches as mpatches  # type: ignore
    import matplotlib.ticker as mtick  # type: ignore
    from matplotlib.ticker import ScalarFormatter  # type: ignore

    diagnostics_dir.mkdir(parents=True, exist_ok=True)

    # QC 2: reads-fraction (matches notebook — one 2×1 figure, same style)
    fig, axes = plt.subplots(2, 1, figsize=(3, 4))
    ax0, ax1 = axes[0], axes[1]
    ax0.hist(df_all["reads_fraction"], bins=50, edgecolor="white")
    ax0.axvline(major_fraction_threshold_molecule, color="red", linestyle="--", linewidth=0.6)
    ax0.set_ylabel("Number of (CR, UR)")
    ax0.yaxis.set_major_formatter(mtick.ScalarFormatter(useMathText=True))
    ax0.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    ax1.scatter(df_all["reads_fraction"], df_all["reads"], s=0.1, alpha=0.1)
    ax1.axvline(major_fraction_threshold_molecule, color="red", linestyle="--", linewidth=0.6)
    ax1.set_xlabel("Reads Fraction")
    ax1.set_ylabel("Reads")
    ax1.set_yscale("log")

    fig.tight_layout()
    fig_path = diagnostics_dir / "qc1_pcr_chimera.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 3))
    if len(cell_summary) > 0:
        k_cat_order = ["≤1", "≤5", "≤10", ">10"]
        k_cat_colors = {
            "≤1": "#4575b4",
            "≤5": "#91bfdb",
            "≤10": "#fee090",
            ">10": "#d73027",
        }
        k_cats = cell_summary["k"].map(_k_category_reads_per_umi)
        for cat in k_cat_order:
            mask = k_cats == cat
            sub = cell_summary.loc[mask]
            if len(sub) == 0:
                continue
            ax.scatter(
                sub["n_reads"],
                sub["n_UR"],
                s=2,
                alpha=0.4,
                color=k_cat_colors[cat],
                label=cat,
            )
        ur_min = float(cell_summary["n_UR"].min())
        ur_max = float(cell_summary["n_UR"].max())
        ax.plot(
            [ur_min, ur_max],
            [ur_min, ur_max],
            linestyle="--",
            color="red",
            linewidth=1,
            label="slope=1",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        legend_handles = [
            mpatches.Patch(color=k_cat_colors[cat], label=f"k {cat}") for cat in k_cat_order
        ]
        ax.legend(handles=legend_handles, title="k = Reads/UMIs", loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_xlabel("Reads per cell")
    ax.set_ylabel("UMIs per cell")
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    fig.savefig(diagnostics_dir / "qc2_capture_oligo_carryover.png", dpi=150)
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
    fig.savefig(diagnostics_dir / "qc2_cells_above_k_cutoff.png", dpi=150)
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
    fig.savefig(diagnostics_dir / "qc3_lineage_barcodes_per_cell.png", dpi=150)
    plt.close(fig)
