from __future__ import annotations

from pathlib import Path

import numpy as np  # type: ignore


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


def _reads_cutoff_grid(max_reads: int) -> np.ndarray:
    max_reads = int(max_reads)
    if max_reads <= 0:
        return np.array([1], dtype=int)

    segments: list[np.ndarray] = [np.arange(1, min(10, max_reads + 1), 1, dtype=int)]
    if max_reads >= 11:
        segments.append(np.arange(11, min(50, max_reads + 1), 3, dtype=int))
    coarse_stop = max_reads // 2
    if max_reads >= 61 and coarse_stop >= 61:
        segments.append(np.arange(61, coarse_stop + 1, 10, dtype=int))
    return np.unique(np.concatenate(segments))


def write_reads_cutoff_plot(df, diagnostics_dir: Path, *, reads_cutoff: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore

    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(4, 4))
    ax_top, ax_bottom = axes[0], axes[1]

    if len(df) == 0 or float(df["reads"].sum()) <= 0:
        for ax, label in ((ax_top, "No molecules"), (ax_bottom, "No molecules")):
            ax.text(0.5, 0.5, label, ha="center", va="center", transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(diagnostics_dir / "qc1_reads_cutoff.png", dpi=150)
        plt.close(fig)
        return

    max_reads = int(df["reads"].max())
    cutoff_values = _reads_cutoff_grid(max_reads)

    num_molecules = [(df[df["reads"] >= int(c)].shape[0]) for c in cutoff_values]
    total_reads = float(df["reads"].sum())
    fraction_reads_retained = [
        float(df[df["reads"] >= int(c)]["reads"].sum()) / total_reads for c in cutoff_values
    ]

    ax_top.plot(cutoff_values, num_molecules, marker="o", markersize=2, linewidth=1)
    ax_top.axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    ax_top.set_xlabel("Reads Cutoff")
    ax_top.set_ylabel("Number of Molecules")
    ax_top.set_yscale("log")
    ax_top.set_xscale("log")
    ax_top.grid(alpha=0.3)

    ax_bottom.plot(
        cutoff_values,
        fraction_reads_retained,
        marker="o",
        markersize=2,
        linewidth=1,
    )
    ax_bottom.axvline(reads_cutoff, color="red", linestyle="--", linewidth=0.6)
    ax_bottom.set_xlabel("Reads Cutoff")
    ax_bottom.set_ylabel("Frac. of Reads\nRetained")
    ax_bottom.set_ylim(0, 1.05)
    ax_bottom.set_xscale("log")
    ax_bottom.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(diagnostics_dir / "qc1_reads_cutoff.png", dpi=150)
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
    fig_path = diagnostics_dir / "qc2_pcr_chimera.png"
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
    fig.savefig(diagnostics_dir / "qc3_capture_oligo_carryover.png", dpi=150)
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
    fig.savefig(diagnostics_dir / "qc3_cells_above_k_cutoff.png", dpi=150)
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
    fig.savefig(diagnostics_dir / "qc4_lineage_barcodes_per_cell.png", dpi=150)
    plt.close(fig)
