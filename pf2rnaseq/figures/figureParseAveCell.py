"""
Plots percentage of Tregs across cytokine conditions
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import read_h5ad

from .common import getSetup, subplotLabel


def plot_treg_percentage_by_condition(
    X,
    ax,
    condition_col: str = "cytokine",
    celltype_col: str = "cell_type",
    treg_label: str = "Treg",
    sort_by: str = "percentage",
    show_counts: bool = True,
    pbs_label: str = "PBS",
    show_pbs_line: bool = True,
    denominator: str = "all", 
    tcell_types: list = None,  
):
    """Plot percentage of Tregs for each cytokine condition

    Parameters
    ----------
    X : anndata.AnnData
        AnnData object
    ax : Axes
        Matplotlib axes
    condition_col : str
        Column name for conditions
    celltype_col : str
        Column name for cell types
    treg_label : str
        Label for Tregs
    sort_by : str
        "percentage" or "alphabetical"
    show_counts : bool
        Show cell counts on bars
    pbs_label : str
        Label for PBS condition
    show_pbs_line : bool
        Show PBS reference line
    denominator : str
        "all" = Tregs out of all cells
        "tcells" = Tregs out of T cells only
    tcell_types : list or None
        List of T cell type labels. If None, uses default:
        ["CD4 Memory", "CD4 Naive", "CD8 Memory", "CD8 Naive", "Treg"]

    Returns
    -------
    stats_df : pd.DataFrame
        Statistics for each condition
    """
    # Default T cell types
    if tcell_types is None:
        tcell_types = ["CD4 Memory", "CD4 Naive", "CD8 Memory", "CD8 Naive", "Treg"]

    # Create DataFrame for easier manipulation
    df = pd.DataFrame(
        {
            "condition": X.obs[condition_col].values,
            "cell_type": X.obs[celltype_col].values,
        }
    )

    # Calculate statistics based on denominator choice
    if denominator == "all":
        # Original: Tregs out of all cells
        total_counts = df.groupby("condition").size()
        ylabel = "Treg Percentage (% of all cells)"
        title_suffix = "(% of All Cells)"
    elif denominator == "tcells":
        # NEW: Tregs out of T cells only
        tcell_df = df[df["cell_type"].isin(tcell_types)]
        total_counts = tcell_df.groupby("condition").size()
        ylabel = "Treg Percentage (% of T cells)"
        title_suffix = "(% of T Cells)"
    else:
        raise ValueError(f"denominator must be 'all' or 'tcells', got '{denominator}'")

    # Count Tregs per condition
    treg_counts = df[df["cell_type"] == treg_label].groupby("condition").size()

    # Create stats dataframe
    stats_df = pd.DataFrame(
        {
            "Condition": total_counts.index,
            "N_Total": total_counts.values,
            "N_Treg": treg_counts.reindex(total_counts.index, fill_value=0).values,
        }
    )

    # Calculate percentage
    stats_df["Pct_Treg"] = stats_df["N_Treg"] / stats_df["N_Total"] * 100

    # Get PBS percentage before sorting
    pbs_pct = None
    if show_pbs_line and pbs_label in stats_df["Condition"].values:
        pbs_pct = stats_df[stats_df["Condition"] == pbs_label]["Pct_Treg"].values[0]

    # Sort
    if sort_by == "percentage":
        stats_df = stats_df.sort_values("Pct_Treg", ascending=False)
    else:
        stats_df = stats_df.sort_values("Condition")

    # Create bar plot
    x_pos = np.arange(len(stats_df))
    colors = plt.cm.RdYlBu_r(stats_df["Pct_Treg"].values / stats_df["Pct_Treg"].max())

    bars = ax.bar(
        x_pos,
        stats_df["Pct_Treg"].values,
        color=colors,
        edgecolor="black",
        linewidth=1.5,
    )

    # Add PBS reference line
    if show_pbs_line and pbs_pct is not None:
        ax.axhline(
            y=pbs_pct,
            color="red",
            linestyle="--",
            linewidth=2.5,
            alpha=0.8,
            label=f"PBS baseline ({pbs_pct:.1f}%)",
            zorder=10,
        )
        ax.legend(loc="upper right", fontsize=10, framealpha=0.9)

    # Add value labels on bars
    for i, (bar, row) in enumerate(zip(bars, stats_df.itertuples(), strict=False)):
        height = bar.get_height()

        # Only show label for top 5 bars
        if i < 5:
            label = (
                f"{row.Pct_Treg:.1f}%\n({row.N_Treg}/{row.N_Total})"
                if show_counts
                else f"{row.Pct_Treg:.1f}%"
            )
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                label,
                ha="center",
                va="bottom",
                fontsize=7,
                fontweight="bold",
                rotation=45,
            )

    # Formatting
    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        stats_df["Condition"].values, rotation=45, ha="right", fontsize=10
    )
    ax.set_xlabel("Cytokine Condition", fontsize=12, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=12, fontweight="bold")
    ax.set_title(
        f"Percentage of Tregs by Cytokine Condition {title_suffix}",
        fontsize=14,
        fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    # Extend y-axis to make room for labels
    y_max = stats_df["Pct_Treg"].max()
    ax.set_ylim(bottom=0, top=y_max * 1.15)

    # Update summary text
    summary_text = f"PBS: {pbs_pct:.1f}%\n"
    summary_text += (
        f"Denominator: {'All cells' if denominator == 'all' else 'T cells only'}"
    )

    ax.text(
        0.98,
        0.98,
        summary_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.7),
        fontsize=9,
    )

    return stats_df


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Create 2 subplots to compare both methods
    ax, f = getSetup((30, 10), (1, 1))

    subplotLabel(ax)

    X = read_h5ad("/opt/data/Parse_10M_PBMC_cytokines.h5ad", backed="r")

    # Panel A: Tregs as % of all cells

    # Panel B: Tregs as % of T cells
    stats_df_tcells = plot_treg_percentage_by_condition(
        X,
        ax[0],
        condition_col="cytokine",
        celltype_col="cell_type",
        treg_label="Treg",
        sort_by="percentage",
        show_counts=True,
        pbs_label="PBS",
        show_pbs_line=True,
        denominator="tcells",  # Tregs out of T CELLS only
        tcell_types=["CD4 Memory", "CD4 Naive", "CD8 Memory", "CD8 Naive", "Treg"],
    )

    print("\n=== Treg Percentage (% of T Cells) ===")
    print(stats_df_tcells.to_string(index=False))

    return f
