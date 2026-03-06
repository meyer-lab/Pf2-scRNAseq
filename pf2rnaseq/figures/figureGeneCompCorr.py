"""Correlation of gene expression with component weights across cytokines"""

import anndata
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from scipy.stats import pearsonr

from ..factorization import correct_conditions
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import (
    cell_count_perc_df,
)


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((10, 10), (1, 1))
    subplotLabel(ax)
    X = anndata.read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X.uns["Pf2_A"] = correct_conditions(X)

    cellDF = cell_count_perc_df(X, "cell_type")

    # Highlight specific cytokines
    # In your makeFigure function
    highlight_cytokines = ["IL-2", "IL-3", "CT-1", "IL-1-beta", "IFN-epsilon"]
    plot_gene_comp_corr(
        X, "FOXP3", "Treg", 38, ax[0], highlight_cytokines=highlight_cytokines
    )

    return f


def plot_gene_comp_corr(
    X: anndata.AnnData,
    gene: str,
    pop: str,
    comp: int,
    ax: Axes,
    unique=None,
    highlight_cytokines=None,
):
    """Plots correlation of gene expression against component weights for each cytokine"""

    # Filter for specific cell type
    cell_mask = X.obs["cell_type"] == pop
    X_filtered = X[cell_mask, :]

    # Calculate average gene expression per cytokine for this cell type
    gene_expr_df = pd.DataFrame()
    cytokines = np.unique(X_filtered.obs["cytokine"])

    for cyt in cytokines:
        cyt_mask = X_filtered.obs["cytokine"] == cyt
        if np.sum(cyt_mask) > 0:  # Check if there are cells for this cytokine
            avg_expr = np.mean(X_filtered[cyt_mask, gene].X)
            gene_expr_df = pd.concat(
                [
                    gene_expr_df,
                    pd.DataFrame({"cytokine": [cyt], f"{gene}_avg": [avg_expr]}),
                ]
            )

    # Get component weights
    comp_df = pd.DataFrame(
        {
            "Comp. " + str(comp): X.uns["Pf2_A"][:, comp - 1],
            "cytokine": np.unique(X.obs["cytokine"]),
        }
    )

    # Merge the dataframes
    newDF = gene_expr_df.merge(comp_df, on="cytokine")

    if unique is not None:
        newDF["cytokine"] = newDF["cytokine"].astype(str)
        newDF.loc[~newDF["cytokine"].isin(unique), "cytokine"] = "Other"

    # Calculate Pearson correlation
    comp_col = "Comp. " + str(comp)
    gene_col = f"{gene}_avg"
    corr_coef, p_value = pearsonr(newDF[comp_col], newDF[gene_col])

    # Create highlight column
    if highlight_cytokines is not None:
        newDF["highlight"] = newDF["cytokine"].isin(highlight_cytokines)

        # Plot non-highlighted points first (in background)
        non_highlight = newDF[~newDF["highlight"]]
        sns.scatterplot(
            non_highlight,
            x=comp_col,
            y=gene_col,
            color="lightgray",
            alpha=0.6,
            ax=ax,
            legend=False,
        )

        # Plot highlighted points on top
        highlight = newDF[newDF["highlight"]]
        sns.scatterplot(
            highlight,
            x=comp_col,
            y=gene_col,
            hue="cytokine",
            s=100,
            ax=ax,
            legend=False,
        )

        # Add labels next to highlighted points
        for _, row in highlight.iterrows():
            ax.text(
                row[comp_col] + 0.02,
                row[gene_col],
                row["cytokine"],
                fontsize=14,
                fontweight="bold",
                va="center",
            )
    else:
        sns.scatterplot(newDF, x=comp_col, y=gene_col, hue="cytokine", ax=ax)

    # Add correlation statistics to the plot
    ax.text(
        0.05,
        0.95,
        f"r = {corr_coef:.3f}\np = {p_value:.3e}",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        verticalalignment="top",
    )

    # Increase font sizes
    ax.set_xlabel(f"Comp. {comp} weights", fontsize=16, fontweight="bold")
    ax.set_ylabel(f"{gene} avg expression ({pop})", fontsize=16, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=14)

    return newDF
