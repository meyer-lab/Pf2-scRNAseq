"""Figure showing correlation of cell percentages against each conditions component value"""

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
    ax, f = getSetup((5, 5), (1, 1))
    subplotLabel(ax)
    X = anndata.read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X.uns["Pf2_A"] = correct_conditions(X)

    cellDF = cell_count_perc_df(X, "cell_type")

    # Highlight specific cytokines
    highlight_cytokines = [
        "IL-2",
        "IL-15",
        "IL-7",
        "IL-1-beta",
    ]  

    plot_cell_perc_comp_corr(
        X, cellDF, "Treg", 56, ax[0], highlight_cytokines=highlight_cytokines
    )

    return f


def plot_cell_perc_comp_corr(
    X: anndata.AnnData,
    cellDF: pd.DataFrame,
    pop: str,
    comp: int,
    ax: Axes,
    unique=None,
    highlight_cytokines=None,
):
    """Plots correlation of cell percentages against each conditions component value"""
    newDF = pd.DataFrame()
    newDF[[pop, "cytokine"]] = cellDF.loc[cellDF["Cell Type"] == pop][
        ["Cell Type Percentage", "cytokine"]
    ]
    newDF2 = pd.DataFrame(
        {
            "Comp. " + str(comp): X.uns["Pf2_A"][:, comp - 1],
            "cytokine": np.unique(X.obs["cytokine"]),
        }
    )
    newDF = newDF.merge(newDF2, on="cytokine")

    if unique is not None:
        newDF["cytokine"] = newDF["cytokine"].astype(str)
        newDF.loc[~newDF["cytokine"].isin(unique), "cytokine"] = "Other"

    # Calculate Pearson correlation
    x_values = newDF["Comp. " + str(comp)]
    y_values = newDF[pop]
    pearson_coeff, p_value = pearsonr(x_values, y_values)

    # Create highlight column
    if highlight_cytokines is not None:
        newDF["highlight"] = newDF["cytokine"].isin(highlight_cytokines)

        # Plot non-highlighted points first (in background)
        non_highlight = newDF[~newDF["highlight"]]
        sns.scatterplot(
            non_highlight,
            x="Comp. " + str(comp),
            y=pop,
            color="lightgray",
            alpha=0.6,
            ax=ax,
            legend=False,
            s=150,
        )

        # Plot highlighted points on top
        highlight = newDF[newDF["highlight"]]
        sns.scatterplot(
            highlight,
            x="Comp. " + str(comp),
            y=pop,
            hue="cytokine",
            s=150,
            ax=ax,
            legend=False,
        )

        # Add labels next to highlighted points
        for _, row in highlight.iterrows():
            ax.text(
                row["Comp. " + str(comp)] + 0.02,
                row[pop],
                row["cytokine"],
                fontsize=18,
                fontweight="bold",
                va="center",
            )
    else:
        sns.scatterplot(newDF, x="Comp. " + str(comp), y=pop, hue="cytokine", ax=ax)

    # Increase font sizes
    ax.set_xlabel(f"Comp. {comp} weights", fontsize=20, fontweight="bold")
    ax.set_ylabel(f"{pop} percentage", fontsize=20, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=18)

    # Add Pearson correlation text to plot
    ax.text(
        0.05,
        0.95,
        f"r = {pearson_coeff:.3f}\np = {p_value:.3e}",
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
