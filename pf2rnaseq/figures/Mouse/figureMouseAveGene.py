"""
Plots factors for a set of genes as well as sum of weights for the genes in the set
"""

from ..imports import import_MouseImmune
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import (
    plot_avegene_per_category,
)

# plots gene component factors for specifc subset of genes


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((25, 15), (2, 2))

    # Add subplot labels
    subplotLabel(ax)

    X = import_MouseImmune(geneThreshold=0.0)
    plot_avegene_per_category(
        ["IL1b", "IL15", "IL2"], "Tnf", X, ax[0], mean=True, cellType="celltype"
    )
    plot_avegene_per_category(
        ["IL1b", "IL15", "IL2"], "Foxp3", X, ax[1], mean=True, cellType="celltype"
    )
    plot_avegene_per_category(
        ["IL1b", "IL15", "IL2"], "Ctla4", X, ax[2], mean=True, cellType="celltype"
    )
    plot_avegene_per_category(
        ["IL1b", "IL15", "IL2"], "Ifng", X, ax[3], mean=True, cellType="celltype"
    )

    # plot_avegene_per_category(["IL10"], "Tnf", X, ax[0], mean=True, cellType="celltype")

    # plot_avegene_per_celltype(X, "IL1R1", ax[0], cellType="cell_type")
    # X_genes = X[:, ["FOXP3", "CTLA4"]].to_memory()
    # X_genes = X_genes[X_genes.obs["cell_type"] == "Treg", :]
    # gene_plot_cells(
    #     X_genes, unique=["IL-1-beta"], hue="cytokine", ax=ax[0], kde=False, cellType="cell_type"
    # )

    return f
