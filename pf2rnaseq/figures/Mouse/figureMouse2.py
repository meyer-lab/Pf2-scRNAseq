"""
Mouse immune dictionary: Plotting factors
"""

from ..imports import import_MouseImmune
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import gene_plot_cells


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((25, 15), (2, 2))

    # Add subplot labels
    subplotLabel(ax)

    X = import_MouseImmune(geneThreshold=0.0)
    X_genes = X[:, ["Foxp3", "Il1r2"]].to_memory()
    X_genes = X_genes[X_genes.obs["celltype"] == "Treg", :]
    gene_plot_cells(
        X_genes,
        unique=["IL1b", "IL1a"],
        hue="cyt",
        ax=ax[0],
        kde=False,
        cellType="celltype",
    )

    return f
