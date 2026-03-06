"""
Weighted projections per component as boxplot of cell types
"""

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import (
    plot_wp_per_celltype,
    plot_wp_pacmap,
)


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((8, 5), (1, 2))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")


    plot_wp_per_celltype(X, 56, ax[0], cellType="cell_type")
    plot_wp_pacmap(X, 56, ax[1], cbarMax=0.25)

    return f
