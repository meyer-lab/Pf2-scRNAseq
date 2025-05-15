"Plot condition pacmap"

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import plot_labels_pacmap


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    # X = read_h5ad("/home/nicoleb/C3TAg_Pf2_30.h5ad")
    X = read_h5ad("/home/nicoleb/'C3TAg_50")

    plot_labels_pacmap(X, "lineage", ax[0])
    # plot_labels_pacmap(X, "sc_lineage", ax[0])

    return f
