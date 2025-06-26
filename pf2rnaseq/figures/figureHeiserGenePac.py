"""Plot gene pacmap"""

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import plot_gene_pacmap


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((5, 5), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/C3TAg_Pf2_30.h5ad")

    plot_gene_pacmap("Folr2", "Pf2", X, ax[0])

    return f
