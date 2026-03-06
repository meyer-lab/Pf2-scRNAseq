"Plot condition pacmap"

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import plot_gene_pacmap


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((15, 15), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad("/home/nicoleb/MouseImmune_pf2_80.h5ad")

    # plot_labels_pacmap(X, "treatment", ax[0])
    plot_gene_pacmap("Il1r2", "pf2", X, ax[0])

    return f
