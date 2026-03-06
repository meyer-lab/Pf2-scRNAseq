"Plot gene pacmap"

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import plot_gene_pacmap


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((5, 5), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")

    # X = import_cytokine()
    # X = pf2(X, 30)

    plot_gene_pacmap("IL2RA", "Pf2", X, ax[0])

    return f
