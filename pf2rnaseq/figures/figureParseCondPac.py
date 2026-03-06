"Plot condition pacmap"

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import plot_labels_pacmap, plot_gene_pacmap


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((15, 15), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    # X=import_Parse(geneThreshold=0.01, doublet=True)
    # X=pf2(X,100)
    # X.write_h5ad("/home/nicoleb/ParsePf2_100_D11_neighbors.h5ad")
    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")

    #plot_labels_pacmap(X, "cytokine", ax[0])
    plot_gene_pacmap("IL2RG", "pf2",X, ax[0])

    return f
