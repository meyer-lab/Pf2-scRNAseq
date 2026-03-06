"""
Highly weighted genes per component
"""

from anndata import read_h5ad

from .common import getSetup
from .commonFuncs.plotFactors import plot_gene_factors_partial


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((5, 5), (1, 1))
    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")

    plot_gene_factors_partial(56, X, ax[0], geneAmount=5, top=True)

    return f
