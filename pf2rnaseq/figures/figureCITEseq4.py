"""
CITEseq: Highly weighted genes per component
"""

from anndata import read_h5ad
from .common import (
    subplotLabel,
    getSetup,
)
from .commonFuncs.plotFactors import plot_gene_factors_partial


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (6, 5))

    X = read_h5ad("/opt/extra-storage/CRC/GSE178341/crc10x_full_50cmp.h5ad", backed="r")
    comps = range(20, 50)

    for i, cmp in enumerate(comps):
        plot_gene_factors_partial(cmp, X, ax[i - 1], geneAmount=15)

    return f
