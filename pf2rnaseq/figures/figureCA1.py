"""
CITEseq: Pf2 factors, weights, PaCMAP labeled by all conditions/leiden clusters,
and ratio of condition components based on days
"""

from anndata import read_h5ad
from matplotlib.axes import Axes
import anndata
from .common import subplotLabel, getSetup
from .commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_eigenstate_factors,
    plot_gene_factors,
    plot_factor_weight,
)
from .commonFuncs.plotPaCMAP import plot_labels_pacmap
import numpy as np


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((30, 8), (2, 4))

    # Add subplot labels
    subplotLabel(ax)


    #Anndata
    plot_condition_factors(X, ax[0])
    plot_eigenstate_factors(X, ax[1])
    plot_gene_factors(X, ax[2])
    plot_factor_weight(X, ax[3])

    plot_labels_pacmap(X, "Condition", ax[4])


    return f

