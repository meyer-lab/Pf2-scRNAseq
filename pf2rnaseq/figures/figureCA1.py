"""
CITEseq: Pf2 factors, weights, PaCMAP labeled by all conditions/leiden clusters,
and ratio of condition components based on days
"""

from anndata import read_h5ad
from matplotlib.axes import Axes
import anndata
from pf2rnaseq.figures.common import subplotLabel, getSetup
from pf2rnaseq.figures.commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_eigenstate_factors,
    plot_gene_factors,
    plot_factor_weight,
)
from pf2rnaseq.figures.commonFuncs.plotPaCMAP import plot_labels_pacmap
import numpy as np
from pf2rnaseq.imports import import_CA, prepare_dataset_bulk

from pf2rnaseq.factorization import pf2

def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((30, 8), (2, 4))

    # Add subplot labels
    subplotLabel(ax)


    #Import the CA data as a DataFrame and put into an Anndata object
    CA_ad = import_CA()
    X = prepare_dataset_bulk(CA_ad, "sample_id", geneThreshold=0.001)

    X = pf2(X, rank=10)

    plot_condition_factors(X, ax[0], condition_name="sample_id")
    plot_eigenstate_factors(X, ax[1])
    plot_gene_factors(X, ax[2])
    plot_factor_weight(X, ax[3])

    # plot_labels_pacmap(X, "time", ax[4])


    return f