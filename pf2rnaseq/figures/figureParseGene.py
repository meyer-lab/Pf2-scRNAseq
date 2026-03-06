"""
Parse data: Plotting factors
"""

import pandas as pd
from anndata import read_h5ad

from ..factorization import correct_conditions
from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import plot_component_top_genes_heatmap


def samples_only(X) -> pd.DataFrame:
    """Obtain samples once only with corresponding observations"""
    samples = X.obs
    df_samples = samples.drop_duplicates(subset="condition_unique_idxs")
    df_samples = df_samples.sort_values("condition_unique_idxs")
    return df_samples


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_01.h5ad")
    # X=import_Parse(geneThreshold=0.06)
    #X=pf2(X,80)

    X.uns["Pf2_A"] = correct_conditions(X)

    plot_component_top_genes_heatmap(X, ax[0], comp=56)

    return f
