"""
Plots factors for a set of genes as well as sum of weights for the genes in the set
"""

import pandas as pd
from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import (
    plot_geneSet_factors,
    plot_geneSetScore,
)


# plots gene component factors for specifc subset of genes
def samples_only(X) -> pd.DataFrame:
    """Obtain samples once only with corresponding observations"""
    samples = X.obs
    df_samples = samples.drop_duplicates(subset="condition_unique_idxs")
    df_samples = df_samples.sort_values("condition_unique_idxs")
    return df_samples


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 20), (2, 1))

    # Add subplot labels
    subplotLabel(ax)

    # X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")

    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_01.h5ad")


    immune_suppressive = [
        "FOXP3",
        "TIGIT",
        "ICOS",
        "IL2RA",
        "PDCD1",
        "TGFB1",
        "SOCS3",
        "CD279",
        "LAG3",
        "TNFRSF18",
        "CTLA4",
        "TNSFSF10",
        "CD274",
        "IL10",
        "IKZF2",
        "GITR",
        "FAS",
    ]



    plot_geneSet_factors(X, ax[0], immune_suppressive, False)
    plot_geneSetScore(X, ax[1], immune_suppressive, False)

    return f
