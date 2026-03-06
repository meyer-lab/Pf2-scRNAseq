"""
Plots factors for a set of genes as well as sum of weights for the genes in the set
"""

import pandas as pd

from ..factorization import correct_conditions, pf2
from ..imports import import_MouseImmune
from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import (
    plot_geneSet_factors,
    plot_geneSetScore,
    plot_geneSetScoreDot,
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
    ax, f = getSetup((15, 15), (2, 2))

    # Add subplot labels
    subplotLabel(ax)
    X = import_MouseImmune()

    X = pf2(X, 20, regularize_A=True, regParam=1e-5)

    X.uns["Pf2_A"] = correct_conditions(X)

    immune_suppressive = [
        "Foxp3",
        "Tigit",
        "Icos",
        "Il2ra",
        "Pdcd1",
        "Tgfb1",
        "Socs3",
        "Pd1",
        "Lag3",
        "Tnfrsf18",
        "Ctla4",
    ]
    immune_activating = ["Gzma", "Gzmb", "Prf1"]
    plot_geneSet_factors(X, ax[0], immune_suppressive, False)
    # plot_geneSet_factors(X,ax[1],["GZMA","GZMB","PRF1",], False )
    plot_geneSetScore(X, ax[1], immune_suppressive, False)
    plot_geneSetScoreDot(X, ax[2], immune_suppressive, False)
    plot_geneSetScore(X, ax[3], immune_activating, False)

    return f
