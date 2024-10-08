"""
Cytokines: Plotting Cytokine factors and weights
"""

import pandas as pd
from .common import subplotLabel, getSetup
from .commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_eigenstate_factors,
    plot_gene_factors,
    plot_factor_weight,
)
from ..factorization import correct_conditions
from ..imports import import_cytokine
from ..factorization import pf2


def samples_only(X) -> pd.DataFrame:
    """Obtain samples once only with corresponding observations"""
    samples = X.obs
    df_samples = samples.drop_duplicates(subset="condition_unique_idxs")
    df_samples = df_samples.sort_values("condition_unique_idxs")
    return df_samples


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((8, 8), (2, 2))

    # Add subplot labels
    subplotLabel(ax)

    X = import_cytokine()
    X = pf2(X, 50, tolerance=1e-6)

    X.uns["Pf2_A"] = correct_conditions(X)
    stimulations = samples_only(X)["Condition"]

    plot_condition_factors(X, ax[0], stimulations, groupConditions=True)
    ax[0].set(yticks=[])
    plot_eigenstate_factors(X, ax[1])
    plot_gene_factors(X, ax[2])
    plot_factor_weight(X, ax[3])

    return f
