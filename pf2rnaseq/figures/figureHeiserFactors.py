"""
Heiser: Plotting factors
"""

import numpy as np
import pandas as pd

from ..factorization import pf2
from ..imports import import_Heiser
from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import (
    plot_condition_factors_groups,
    plot_eigenstate_factors,
    plot_factor_weight,
    plot_gene_factors,
)


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
    X = import_Heiser(deviance=True)
    print(f"Data shape: {X.X.shape}")
    print(f"Data type: {X.X.dtype}")
    print(f"Contains NaN: {np.isnan(X.X).any()}")
    print(f"Contains inf: {np.isinf(X.X).any()}")
    print(f"Min value: {X.X.min()}")
    print(f"Max value: {X.Xmax()}")

    X = pf2(X, 20)

    stimulations = samples_only(X)["treatment"]
    tumors = samples_only(X)["expBatch"]
    print(stimulations)

    plot_condition_factors_groups(
        X, ax[0], stimulations, tumors, cond="sample_id", groupConditions=True
    )
    plot_eigenstate_factors(X, ax[1])
    plot_gene_factors(X, ax[2])
    plot_factor_weight(X, ax[3])

    return f
