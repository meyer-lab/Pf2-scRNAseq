"""
Heiser: Plotting factors
"""

import pandas as pd
import numpy as np
import scipy.sparse as sp
from anndata import read_h5ad, AnnData
import anndata as an
from .common import subplotLabel, getSetup
from .commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_condition_factors_groups,
    plot_eigenstate_factors,
    plot_gene_factors,
    plot_factor_weight,
)
from ..factorization import correct_conditions
from ..imports import import_Heiser
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
    ax, f = getSetup((15, 15), (2, 2))

    # Add subplot labels
    subplotLabel(ax)
    X= read_h5ad('/home/nicoleb/Heiser_minimal.h5ad')
    #X = import_Heiser()
    
    #X = pf2(X, 15)
    
    #X.uns["Pf2_A"] = correct_conditions(X)
    
    



    
    
    stimulations = samples_only(X)["treatment"]
    tumors = samples_only(X)["tumorType"]
    print(stimulations)
    
    #plot_condition_factors(X, ax[0], stimulations, cond="sample_id", groupConditions=True)
    plot_condition_factors_groups(X, ax[0], stimulations, tumors, cond="sample_id", groupConditions=True)
    plot_eigenstate_factors(X, ax[1])
    plot_gene_factors(X, ax[2])
    plot_factor_weight(X, ax[3])

    return f




