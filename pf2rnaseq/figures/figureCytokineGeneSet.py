"""
Plots factors for a set of genes as well as sum of weights for the genes in the set
"""

import pandas as pd
from anndata import read_h5ad
from .common import subplotLabel, getSetup
from .commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_eigenstate_factors,
    plot_geneSet_factors,
    plot_factor_weight,
    plot_geneSetScore,
    plot_geneSetScoreDot
)
from ..factorization import correct_conditions
from ..imports import import_cytokine,import_pf2Cytokine30 
from ..factorization import pf2

#plots gene component factors for specifc subset of genes
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
    X = import_cytokine()
    

    X = pf2(X, 15,regularize_A=True, regParam=5e-1)
    
    X.uns["Pf2_A"] = correct_conditions(X)
    
    immune_suppressive = ["FOXP3","TIGIT","ICOS","IL2RA","PDCD1","TGFB1","SOCS3","PD1", "LAG3","TNFRSF18","CTLA4"]  
    immune_activating = ["GZMA","GZMAB", "PRF1"] 
    plot_geneSet_factors(X,ax[0],immune_suppressive, False )
    #plot_geneSet_factors(X,ax[1],["GZMA","GZMB","PRF1",], False )
    plot_geneSetScore(X,ax[1],immune_suppressive, False )
    plot_geneSetScoreDot(X,ax[2],immune_suppressive, False )
    plot_geneSetScore(X,ax[3],immune_activating, False )

    

    return f