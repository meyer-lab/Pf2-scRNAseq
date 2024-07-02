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
import pandas as pd

def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((30, 8), (2, 4))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/opt/andrew/lupus/lupus_fitted_ann.h5ad")
    bulk_matrix, bulk_tensor = pseudobulk_lupus(X)
    
    print(bulk_matrix)

    
        


    return f

def pseudobulk_lupus(X, cellType="Cell Type"):
    """Average gene expression for each condition and cell type; 
    creates matrix and tensor version"""
    X_df = X.to_df()
    X_df = X_df.subtract(X.var["means"].values)
    X_df["Condition"] = X.obs["Condition"].values
    X_df["Cell Type"] = X.obs[cellType].values
    
    X_matrix = X_df.groupby(["Condition", "Cell Type"], observed=False).mean().reset_index()
    
    
    # connet status too
    
    conds = pd.unique(X_matrix["Condition"])
    celltypes = pd.unique(X_matrix["Cell Type"])
    genes = X.var_names.values            
    
    X_tensor = np.empty((len(conds), len(celltypes), len(genes)))
    X_tensor[:] = np.nan
    
    for i, cond in enumerate(conds):
        for j, celltype in enumerate(celltypes):
            specific_df = X_matrix.loc[(X_matrix["Condition"] == cond) & (X_matrix["Cell Type"] == celltype)] 
            X_tensor[i, j, :] = specific_df.iloc[0, 2:].to_numpy()
          
    return X_df, X_tensor



