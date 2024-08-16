"""
CC: Pf2 factors
"""

from anndata import read_h5ad
from matplotlib.axes import Axes
import pandas as pd
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
    ax, f = getSetup((6, 6), (2, 2))
    subplotLabel(ax)

    X = anndata.read_h5ad("/opt/extra-storage/CRC/GSE178341/crc10x_full_50cmp.h5ad")
    
    samples_names = sample_names_only(X, "HistologicGradeSimple")


    plot_condition_factors(X, ax[0], condition_label="PID", cond_group_labels=pd.Series(samples_names), groupConditions=True)
    ax[0].yaxis.set_ticklabels([])
    plot_eigenstate_factors(X, ax[1])
    plot_gene_factors(X, ax[2])
    ax[2].yaxis.set_ticklabels([])


    return f


def sample_names_only(X: anndata.AnnData, label: str):
    """Obtain samples once only with corresponding observations"""
    samples = X.obs
    unique_idx = np.unique(samples["condition_unique_idxs"])
    label_samples = []
    
    for i in range(len(unique_idx)): 
        samples_idx = samples.loc[samples["condition_unique_idxs"] == i]
        
        if pd.isna(samples_idx[label].to_numpy()).any() == True:
            samples_idx_np = samples_idx[label].to_numpy()
            label_wo_nan = np.unique(samples_idx_np[~pd.isna(samples_idx_np)])
            label_w_nan = label_wo_nan + "-NaN"
            label_samples.append(label_w_nan[0])
        else:
            label_no_nan = np.unique(samples_idx[label])
            label_samples.append(label_no_nan[0])

    
    return label_samples