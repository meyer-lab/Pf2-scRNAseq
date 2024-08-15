"""
CITEseq: Pf2 factors, weights, PaCMAP labeled by all conditions/leiden clusters,
and ratio of condition components based on days
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
    ax, f = getSetup((20, 4), (1, 3))
    subplotLabel(ax)

    X = anndata.read_h5ad("/opt/extra-storage/CRC/GSE178341/crc10x_full_50cmp.h5ad")
    
    samples_only_df = sample_names_only(X, "HistologicGradeSimple")
    

    # grouping_hgsgs = pd.Series(samples_only_df["HistologicGradeSimple"].to_numpy())

    # for i in grouping_hgsgs:
    #     print(i)
    # print(grouping_hgsgs)
    # print(np.unique(grouping_hgsgs))
    # plot_condition_factors(X, ax[0], condition_label="PID", cond_group_labels=grouping_hgsgs, groupConditions=True)


    
    # plot_condition_factors(X, ax[0], condition_label="HistologicGradeSimpleGradeSimple"), X.obs["MMRStatus"], groupConditions=True)
    # plot_eigenstate_factors(X, ax[1])
    # plot_gene_factors(X, ax[2])
    # plot_factor_weight(X, ax[3])

    # plot_labels_pacmap(X, "time", ax[4])


    return f


def sample_names_only(X: anndata.AnnData, label: str):
    """Obtain samples once only with corresponding observations"""
    samples = X.obs
    print(samples)
    
    unique_idx = np.unique(samples["condition_unique_idxs"])
    
    label_samples = np.empty(len(unique_idx), dtype=str)
    for i in range(20): 
        samples_idx = samples.loc[samples["condition_unique_idxs"] == i]
        print(samples_idx[label])
        print(np.unique(samples_idx[label]))
        label_samples[i] = str(np.unique(samples_idx[label]))
    
    print(label_samples)
        
        
    
    # df_samples = samples.drop_duplicates(subset="condition_unique_idxs")
    # df_samples = df_samples.sort_values("condition_unique_idxs")

    # return df_samples