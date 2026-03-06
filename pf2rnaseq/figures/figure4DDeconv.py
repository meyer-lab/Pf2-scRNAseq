"""
Parse data: Plotting factors
"""

import numpy as np
import pandas as pd
import seaborn as sns
from anndata import read_h5ad
from matplotlib import pyplot as plt

from ..factorization import correct_conditions, deconvolution_cytokine_admm
from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_comp_weights
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
    ax, f = getSetup((25, 15), (1, 3))

    # Add subplot labels
    subplotLabel(ax)

    # Load data
    A = np.load("/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/rank200_tensor_cytokine_factor_final.npy")
    samples_df = np.load("/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/cytokines_order.npy", allow_pickle=True)

    # Center A by cytokine medians
    cytokine_medians = np.median(A, axis=1, keepdims=True)
    A_centered = A - cytokine_medians

    W, H, _ = deconvolution_cytokine_admm(A_centered, alpha_h=0.01, alpha_w=0.001, rho=0.5, non_negative_w=True, adaptive_rho=True)

    # Get cytokine names from samples_df
    if isinstance(samples_df, np.ndarray):
        # If it's a numpy array, convert to appropriate format
        cytokine_names = samples_df
    else:
        cytokine_names = samples_df["cytokine"].values
    
    # Plot H heatmap
    n_components = H.shape[1]
    component_labels = [f"C{i}" for i in range(n_components)]
    
    sns.heatmap(
        H,
        ax=ax[0],
        cmap="RdBu_r",
        robust=True,
        center=0,
        cbar_kws={"label": "Weight"},
        xticklabels=component_labels,
        yticklabels=cytokine_names,
    )
    ax[0].set_title("Deconvolved Effects (H)", fontsize=12, fontweight="bold")
    ax[0].set_xlabel("Component", fontsize=10)
    ax[0].set_ylabel("Cytokine", fontsize=10)
    plt.setp(ax[0].get_xticklabels(), rotation=90, ha="center", fontsize=8)
    plt.setp(ax[0].get_yticklabels(), rotation=0, fontsize=8)

    # Plot A_centered heatmap for comparison
    sns.heatmap(
        A_centered,
        ax=ax[1],
        cmap="RdBu_r",
        robust=False,
        center=0,
        cbar_kws={"label": "Weight"},
        xticklabels=component_labels,
        yticklabels=cytokine_names,
    )
    ax[1].set_title("Original Effects (A centered)", fontsize=12, fontweight="bold")
    ax[1].set_xlabel("Component", fontsize=10)
    ax[1].set_ylabel("Cytokine", fontsize=10)
    plt.setp(ax[1].get_xticklabels(), rotation=90, ha="center", fontsize=8)
    plt.setp(ax[1].get_yticklabels(), rotation=0, fontsize=8)

    # Plot W heatmap (cytokine interactions)
    sns.heatmap(
        W,
        ax=ax[2],
        cmap="RdBu_r",
        robust=False,
        center=0,
        square=True,
        cbar_kws={"label": "Signaling Strength"},
        xticklabels=cytokine_names,
        yticklabels=cytokine_names,
    )
    ax[2].set_title("Cytokine Signaling (W)", fontsize=12, fontweight="bold")
    ax[2].set_xlabel("Induced Cytokine", fontsize=10)
    ax[2].set_ylabel("Inducing Cytokine", fontsize=10)
    plt.setp(ax[2].get_xticklabels(), rotation=90, ha="center", fontsize=6)
    plt.setp(ax[2].get_yticklabels(), rotation=0, fontsize=6)
    f.savefig(f"figure4D.png", dpi=300, bbox_inches="tight")

    return f