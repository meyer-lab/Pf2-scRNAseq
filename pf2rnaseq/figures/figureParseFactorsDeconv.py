"""
Parse data: Plotting factors
"""

import numpy as np
import pandas as pd
import seaborn as sns
from anndata import read_h5ad
from matplotlib import pyplot as plt

from ..factorization import correct_conditions, deconvolution_cytokine
from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import (
    plot_condition_factors,
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
    ax, f = getSetup((22, 15), (1, 3))

    # Add subplot labels
    subplotLabel(ax)

    # Load data
    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X.uns["Pf2_A"] = correct_conditions(X)

    W, H = deconvolution_cytokine(X.uns["Pf2_A"], alpha=9e-5, max_iter=5000)


    # Get cytokine names in correct order
    samples_df = samples_only(X)

    # Create deconvolved version for plotting
    X_deconv = X.copy()
    X_deconv.uns["Pf2_A"] = H  # Use primary effects only

    plot_condition_factors(
        X_deconv,
        ax[0],
        samples_df["cytokine"],
        groupConditions=True,
        cond="cytokine",
        log_scale=False,
        centering=False
    )
    ax[0].set_title("Deconvolved matrix (H)", fontsize=12, fontweight="bold")

    plot_condition_factors(
        X,
        ax[1],
        samples_df["cytokine"],
        groupConditions=True,
        cond="cytokine",
        log_scale=False,
        centering=False
    )
    ax[1].set_title("Original Effects (A)", fontsize=12, fontweight="bold")

    cytokine_names = samples_df["cytokine"].values

    # Plot 2: W heatmap (primary effects)
    sns.heatmap(
        W,
        ax=ax[2],
        cmap="YlOrRd",
        robust=True,
        square=True,
        cbar_kws={"label": "Signaling Strength"},
        xticklabels=cytokine_names,
        yticklabels=cytokine_names,
    )
    ax[2].set_title("Cytokine Signaling (W)", fontsize=12, fontweight="bold")
    ax[2].set_xlabel("Inducing Cytokine →", fontsize=10)
    ax[2].set_ylabel("← Induced Cytokine", fontsize=10)
    plt.setp(ax[2].get_xticklabels(), rotation=90, ha="center", fontsize=6)
    plt.setp(ax[2].get_yticklabels(), rotation=0, fontsize=6)

    return f
