"""
Parse data: Plotting factors
"""

import numpy as np
import pandas as pd
import seaborn as sns
from anndata import read_h5ad
from matplotlib import pyplot as plt

from ..factorization import correct_conditions
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
    ax, f = getSetup((18, 28), (2, 2))

    # Add subplot labels
    subplotLabel(ax)
    """
    Load and plot saved W and H matrices from CSV files.
    
    Parameters
    ----------
    w_csv : str
        Path to W matrix CSV file
    h_csv : str
        Path to H matrix CSV file
    figsize : tuple
        Figure size (width, height)
    
    Returns
    -------
    f : matplotlib.figure.Figure
        The figure object
    """
    # Load the CSV files
    W_df = pd.read_csv(
        "/home/nicoleb/Pf2-scRNAseq-1/cytokine_crosstalk_W.csv", index_col=0
    )
    H_df = pd.read_csv(
        "/home/nicoleb/Pf2-scRNAseq-1/cytokine_primary_effects_H.csv", index_col=0
    )

    # X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X = read_h5ad("/home/nicoleb/ParsePf2_80.h5ad")
    X.uns["Pf2_A"] = correct_conditions(X)
    A = X.uns["Pf2_A"]
    has_negative = np.any(A < 0)
    min_val = A.min()

    n_negative = np.sum(A < 0)
    pct_negative = (n_negative / A.size) * 100

    print("\n=== A Matrix Statistics ===")
    print(f"  Has negative values: {has_negative}")
    print(f"  Min value: {min_val:.6f}")
    print(f"  Max value: {pct_negative:.6f}")

    X_deconv = X.copy()
    X_deconv.uns["Pf2_A"] = H_df  # Use primary effects only
    samples_df = samples_only(X)

    plot_condition_factors(
        X_deconv,
        ax[0],
        samples_df["cytokine"],
        groupConditions=True,
        cond="cytokine",
        log_scale=False,
    )
    plot_condition_factors(
        X,
        ax[1],
        samples_df["cytokine"],
        groupConditions=True,
        cond="cytokine",
        log_scale=False,
    )

    # Plot 2: W heatmap (primary effects)
    sns.heatmap(
        W_df,
        ax=ax[2],
        cmap="YlOrRd",
        robust=True,
        square=True,
        cbar_kws={"label": "Interaction Strength"},
        xticklabels=True,
        yticklabels=True,
    )
    ax[2].set_title("Cytokine Cross-talk (W)", fontsize=12, fontweight="bold")
    ax[2].set_xlabel("Inducing Cytokine →", fontsize=10)
    ax[2].set_ylabel("← Induced Cytokine", fontsize=10)
    plt.setp(ax[2].get_xticklabels(), rotation=90, ha="center", fontsize=6)
    plt.setp(ax[2].get_yticklabels(), rotation=0, fontsize=6)

    return f
