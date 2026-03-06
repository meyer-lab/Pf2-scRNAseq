"""
Figure A4: Correlation Matrix
"""

import anndata
import numpy as np
import pandas as pd
import seaborn as sns

from .common import getSetup


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((30, 30), (2, 2))

    X = anndata.read_h5ad("/home/nicoleb/C3TAg_Pf2_30.h5ad")
    # X.uns["Pf2_A"] = correct_conditions(X)
    X.uns["Pf2_A"] -= np.min(X.uns["Pf2_A"], axis=0)
    X.uns["Pf2_A"] += np.median(X.uns["Pf2_A"], axis=0)
    X.uns["Pf2_A"] = np.log(X.uns["Pf2_A"])

    condition_factors_df = pd.DataFrame(
        data=X.uns["Pf2_A"],
        columns=[f"Cmp. {i}" for i in np.arange(1, X.uns["Pf2_A"].shape[1] + 1)],
    )

    pc_df = partial_correlation_matrix(condition_factors_df)

    f = plot_partial_correlation_matrix(pc_df, f)

    return f


def partial_correlation_matrix(df: pd.DataFrame):
    """Calculates partial correlation matrix"""
    cov_df = df.cov()
    vi = np.linalg.pinv(cov_df, hermitian=True)
    vi_diag = vi.diagonal()
    D = np.diag(np.sqrt(1 / vi_diag))
    pcor = -1 * (D @ vi @ D)
    pcor[np.diag_indices_from(pcor)] = 1
    pcor_df = pd.DataFrame(pcor, columns=cov_df.columns, index=cov_df.columns)

    return pcor_df


def plot_partial_correlation_matrix(df: pd.DataFrame, f):
    """Plots partial correlation matrix"""

    cmap = sns.color_palette("vlag", as_cmap=True)
    f = sns.clustermap(
        df,
        robust=True,
        vmin=-1,
        vmax=1,
        #    row_cluster=True,
        #    col_cluster=True,
        #    annot=True,
        cmap=cmap,
        figsize=(25, 25),
    )

    return f
