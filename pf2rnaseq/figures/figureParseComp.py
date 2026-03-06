"""
Parse data: Plotting factors
"""

import pandas as pd
from anndata import read_h5ad

from ..factorization import correct_conditions
from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import (
    plot_gene_expression_top_cells_by_condition,
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
    ax, f = getSetup((20, 14), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X.uns["Pf2_A"] = correct_conditions(X)
    # X = import_Parse(geneThreshold=0.01)
    # X = pf2(X, 80)
    # X.write_h5ad("/home/nicoleb/ParsePf2_80.h5ad")

    # plot_single_component_factors(X, ax[0], 56, cond="cytokine")
    # plot_single_cytokine_factors(X, ax[0], "IL-2", (1,100), cond="cytokine")

    # X_treg = X[X.obs["cell_type"] == "Treg", :].copy()
    plot_gene_expression_top_cells_by_condition(
        X, "IL2RA", cmp=56, ax=ax[0], top_percentile=1, aggregation="median"
    )
    # plot_gene_expression_by_wp_percentile(X_treg_il15, 'BCL2', cmp=56, ax=ax[0])
    # plot_eigenstate_factors(X, ax[0])
    # plot_gene_factors(X, ax[0])
    # plot_comp_weights(X, ax[2], 18, cond="cytokine")

    return f
