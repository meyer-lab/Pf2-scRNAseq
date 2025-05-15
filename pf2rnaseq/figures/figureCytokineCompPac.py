"""
Weighted projections per component in PaCMAP and boxplot
"""

import numpy as np
from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotPaCMAP import (
    plot_wp_pacmap,
    plot_wp_per_celltype,
)


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((35, 35), (10, 10))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")
    # X = import_cytokine()
    # X = pf2(X, 30, tolerance=1e-6)
    # X.uns["Pf2_A"] = correct_conditions(X)

    comps = np.arange(1, 31)

    for i, cmp in enumerate(comps):
        plot_wp_per_celltype(X, cmp, ax[2 * i], cellType="CellType2")
        plot_wp_pacmap(X, cmp, ax[2 * i + 1], cbarMax=0.25)

    return f
