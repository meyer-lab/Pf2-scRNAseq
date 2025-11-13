"""
Weighted projections per component in PaCMAP and boxplot of cell types
"""

import numpy as np

from ..factorization import correct_conditions, pf2
from ..imports import import_Heiser
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

    X = import_Heiser()
    X = pf2(X, 30)
    X.uns["Pf2_A"] = correct_conditions(X)

    comps = np.arange(1, 31)

    for i, cmp in enumerate(comps):
        plot_wp_per_celltype(X, cmp, ax[2 * i], cellType="cell_type_1")
        plot_wp_pacmap(X, cmp, ax[2 * i + 1], cbarMax=0.25)

    return f
