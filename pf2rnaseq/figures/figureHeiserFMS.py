"""
factorization score
factorization score
"""

import numpy as np

from ..imports import import_Heiser
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import plot_fms_diff_ranks, plot_fms_percent_drop


def makeFigure():
    ax, f = getSetup((6, 3), (1, 2))
    subplotLabel(ax)

    X = import_Heiser()
    percentList = np.arange(0.0, 55.0, 5.0)
    plot_fms_percent_drop(X, ax[0], percentList=percentList, runs=2)

    ranks = np.arange(10, 101, 10)
    plot_fms_diff_ranks(X, ax[1], ranksList=ranks, runs=2)

    return f
