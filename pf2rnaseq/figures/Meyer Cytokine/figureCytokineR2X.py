"""
Plotting the R2X
"""

from ..factorization import pf2
from ..imports import import_cytokine
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import plot_r2x


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    X = import_cytokine()

    X = pf2(X, 30)
    ranks = list(range(1, 31))
    plot_r2x(X, ranks, ax[0])

    return f
