"""
Plotting the R2X
"""

from ..imports import import_Heiser
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import plot_r2x


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 10), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    X = import_Heiser()

    ranks = list(range(1, 31))
    plot_r2x(X, ranks, ax[0])

    return f
