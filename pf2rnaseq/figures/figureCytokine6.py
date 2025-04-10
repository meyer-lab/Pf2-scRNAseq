"""
Plots pacmap for a single component
"""

from anndata import read_h5ad
from ..imports import import_cytokine,import_pf2Cytokine30 
from .common import subplotLabel, getSetup
from .commonFuncs.plotPaCMAP import plot_wp_pacmap
import seaborn as sns
from .commonFuncs.plotGeneral import cell_count_perc_df


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((20, 3), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    X = import_pf2Cytokine30 ()

    cmp=16

    plot_wp_pacmap(X, cmp, ax[0], 1.0)

    return f
