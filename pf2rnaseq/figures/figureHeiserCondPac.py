"Plot condition pacmap"

import pandas as pd
from anndata import read_h5ad
from .common import subplotLabel, getSetup
from .commonFuncs.plotPaCMAP import (
   plot_labels_pacmap
)
from ..factorization import correct_conditions

from ..factorization import pf2

def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((15, 15), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad('/home/nicoleb/C3TAg_Pf2_30.h5ad')
    #X = import_cytokine()
    #X = pf2(X, 30, tolerance=1e-6)
    
    
    plot_labels_pacmap(X,'', ax[0])

    return f