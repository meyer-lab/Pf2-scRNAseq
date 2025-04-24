"""
Weighted projections per component in PaCMAP and boxplot
"""

import pandas as pd

import numpy as np
from .common import subplotLabel, getSetup
from anndata import read_h5ad
from .commonFuncs.plotPaCMAP import (
    plot_wp_per_celltype,
    plot_wp_pacmap,
   
)
from ..factorization import correct_conditions
from ..imports import import_Heiser
from ..factorization import pf2



def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((35, 35), (10, 10))

    # Add subplot labels
    subplotLabel(ax)
   
    
    X=import_Heiser()
    X = pf2(X, 30)
    X.uns["Pf2_A"] = correct_conditions(X)
    
    

    
    comps = np.arange(1,31)

    for i, cmp in enumerate(comps):
        plot_wp_per_celltype(X, cmp, ax[2 * i], cellType="cell_type_1")
        plot_wp_pacmap(X, cmp, ax[2 * i + 1], cbarMax=0.25)

    return f
