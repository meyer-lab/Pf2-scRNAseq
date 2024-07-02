"""
XX:
"""

from anndata import read_h5ad
from ..imports import pseudobulk_lupus
from .common import subplotLabel, getSetup
import numpy as np

def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((30, 8), (2, 4))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/opt/andrew/lupus/lupus_fitted_ann.h5ad")
    bulk_matrix, bulk_tensor = pseudobulk_lupus(X)
    
    print(bulk_matrix)
    print(np.shape(bulk_tensor))


    return f



