"""plot average gene expression for a set of genes for each condition"""""
import pandas as pd
from anndata import read_h5ad
from .common import subplotLabel, getSetup
from .commonFuncs.plotGeneral import (
   plot_gene_set_expression
  
)
from ..factorization import correct_conditions
from ..imports import import_pf2Cytokine30, import_cytokine
from ..factorization import pf2



def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((15, 15), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    
    X = import_pf2Cytokine30()
    
    immune_suppressive = ["FOXP3","TIGIT","ICOS","IL2RA","INDO","TGFB1","SOCS3","SOCS1" "MCL1","BCL2", "TIM-3","TNFRSF18"]  
    plot_gene_set_expression(X, immune_suppressive,ax[0])
    

    return f