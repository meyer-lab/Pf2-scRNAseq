"""plot average gene expression per cell type for all conditions"""""
import pandas as pd
from anndata import read_h5ad
from .common import subplotLabel, getSetup
from .commonFuncs.plotGeneral import (
    
    plot_avegene_per_celltype
)
from ..factorization import correct_conditions
from ..imports import import_pf2Cytokine30, import_cytokine
from ..factorization import pf2
import scanpy as sc





def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((15, 15), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")
    #X = read_h5ad('/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Cytokine_Pf2_annotated_NB_new.h5ad')
    #X = import_cytokine()
    #X = pf2(X, 30)
    #X.write_h5ad("Cytokine_Pf2__NB30_022725.h5ad")
    immune_suppressive = ["FOXP3","TIGIT"]  
    #X.var_names_make_unique()
    conds=['IL10_2000nM', 'IL10_500nM', 'IL2_100pM', 'IL2_10nM', 'IL2_1nM',
       'IL2_200nM', 'IL2_50nM', 'IL7_100nM', 'TGFB_10nM', 'TGFB_50nM',
       'control']
    plot_avegene_per_celltype(X, "FOXP3", ax[0], 'CellType2')
    

    return f