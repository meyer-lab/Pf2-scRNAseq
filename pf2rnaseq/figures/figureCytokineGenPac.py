"Plot gene pacmap"

import pandas as pd
from anndata import read_h5ad
from .common import subplotLabel, getSetup
from .commonFuncs.plotPaCMAP import (
   plot_gene_pacmap,
   plot_wp_pacmap,
   plot_labels_pacmap
)
from ..factorization import correct_conditions
from ..imports import import_pf2Cytokine30, import_cytokine
from ..factorization import pf2

def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((5, 5), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad("/home/nicoleb/Cytokine_Pf2__NB30_noreg.h5ad")
    #X = read_h5ad('/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Cytokine_Pf2_annotated_NB_new.h5ad')
    #X = import_cytokine()
    #X = pf2(X, 30, tolerance=1e-6)
    #plot_labels_pacmap(X, "Condition",ax[0], condition='control' )
    #plot_wp_pacmap(X, 9, ax[0])
    plot_gene_pacmap('IL2RA',"Pf2",X, ax[0])

    return f