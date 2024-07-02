"""
CITEseq: Pf2 factors, weights, PaCMAP labeled by all conditions/leiden clusters,
and ratio of condition components based on days
"""

from anndata import read_h5ad
from matplotlib.axes import Axes
import anndata
from .common import subplotLabel, getSetup
from .commonFuncs.plotFactors import (
    plot_condition_factors,
    plot_eigenstate_factors,
    plot_gene_factors,
    plot_factor_weight,
)
from .commonFuncs.plotPaCMAP import plot_labels_pacmap
import numpy as np


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((30, 8), (2, 4))

    # Add subplot labels
    subplotLabel(ax)

    X = read_h5ad("/opt/andrew/lupus/lupus_fitted_ann.h5ad")
    
    
    # bulk_de_genes_up_list = bulk_de_genes['Gene'].tolist()
    # # Subset the data based on the list of genes
    # adata2 = adata[:, adata.var_names.isin(bulk_de_genes_up_list)]
    # average_expression = adata2.X.mean(axis=1)
    # adata2.obs['bulk_de_gene_average'] = average_expression
    # sc.pl.umap(adata2, color='bulk_de_gene_average', cmap='viridis')
        
    

  

    return f



