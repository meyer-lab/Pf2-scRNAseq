"""
Plots average gene expression for a specific gene across cytokine conditions, faceted by cell type
"""


from ..imports import import_Parse
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import plot_avegene_per_category
from anndata import read_h5ad

# plots gene component factors for specifc subset of genes


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((15, 10), (2, 2))

    # Add subplot labels
    subplotLabel(ax)
    # X = read_h5ad("/home/nicoleb/ParsePf2_100_D11.h5ad")
    #X =read_h5ad("/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Parse_Donor11.h5ad")    
    #X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_01.h5ad")

    
    plot_avegene_per_category(
        ["IL-15","IL-1-beta","IFN-epsilon","IL-2","PBS"],
        ["FOXP3"],
        X,
        ax[0],
        mean=False,
        cellType="cell_type",
        center_data=False,
    )

    plot_avegene_per_category(
        ["IL-15","IL-1-beta","IFN-epsilon","IL-2","PBS"],
        ["CTLA4"],
        X,
        ax[1],
        mean=False,
        cellType="cell_type",
        center_data=False,
    )
    

   
    return f
