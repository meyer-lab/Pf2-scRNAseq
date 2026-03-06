"""
Plots factors for a set of genes as well as sum of weights for the genes in the set
"""


from ...imports import import_Parse
from ..common import getSetup, subplotLabel
from ..commonFuncs.plotGeneral import plot_fc, plot_fc_heatmap
from anndata import read_h5ad

# plots gene component factors for specifc subset of genes


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((10, 5), (1, 1))

    # Add subplot labels
    subplotLabel(ax)
 
    #X = read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")
    X=import_Parse(geneThreshold=0.0001,doublet=True)

    
    plot_fc_heatmap(
       X,
       "IL32",
       "PBS",
        ax[0],
        cellType="cell_type"
    )
    # plot_fc(
    #     X,
    #     "IRF1",
    #     "PBS",
    #     ax[0]
    # )
    

    return f
