"""
Plots factors for a set of genes as well as sum of weights for the genes in the set
"""

from ..imports import import_Parse
from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import plot_avegene_per_celltype

# plots gene component factors for specifc subset of genes


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((25, 10), (2, 1))

    # Add subplot labels
    subplotLabel(ax)
    # X = read_h5ad("/home/nicoleb/ParsePf2_100_D11.h5ad")
    X = import_Parse(geneThreshold=0.001, doublet=True)

    # X_filt = X[X.obs["cytokine"] == "IL-1-beta", :].copy()
    X_filt = X[X.obs["cytokine"] == "IL-12", :].copy()
    X_filt15 = X[X.obs["cytokine"] == "IL-15", :].copy()


    plot_avegene_per_celltype(X_filt, "IFNG", ax[0], cellType="cell_type", center_data=True)
    plot_avegene_per_celltype(X_filt15, "IFNG", ax[1], cellType="cell_type", center_data=True)

    # X_genes = X[:, ["FOXP3", "CTLA4"]].to_memory()
    # X_genes = X_genes[X_genes.obs["cell_type"] == "Treg", :]
    # gene_plot_cells(
    #    X_genes, unique=["IL-1-beta"], hue="cytokine", ax=ax[0], kde=False, cellType="cell_type"
    # )

    return f
