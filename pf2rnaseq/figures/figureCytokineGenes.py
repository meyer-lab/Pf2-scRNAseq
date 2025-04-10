"""
Cytokine: Highly weighted genes per component
"""
from anndata import read_h5ad
from .common import getSetup
from .commonFuncs.plotFactors import plot_gene_factors_partial
from ..imports import import_cytokine,import_pf2Cytokine30
from ..factorization import correct_conditions
from ..factorization import pf2


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((21, 24), (10, 10))
    X = import_cytokine()
    #X=read_h5ad("/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Cytokine_Pf2__NB30.h5ad")
    #X = import_pf2Cytokine30()
    X = pf2(X, 30, regParam=1e-5, regularize_A=True)

    X.uns["Pf2_A"] = correct_conditions(X)

    for i in range(X.uns["Pf2_A"].shape[1]):
        plot_gene_factors_partial(i + 1, X, ax[2 * i], geneAmount=10, top=True)
        plot_gene_factors_partial(i + 1, X, ax[(2 * i) + 1], geneAmount=5, top=False)

    return f
