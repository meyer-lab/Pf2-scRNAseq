"""
Cytokine: Highly weighted genes per component
"""

from ..factorization import correct_conditions, pf2
from ..imports import import_MouseImmune
from .common import getSetup
from .commonFuncs.plotFactors import plot_gene_factors_partial


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((21, 24), (10, 10))
    X = import_MouseImmune()
    X = pf2(X, 20, regularize_A=True, regParam=1e-5)

    X.uns["Pf2_A"] = correct_conditions(X)

    for i in range(X.uns["Pf2_A"].shape[1]):
        plot_gene_factors_partial(i + 1, X, ax[2 * i], geneAmount=10, top=True)
        plot_gene_factors_partial(i + 1, X, ax[(2 * i) + 1], geneAmount=5, top=False)

    return f
