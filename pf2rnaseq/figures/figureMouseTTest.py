"""
Mouse cytokine Component ANOVA/TTest: Heatmap showing cytokine significance across components
"""

from anndata import read_h5ad

from ..factorization import correct_conditions
from .common import getSetup
from .commonFuncs.plotFactors import plot_ttest


def makeFigure():
    """Create heatmap figure showing dominant cytokines across components."""
    ax, f = getSetup((12, 8), (1, 1))

    X = read_h5ad("/home/nicoleb/Pf2-scRNAseq-1/MouseImmune_pf2_45dev.h5ad")
    X.uns["Pf2_A"] = correct_conditions(X)

    plot_ttest(X, ax[0])

    return f
