"""
factorization score

"""

from anndata import read_h5ad

from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import plot_fms_percent_drop_counts


def makeFigure():
    ax, f = getSetup((6, 3), (1, 2))
    subplotLabel(ax)
    # Using our cytokine dataset
    X = read_h5ad("/opt/extra-storage/Treg_h5ads/Treg_raw.h5ad")

    # Remove multiplexing identifiers
    X = X[:, ~X.var_names.str.match("^CMO3[0-9]{2}$")]  # type: ignore
    # Remove genes with too few reads now
    X = X[X.X.sum(axis=1) > 10, X.X.mean(axis=0) > 0.1]
    X = X.copy()
    percentList = [0.0, 25.0, 50.0, 75.0, 95.0, 99.0]
    plot_fms_percent_drop_counts(X, ax[0], percentList, runs=2, rank=20, deviance=True)

    return f
