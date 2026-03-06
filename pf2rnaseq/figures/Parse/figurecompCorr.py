import anndata
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from scipy.stats import pearsonr

from ..common import getSetup, subplotLabel
from ..commonFuncs.plotGeneral import (
    cell_count_perc_df,
    rotate_xaxis,
)


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((18, 16), (5, 4))
    subplotLabel(ax)

    X = anndata.read_h5ad("/home/nicoleb/ParsePf2_100_D11_filt.h5ad")

    celltype_count_perc_df = cell_count_perc_df(X, celltype="cell_type")
    cmps = [56]
    for i, cmp in enumerate(cmps):
        plot_correlation_cmp_cell_count_perc(
            X, cmp, celltype_count_perc_df, ax[i + 2], cellPerc=False
        )

    return f


def plot_correlation_cmp_cell_count_perc(
    X: anndata, cmp: int, cellcountDF: pd.DataFrame, ax: Axes, cellPerc=True
):
    """Plot component weights by cell type count or percentage for a cell type"""
    yt = np.unique(X.obs["cytokine"])
    factorsA = np.array(X.uns["Pf2_A"])
    factorsA = factorsA[:, cmp - 1]
    cellPerc = "Cell Type Percentage" if cellPerc is True else "Cell Count"
    totaldf = pd.DataFrame([])
    correlationdf = pd.DataFrame([])
    cellcountDF["cytokine"] = pd.Categorical(cellcountDF["cytokine"], yt)
    for celltype in np.unique(cellcountDF["Cell Type"]):
        for j, cond in enumerate(np.unique(cellcountDF["cytokine"])):
            smalldf = cellcountDF.loc[
                (cellcountDF["cytokine"] == cond)
                & (cellcountDF["Cell Type"] == celltype)
            ]
            if smalldf.empty is False:
                smalldf = smalldf.assign(Cmp=factorsA[j])
            else:
                smalldf = pd.DataFrame(
                    {
                        "cytokine": [cond],
                        "Cell Type": [celltype],
                        cellPerc: [0],
                        "Cmp": [factorsA[j]],
                    }
                )

            totaldf = pd.concat([totaldf, smalldf])

        df = totaldf.loc[totaldf["Cell Type"] == celltype]
        pearson = pearsonr(df["Cmp"], df[cellPerc])[0]

        correlationdf = pd.concat(
            [
                correlationdf,
                pd.DataFrame(
                    {
                        "Cell Type": [celltype],
                        "Correlation": ["Pearson"],
                        "Value": [pearson],
                    }
                ),
            ]
        )

    sns.swarmplot(data=correlationdf, y="Value", hue="Correlation", ax=ax)
    rotate_xaxis(ax)
    ax.set(title=f"Cmp. {cmp} V. {cellPerc}")
