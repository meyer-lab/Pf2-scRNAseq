"""
Plot the score for each cell type using scanpy score genes
"""

import anndata
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import read_h5ad
from matplotlib.axes import Axes

from .common import getSetup, subplotLabel
from .commonFuncs.plotGeneral import rotate_xaxis

# from .figure4e_k import plot_correlation_cmp_cell_count_perc


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((8, 8), (1, 1))
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")
    X = suppressive_score(X)
    plot_score(X, ax[0], cellType="CellType2")
    ax[0].set(ylabel="Cytotoxic Score")

    # plot_gene_pacmap("RETN", "Pf2", X, ax[4])

    # celltype_count_perc_df = cell_count_perc_df(X, celltype="leiden", status=True)

    return f


def suppressive_score(X: anndata.AnnData):
    """Scanpy average gene score for all cells"""
    X.var_names_make_unique()
    immune_suppressive = [
        "FOXP3",
        "TIGIT",
        "ICOS",
        "IL2RA",
        "TGFB1",
        "SOCS3",
        "BCL2",
        "TNFRSF18",
    ]

    X = sc.tl.score_genes(
        adata=X, gene_list=immune_suppressive, copy=True, use_raw=False
    )

    return X


def plot_score(X: anndata.AnnData, ax: Axes, cellType: str = "Cell Type"):
    """Plots distribution of scores across cell types and conditions"""
    # Create DataFrame with individual cell scores
    df = pd.DataFrame(
        {
            "Score": X.obs["score"].values,
            "Condition": X.obs["Condition"].values,
            "Cell Type": X.obs[cellType].values,
        }
    )

    # Optional: You can filter for specific cell types if needed
    # df = df[df["Cell Type"].isin(["CD8+ T cells", "T reg", "CD4+ T cells"])]

    # Create the boxplot with all individual values
    sns.boxplot(
        data=df,
        x="Cell Type",
        y="Score",
        hue="Condition",
        order=sorted(df["Cell Type"].unique()),
        ax=ax,
        showfliers=False,
    )

    rotate_xaxis(ax)

    # Improve legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[: len(df["Condition"].unique())],
        labels[: len(df["Condition"].unique())],
    )
