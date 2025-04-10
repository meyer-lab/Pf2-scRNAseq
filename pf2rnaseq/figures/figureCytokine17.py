"""plotting average gene expression cell types"""  

import anndata
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from anndata import read_h5ad
from .common import getSetup, subplotLabel



def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((8, 8), (1, 1))
    subplotLabel(ax)

    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")


    immune_suppressive = ["FOXP3","TIGIT"]  

    plot_avegene_per_status_per_cluster( X, "CTLA4", ax[0], clusterName1='CD4 T', cellType="CellType2" )

    

    return f


def plot_avegene_per_status_per_cluster(
    X: anndata.AnnData,
    gene: str,
    ax: Axes,
    clusterName1: str,
    clusterName2=None,
    cellType: str = "Cell Type",
):
    """Plots average gene expression across cell types for a category of drugs"""
    genesV = X[:, gene]
    dataDF = genesV.to_df()
    dataDF = dataDF.subtract(genesV.var["means"].values)
    
    dataDF["Condition"] = genesV.obs["Condition"].values
    dataDF["Cell Type"] = genesV.obs[cellType].values

    df = pd.melt(
        dataDF, id_vars=["Cell Type", "Condition"], value_vars=gene
    ).rename(columns={"variable": "Gene", "value": "Value"})

    df = df.groupby(["Cell Type", "Gene", "Condition"], observed=False).mean()
    df = df.rename(columns={"Value": "Average Gene Expression"}).reset_index()

    if clusterName2 is None:
        dfClust = df.loc[df["Cell Type"] == clusterName1]
        clust_list = dfClust["Cell Type"].to_numpy()
        dfOther = df.loc[df["Cell Type"] != clusterName1]
        other_list = np.repeat("Other", dfOther.shape[0])

        dfClust = pd.concat([dfClust, dfOther]).reset_index(drop=True)
        dfClust["Cell Type"] = np.concatenate([clust_list, other_list])

    else:
        dfClust = df.loc[
            (df["Cell Type"] == clusterName1) & (df["Cell Type"] == clusterName2)
        ]
        clust_list = dfClust["Cell Type"].to_numpy()
        dfOther = df.loc[
            (df["Cell Type"] != clusterName1) & (df["Cell Type"] != clusterName2)
        ]
        other_list = np.repeat("Other", dfOther.shape[0])

        dfClust = pd.concat([dfClust, dfOther]).reset_index(drop=True)
        dfClust["Cell Type"] = np.concatenate([clust_list, other_list])

    sns.boxplot(
        data=dfClust,
        x="Cell Type",
        y="Average Gene Expression",
        hue="Condition",
        ax=ax,
        showfliers=False,
    )

    ax.set(
        title=gene,
        yticks=np.linspace(
            0, np.max(dfClust["Average Gene Expression"]) + 0.00005, num=5
        ),
    )
