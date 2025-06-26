import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from matplotlib.axes import Axes
from tensorly.cp_tensor import CPTensor
from tlviz.factor_tools import factor_match_score as fms

from ...factorization import pf2, pf2_pca_r2x


def plot_r2x(data, rank_vec, ax: Axes):
    """Creates R2X plot for parafac2 tensor decomposition and pca"""

    r2xError = pf2_pca_r2x(data, rank_vec)
    labelNames = ["Fit: Pf2", "Fit: PCA"]
    colorDecomp = ["r", "b"]
    markerShape = ["o", "o"]
    for i in range(2):
        ax.scatter(
            rank_vec,
            r2xError[i],
            label=labelNames[i],
            marker=markerShape[i],
            c=colorDecomp[i],
            s=30.0,
        )
    ax.set(
        ylabel="Variance Explained",
        xlabel="Number of Components",
        xticks=np.linspace(0, rank_vec[-1], num=6, dtype=int),
        yticks=np.linspace(
            0, np.max(np.append(r2xError[0], r2xError[1])) + 0.01, num=5
        ),
    )
    ax.legend()


def plot_avegene_per_celltype(adata, genes, ax, cellType="Cell Type"):
    """Plots average gene expression across cell types for all conditions"""
    genesV = adata[:, genes]
    dataDF = genesV.to_df()
    dataDF = dataDF.subtract(genesV.var["means"].values)
    dataDF["Condition"] = genesV.obs["Condition"].values
    dataDF["Cell Type"] = genesV.obs[cellType].values

    data = pd.melt(dataDF, id_vars=["Condition", "Cell Type"], value_vars=genes).rename(
        columns={"variable": "Gene", "value": "Value"}
    )
    df = data.groupby(["Condition", "Cell Type", "Gene"], observed=False).mean()
    df = df.rename(columns={"Value": "Average Gene Expression"})
    sns.boxplot(
        data=df,
        x="Gene",
        y="Average Gene Expression",
        hue="Cell Type",
        ax=ax,
        fliersize=0,
    )


def plot_avegene_per_category(conds, gene, adata, ax, mean=True, cellType="Cell Type"):
    """Plots average gene expression across cell types for a category of drugs"""
    genesV = adata[:, gene]
    dataDF = genesV.to_df()
    dataDF = dataDF.subtract(genesV.var["means"].values)
    dataDF["Condition"] = genesV.obs["Condition"].values
    dataDF["Cell Type"] = genesV.obs[cellType].values

    df = pd.melt(dataDF, id_vars=["Condition", "Cell Type"], value_vars=gene).rename(
        columns={"variable": "Gene", "value": "Value"}
    )
    if mean is True:
        df = df.groupby(["Condition", "Cell Type", "Gene"], observed=False).mean()

    df = df.rename(columns={"Value": "Average Gene Expression For Drugs"}).reset_index()
    df = df[df["Condition"].isin(conds)]

    # df["Condition"] = np.where(df["Condition"].isin(conds), df["Condition"], "Other")
    # df["Condition"] = df[df["Condition"]==conds]
    # for i in conds:
    # df = df.replace({"Condition": {i: categoryCond}})

    sns.boxplot(
        data=df.loc[df["Gene"] == gene],
        x="Cell Type",
        y="Average Gene Expression For Drugs",
        hue="Condition",
        ax=ax,
        showfliers=False,
    )
    ax.set(title=gene)
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=45)


def heatmapGeneFactors(
    cmps: list, dataIn: anndata.AnnData, ax: Axes, geneAmount: int = 20
):
    """Plotting weights for gene factors for both most negatively/positively weighted terms"""
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    df = pd.DataFrame(
        data=dataIn.varm["Pf2_C"],
        index=dataIn.var_names,
        columns=range(1, dataIn.varm["Pf2_C"].shape[1] + 1),
    )
    df = df.reset_index(names="Gene")

    genes = np.array([])
    for cmp in cmps:
        sortDF = df.sort_values(by=cmp)
        top = sortDF.iloc[-geneAmount:, :].Gene.values
        bottom = sortDF.iloc[:geneAmount:, :].Gene.values
        genes = np.concatenate((genes, np.flip(top)))
        genes = np.concatenate((genes, bottom))

    heatmapDF = df.loc[df.Gene.isin(genes)][cmps + ["Gene"]].set_index("Gene")
    vmax = np.abs(heatmapDF.values).max()

    sns.heatmap(
        data=heatmapDF.transpose()[genes], ax=ax, cmap=cmap, vmin=-vmax, vmax=vmax
    )


def cell_comp_hist(X, category: str, comp: int, unique, ax: Axes):
    """Plots weighted projections of each cell according to category"""
    w_proj = X.obsm["weighted_projections"][:, comp - 1]
    obsDF = pd.DataFrame({category: X.obs[category].astype(str).values})
    w_proj = X.obsm["weighted_projections"][:, comp - 1]
    if category is not None:
        if unique is not None:
            obsDF.loc[obsDF[category] != unique, category] = "Other"
        labels = obsDF[category]
        # obsDF[category] = obsDF[category].astype(str)
        histDF = pd.DataFrame({"Component " + str(comp): w_proj, category: labels})
        sns.histplot(
            data=histDF,
            x="Component " + str(comp),
            hue=category,
            kde=True,
            ax=ax,
            stat="density",
            common_norm=False,
        )


def gene_plot_cells(
    X: anndata.AnnData,
    hue: str,
    ax: Axes,
    unique=None,
    average=False,
    kde=False,
    cellType="Cell Type",
):
    """Plots two genes on either a per cell or per cell type basis"""
    assert X.shape[1] == 2
    genes = X.var_names
    dataDF = X.to_df()
    dataDF = dataDF.subtract(X.var["means"].values)
    dataDF[hue] = X.obs[hue].values
    dataDF["Cell Type"] = X.obs[cellType].values
    alpha = 1

    if average:
        dataDF = dataDF.groupby([hue], observed=True).mean().reset_index()
        alpha = 1

    if unique is not None:
        dataDF[hue] = dataDF[hue].astype(str)
        dataDF.loc[~dataDF[hue].isin(unique), hue] = "Other"

    sns.scatterplot(data=dataDF, x=genes[0], y=genes[1], hue=hue, ax=ax, alpha=alpha)
    if kde:
        sns.kdeplot(
            data=dataDF,
            x=genes[0],
            y=genes[1],
            hue=hue,
            levels=5,
            fill=True,
            alpha=0.3,
            cut=2,
            ax=ax,
        )


def gene_plot_conditions(X, condition: str, genes, ax: Axes, hue=None, unique=None):
    """Plots two genes on either a per cell or per cell type basis"""
    adata = X[:, [genes[0], genes[1]]]
    sc.pp.subsample(adata, fraction=0.01, random_state=0)

    dataDF = pd.DataFrame(columns=genes, data=adata.X)
    dataDF[condition] = adata.obs[condition].values
    dataDF[condition] = dataDF[condition].astype("str")
    if hue:
        dataDF[hue] = adata.obs[hue].values
        dataDF[hue] = dataDF[hue].astype("str")
        dataDF = dataDF.groupby([condition, hue]).mean()
    else:
        dataDF = dataDF.groupby([condition]).mean()
    if unique is not None:
        dataDF[condition] = dataDF[condition].astype(str)
        dataDF.loc[dataDF[condition] != unique, condition] = "Other"
    if hue is not None:
        sns.scatterplot(data=dataDF, x=genes[0], y=genes[1], hue=hue, ax=ax, alpha=5)
    else:
        sns.scatterplot(data=dataDF, x=genes[0], y=genes[1], ax=ax, alpha=0.2)


def geneSig_plot_cells(
    X, comps: list[int], hue: str, ax: Axes, unique=None, average=False, kde=False
):
    """Plots two genes on either a per cell or per cell type basis"""

    geneSigDF = pd.DataFrame()
    geneVecs = X.varm["Pf2_C"][:, comps]
    for i, _ in enumerate(comps):
        geneSigDF[str(comps[i])] = np.matmul(X.X, geneVecs[:, i])

    geneSigDF[hue] = X.obs[hue].values
    geneSigDF = geneSigDF.sample(n=10000)
    alpha = 0.3

    if unique is not None:
        geneSigDF[hue] = geneSigDF[hue].astype(str)
        geneSigDF.loc[geneSigDF[hue] != unique, hue] = "Other"
    if average:
        geneSigDF = geneSigDF.groupby([hue], observed=True).mean()
        alpha = 1
    sns.scatterplot(
        data=geneSigDF,
        x=str(comps[0]),
        y=str(comps[1]),
        hue=hue,
        ax=ax,
        alpha=alpha,
    )
    if kde:
        sns.kdeplot(
            data=geneSigDF,
            x=str(comps[0]),
            y=str(comps[1]),
            hue=hue,
            levels=5,
            fill=True,
            alpha=0.3,
            cut=2,
            ax=ax,
        )
    ax.set(
        xlabel="Comp. " + str(comps[0]) + " Signature",
        ylabel="Comp. " + str(comps[1]) + " Signature",
    )


def plot_cell_gene_corr(
    X: anndata.AnnData,
    hue: str,
    cells: list,
    ax: Axes,
    unique=None,
    cellType="Cell Type",
):
    """Plots two genes on either a per cell or per cell type basis"""
    assert X.shape[1] == 2
    genes = X.var_names
    dataDF = X.to_df()
    dataDF = dataDF.subtract(X.var["means"].values)
    dataDF[hue] = X.obs[hue].values
    dataDF["Cell Type"] = X.obs[cellType].values
    alpha = 0.3

    dataDF = dataDF.groupby([hue, "Cell Type"], observed=True).mean().reset_index()
    alpha = 1

    corrDF = pd.DataFrame()
    for cond in dataDF[hue].unique():
        cell_gene1 = dataDF.loc[
            (dataDF[hue] == cond) & (dataDF["Cell Type"] == cells[0])
        ][genes[0]].values
        cell_gene2 = dataDF.loc[
            (dataDF[hue] == cond) & (dataDF["Cell Type"] == cells[1])
        ][genes[1]].values
        corrDF = pd.concat(
            [
                corrDF,
                pd.DataFrame(
                    {
                        hue: cond,
                        cells[0] + " " + genes[0]: cell_gene1,
                        cells[1] + " " + genes[1]: cell_gene2,
                    }
                ),
            ]
        )

    if unique is not None:
        corrDF[hue] = corrDF[hue].astype(str)
        corrDF.loc[~corrDF[hue].isin(unique), hue] = "Other"

    sns.scatterplot(
        data=corrDF,
        x=cells[0] + " " + genes[0],
        y=cells[1] + " " + genes[1],
        hue=hue,
        ax=ax,
        alpha=alpha,
    )


def cell_count_perc_df(X, celltype="Cell Type", status=False, grouping="Condition"):
    """Returns DF with cell counts and percentages for experiment"""
    if status is False:
        grouping = [celltype, grouping]
    else:
        grouping = [celltype, "Condition", "SLE_status"]

    df = X.obs[grouping].reset_index(drop=True)

    dfCond = df.groupby([grouping], observed=True).size().reset_index(name="Cell Count")
    dfCellType = (
        df.groupby(grouping, observed=True).size().reset_index(name="Cell Count")
    )
    dfCellType["Cell Count"] = dfCellType["Cell Count"].astype("float")

    dfCellType["Cell Type Percentage"] = 0.0
    for cond in np.unique(df[grouping]):
        dfCellType.loc[dfCellType[grouping] == cond, "Cell Type Percentage"] = (
            100
            * dfCellType.loc[dfCellType[grouping] == cond, "Cell Count"].to_numpy()
            / dfCond.loc[dfCond[grouping] == cond]["Cell Count"].to_numpy()
        )

    dfCellType.rename(columns={celltype: "Cell Type"}, inplace=True)

    return dfCellType


def plot_gene_set_expression(adata, gene_set, ax: Axes):
    """
    Plots the average gene expression level for a given gene set per condition.

    Parameters:
    adata (anndata.AnnData): The AnnData object containing the data.
    gene_set (list): A list of genes to include in the gene set.
    ax (Axes): The matplotlib axes to plot on.
    cellType (str): The cell type to group by (default is "Cell Type").
    """
    # Filter the genes in the gene set that are present in the data
    valid_genes = [gene for gene in gene_set if gene in adata.var_names]
    if not valid_genes:
        raise ValueError("None of the genes in the gene set are present in the data.")

    # Extract the data for the valid genes
    genesV = adata[:, valid_genes]
    dataDF = genesV.to_df()
    dataDF = dataDF.subtract(genesV.var["means"].values)
    dataDF["Condition"] = genesV.obs["Condition"].values

    # Calculate the average expression for the gene set
    dataDF["Average Gene Expression"] = dataDF[valid_genes].mean(axis=1)

    # Create the plot
    sns.boxplot(
        data=dataDF,
        x="Condition",
        y="Average Gene Expression",
        ax=ax,
        showfliers=False,
    )
    ax.set_title("Average Gene Expression for Gene Set")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Average Gene Expression")
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=45)


def rotate_xaxis(ax, rotation=90):
    """Rotates text by 90 degrees for x-axis"""
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=rotation)


def rotate_yaxis(ax, rotation=90):
    """Rotates text by 90 degrees for y-axis"""
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(labels=ax.get_yticklabels(), rotation=rotation)


def plot_boxplot_gene_celltype(
    conds, gene, adata, ax, mean=False, cellType="Cell Type", cells=["T reg"]
):
    """Boxplot of gene expression for a specific cell type across conditions"""
    grouping = [cellType, "Condition"]

    df = adata.obs[grouping].reset_index(drop=True)
    grouped_df = (
        adata.obs.groupby(["CellType2", "Condition"], observed=False)
        .size()
        .reset_index(name="Cell Count")
    )

    df = df[df["Condition"].isin([conds[0]])]
    df = df[df["CellType2"].isin(cells)]

    # print(df)
    genesV = adata[:, gene]
    if scipy.sparse.issparse(genesV.X):
        # If the data is sparse, convert to dense array first
        expr_values = genesV.X.toarray()
    else:
        expr_values = genesV.X
    gene_mean = genesV.var["means"].values[0]  # Index 0 since it's a single gene
    print(gene_mean)
    # Subtract mean from expression values
    expr_values = expr_values - gene_mean
    dataDF = pd.DataFrame(expr_values, columns=[gene])
    dataDF["Condition"] = genesV.obs["Condition"].values
    dataDF["Cell Type"] = genesV.obs["CellType2"].values

    df = dataDF[dataDF["Condition"].isin(conds)]
    df = df[df["Cell Type"].isin(cells)]

    if mean is True:
        df = df.groupby(["Condition", "Cell Type"], observed=False).mean()

    df.rename(columns={gene: "Gene Expression"}, inplace=True)

    print(df)

    sns.boxplot(
        data=df,
        x="Condition",
        y="Gene Expression",
        hue="Condition",
        ax=ax,
        showfliers=False,
    )
    ax.set(title=gene)
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=45)


def calculateFMS(A: anndata.AnnData, B: anndata.AnnData):
    """Calculates FMS between 2 factors"""
    factors = [A.uns["Pf2_A"], A.uns["Pf2_B"], A.varm["Pf2_C"]]
    A_CP = CPTensor(
        (
            A.uns["Pf2_weights"],
            factors,
        )
    )

    factors = [B.uns["Pf2_A"], B.uns["Pf2_B"], B.varm["Pf2_C"]]
    B_CP = CPTensor(
        (
            B.uns["Pf2_weights"],
            factors,
        )
    )

    return fms(A_CP, B_CP, consider_weights=False, skip_mode=1)  # type: ignore


def plot_fms_percent_drop(
    X: anndata.AnnData,
    ax: Axes,
    percentList: np.ndarray,
    runs: int,
    rank: int = 30,
):
    # Plots FMS score when percentage is removed from data
    dataX = pf2(X, rank, doEmbedding=False)

    fmsLists = []

    for j in range(0, runs, 1):
        scores = [1.0]

        for i in percentList[1:]:
            sampled_data: anndata.AnnData = sc.pp.subsample(
                X, fraction=1 - (i / 100), random_state=j, copy=True
            )  # type: ignore
            sampledX = pf2(sampled_data, rank, random_state=j + 2, doEmbedding=False)

            fmsScore = calculateFMS(dataX, sampledX)
            scores.append(fmsScore)

        fmsLists.append(scores)

    runsList_df = []
    for i in range(0, runs):
        for j in range(0, len(percentList)):
            runsList_df.append(i)
    percentList_df = []
    for i in range(0, runs):
        for j in range(0, len(percentList)):
            percentList_df.append(percentList[j])
    fmsList_df = []
    for sublist in fmsLists:
        fmsList_df += sublist
    df = pd.DataFrame(
        {
            "Run": runsList_df,
            "Percentage of Data Dropped": percentList_df,
            "FMS": fmsList_df,
        }
    )

    sns.lineplot(data=df, x="Percentage of Data Dropped", y="FMS", ax=ax)
    ax.set_ylim(0, 1)


def resample(data: anndata.AnnData) -> anndata.AnnData:
    """Bootstrapping dataset"""
    indices = np.random.randint(0, data.shape[0], size=(data.shape[0],))
    data = data[indices].copy()
    return data


def plot_fms_diff_ranks(
    X: anndata.AnnData,
    ax: Axes,
    ranksList: list[int],
    runs: int,
):
    # Plots FMS when using different Pf2 components
    fmsLists = []

    for j in range(0, runs, 1):
        scores = []
        for i in ranksList:
            dataX = pf2(X, rank=i, random_state=j, doEmbedding=False)

            sampledX = pf2(resample(X), rank=i, random_state=j, doEmbedding=False)

            fmsScore = calculateFMS(dataX, sampledX)
            scores.append(fmsScore)
        fmsLists.append(scores)

    runsList_df = []
    for i in range(0, runs):
        for j in range(0, len(ranksList)):
            runsList_df.append(i)
    ranksList_df = []
    for i in range(0, runs):
        for j in range(0, len(ranksList)):
            ranksList_df.append(ranksList[j])
    fmsList_df = []
    for sublist in fmsLists:
        fmsList_df += sublist
    df = pd.DataFrame(
        {"Run": runsList_df, "Component": ranksList_df, "FMS": fmsList_df}
    )

    sns.lineplot(data=df, x="Component", y="FMS", ax=ax)
    ax.set_ylim(0, 1)
