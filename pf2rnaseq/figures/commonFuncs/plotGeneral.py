import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import seaborn as sns
from matplotlib.axes import Axes

from ...factorization import fms_diff_ranks, fms_percent_drop, pf2_pca_r2x


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
    # Increase font sizes
    ax.set_xlabel("Number of Components", fontsize=18)
    ax.set_ylabel("Variance Explained", fontsize=18)
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.legend()


def plot_avegene_per_celltype(
    adata,
    genes,
    ax,
    cellType="Cell Type",
    condition="cytokine",
    mean=False,
    center_data=False,
):
    """Plots average gene expression across cell types for all conditions"""
    genesV = adata[:, genes]
    dataDF = genesV.to_df()
    if center_data:
        dataDF = dataDF.subtract(genesV.var["means"].values)
        gene_std = dataDF.std(axis=0)  # Standard deviation for each gene
        # Z-score: data is already mean-centered, just divide by std
        dataDF = dataDF / gene_std
    dataDF["Condition"] = genesV.obs[condition].values
    dataDF["Cell Type"] = genesV.obs[cellType].values

    # Melt the data
    data = pd.melt(dataDF, id_vars=["Condition", "Cell Type"], value_vars=genes).rename(
        columns={"variable": "Gene", "value": "Gene Expression"}
    )

    # Apply grouping if mean=True
    if mean is True:
        data = (
            data.groupby(["Condition", "Cell Type", "Gene"], observed=False)
            .mean()
            .reset_index()
        )

    sns.boxplot(
        data=data,
        x="Cell Type",
        y="Gene Expression",
        hue="Condition",
        ax=ax,
        fliersize=0,
    )
    ax.set_title(genes, fontsize=25)
    ax.set_xlabel("Cell Type", fontsize=22)
    ax.set_ylabel("Gene Expression", fontsize=22)
    ax.tick_params(axis="x", rotation=45, labelsize=20)
    ax.tick_params(axis="y", labelsize=20)


def plot_avegene_per_category(
    conds,
    gene,
    adata,
    ax,
    mean=True,
    cellType="Cell Type",
    swarm=False,
    center_data=False,
    condition="cytokine",
):
    """Plots average gene expression across cell types for a specified condition"""
    genesV = adata[:, gene]
    dataDF = genesV.to_df()

    if center_data:
        gene_means = dataDF.mean(axis=0)
        dataDF = dataDF - gene_means
        gene_std = dataDF.std(axis=0)  # Standard deviation for each gene
        # Z-score: data is already mean-centered, just divide by std
        dataDF = dataDF / gene_std

    dataDF["Condition"] = genesV.obs[condition].values
    dataDF["Cell Type"] = genesV.obs[cellType].values

    df = pd.melt(dataDF, id_vars=["Condition", "Cell Type"], value_vars=[gene]).rename(
        columns={"variable": "Gene", "value": "Gene Expression"}
    )

    if mean is True:
        df = (
            df.groupby(["Condition", "Cell Type", "Gene"], observed=False)
            .mean()
            .reset_index()
        )

    df["Condition"] = np.where(df["Condition"].isin(conds), df["Condition"], "Other")

    if swarm is False:
        sns.boxplot(
            data=df,
            x="Cell Type",
            y="Gene Expression",
            hue="Condition",
            ax=ax,
            showfliers=False,
        )
    else:
        sns.stripplot(
            data=df,
            x="Condition",
            y="Gene Expression",
            hue="Condition",
            ax=ax,
            alpha=0.6,
        )

    ax.set(title=gene)
    ax.tick_params(axis="x", rotation=45)


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


def cell_count_perc_df(X, celltype="Cell Type", condition="cytokine"):
    """Returns DF with cell counts and percentages for experiment"""
    grouping_all = [
        celltype,
        condition,
        "condition_unique_idxs",
    ]
    grouping = [celltype, condition]

    df = X.obs[grouping_all].reset_index(drop=True)

    idx_mapping = X.obs.groupby(condition, observed=False)[
        "condition_unique_idxs"
    ].first()

    dfCond = (
        df.groupby([condition], observed=True).size().reset_index(name="Cell Count")
    )
    dfCellType = (
        df.groupby(grouping, observed=True).size().reset_index(name="Cell Count")
    )
    dfCellType["Cell Count"] = dfCellType["Cell Count"].astype("float")

    dfCellType["Cell Type Percentage"] = 0.0
    for cond in np.unique(df[condition]):
        dfCellType.loc[dfCellType[condition] == cond, "Cell Type Percentage"] = (
            100
            * dfCellType.loc[dfCellType[condition] == cond, "Cell Count"].to_numpy()
            / dfCond.loc[dfCond[condition] == cond]["Cell Count"].to_numpy()
        )

    dfCellType["condition_unique_idxs"] = dfCellType[condition].map(idx_mapping)
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


def plot_fms_diff_ranks(
    X: anndata.AnnData,
    ax: Axes,
    ranksList: list[int],
    runs=3,
):
    """Plots FMS when using different Pf2 components"""
    df = fms_diff_ranks(X, ranksList, runs)
    sns.lineplot(data=df, x="Component", y="FMS", ax=ax)
    ax.set_ylim(0, 1)


def plot_fms_percent_drop(
    X: anndata.AnnData, ax: Axes, percentList: np.ndarray, runs=3, rank: int = 30
):
    """Plots FMS when dropping different percentages of data"""
    df = fms_percent_drop(X, percentList, runs, rank)
    sns.lineplot(data=df, x="Percentage of Data Dropped", y="FMS", ax=ax)
    ax.set_ylim(0, 1)
