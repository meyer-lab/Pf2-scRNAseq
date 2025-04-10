"""
Plots ave gene expression, plots each cell
"""

import anndata
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.axes import Axes

from ..imports import import_pf2Cytokine30, import_cytokine
from ..factorization import pf2
from anndata import read_h5ad
from .common import getSetup, subplotLabel
from .commonFuncs.plotFactors import plot_gene_factors
from .commonFuncs.plotGeneral import cell_count_perc_df, rotate_xaxis
from .commonFuncs.plotPaCMAP import plot_gene_pacmap, plot_wp_pacmap



def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((8, 8), (1, 1))
    subplotLabel(ax)

    
    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")

    plot_gene_expression(X, "FOXP3", ax[0], cellType="CellType2", cell_types=["CD4+ T cells", "T reg"], jitter=False)
    

    #plot_gene_pacmap("RETN", "Pf2", X, ax[4])

    #celltype_count_perc_df = cell_count_perc_df(X, celltype="leiden", status=True)
    
    return f


def plot_gene_expression(X: anndata.AnnData, gene: str, ax: Axes, cellType: str = "CellType2", 
                        cell_types=None, conditions=None, jitter=True):
    """
    Plots individual cell gene expression values with average bars.
    
    Parameters:
    -----------
    X: anndata.AnnData
        AnnData object containing gene expression data
    gene: str
        Name of the gene to plot
    ax: Axes
        Matplotlib axes to plot on
    cellType: str
        Name of the column in X.obs containing cell type information
    cell_types: list, optional
        List of cell types to include (if None, include all)
    conditions: list, optional
        List of conditions to include (if None, include all)
    jitter: bool, optional
        Whether to jitter points (True) or use swarm plot (False)
    """
    # Make gene names unique if needed
    if len(X.var_names) != len(set(X.var_names)):
        X = X.copy()
        X.var_names_make_unique()
    
    # Extract gene expression values
    gene_data = X[:, gene]
    
    # Convert to array if sparse
    if scipy.sparse.issparse(gene_data.X):
        expr_values = gene_data.X.toarray().flatten()
    else:
        expr_values = gene_data.X.flatten()
    
    # Create DataFrame with gene expression values
    df = pd.DataFrame({
        "Expression": expr_values,
        "Condition": gene_data.obs["Condition"].values,
        "Cell Type": gene_data.obs[cellType].values
    })
    
    # Filter for specific cell types and conditions if provided
    if cell_types is not None:
        df = df[df["Cell Type"].isin(cell_types)]
    
    if conditions is not None:
        df = df[df["Condition"].isin(conditions)]

    # Define cell type order
    cell_type_order = sorted(df["Cell Type"].unique())
        
    # Plot individual points
    if jitter:
        # Use stripplot for individual points with jitter
        sns.stripplot(
            data=df,
            x="Cell Type",
            y="Expression",
            hue="Condition",
            dodge=True,
            order=cell_type_order,
            alpha=0.4,
            size=3,
            jitter=0.3,
            ax=ax
        )
    else:
        # Use swarmplot for individual points (no jitter)
        sns.swarmplot(
            data=df,
            x="Cell Type",
            y="Expression",
            hue="Condition",
            dodge=True,
            order=cell_type_order,
            alpha=0.4,
            size=3,
            ax=ax
        )
    
    # Add bar for the average (using pointplot with line segments)
    sns.pointplot(
        data=df,
        x="Cell Type",
        y="Expression",
        hue="Condition",
        order=cell_type_order,
        dodge=0.5,
        join=False,
        scale=1.2,
        ci=None,
        markers="_",
        palette="dark",
        ax=ax
    )
    
    # Fix legend (only keep one set of entries)
    handles, labels = ax.get_legend_handles_labels()
    n_conditions = len(df["Condition"].unique())
    ax.legend(handles[:n_conditions], 
              labels[:n_conditions],
              title="Condition",
              bbox_to_anchor=(1.05, 1),
              loc="upper left")
    
    # Set title and labels
    ax.set_title(f"{gene} Expression")
    ax.set_ylabel(f"Expression Level")
    
    # Rotate x-axis labels for readability
    rotate_xaxis(ax)
    
    return df
