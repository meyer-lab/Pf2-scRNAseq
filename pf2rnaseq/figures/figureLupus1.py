"""
XX:
"""

from matplotlib.pylab import eig
import plotly.express as px
from anndata import read_h5ad
from ..imports import pseudobulk_lupus
from .common import subplotLabel, getSetup
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from os.path import join, dirname
from tensorly.decomposition import parafac
from tensorly.cp_tensor import CPTensor
from tensorly.tucker_tensor import TuckerTensor, tucker_to_tensor
from tensorly.tenalg import multi_mode_dot

def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    # Makes the grpahs, First order pair deals with sizing of the actual graph and second order pair makes the amount of rows and colomuns of graphs.
    ax, f = getSetup((20, 20), (4, 4)) 
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad("/opt/andrew/lupus/lupus_fitted_ann.h5ad")

    bulk_matrix, bulk_tensor = pseudobulk_lupus(X)
    # cell_type_names = pd.unique(bulk_matrix["Cell Type"])
    # #genes = pd.unique(bulk_matrix[])
    # condition_type = pd.unique(bulk_matrix["Condition"])
    # status_type = pd.unique(bulk_matrix["Status"])

    # print(bulk_matrix)

    # cp = parafac(bulk_tensor, 3)
    # print(cp[1])
    # sns.heatmap(data=cp[1][1], ax=ax[0],
    # xticklabels=[str(ii+1) for ii in range(cp[1][1].shape[1])], yticklabels=cell_type_names)
    # ax[0].set_xlabel("Components")
    # ax[0].set_ylabel("Cell Type")
    # sns.heatmap(data=cp[1][1], ax=ax[1],
    # xticklabels=[str(ii+1) for ii in range(cp[1][1].shape[1])], yticklabels="Condition")
    # ax[1].set_xlabel("Components")
    # ax[1].set_ylabel("Condition") 
    # sns.heatmap(data=cp[1][1], ax=ax[2],
    # xticklabels=[str(ii+1) for ii in range(cp[1][1].shape[1])], yticklabels=status_type)
    # ax[2].set_xlabel("Components")
    # ax[2].set_ylabel("Status")
    # sns.heatmap(data=cp[1][1], ax=ax[3],
    # xticklabels=[str(ii+1) for ii in range(cp[1][1].shape[1])], yticklabels= "genes")
    # ax[3].set_xlabel("Components")
    # ax[3].set_ylabel("genes")

    # errors = pd.DataFrame(columns=["Num of Comps.", "Error"])
    # for bulk_matirx in np.arange(1, 2):
    #     errors.loc[len(errors)+1] = [3, errors(bulk_tensor, parafac(bulk_tensor, 3))]

    # sns.pointplot(data=errors, x="Num of Comps.", y="Error" , ax=ax[4])






    #print(np.shape(bulk_matrix))

    # bulk_matrix = bulk_matrix.iloc[0:, 2:2000]
    # #print(bulk_matrix)
    # bulk_matrix = bulk_matrix.iloc[0:, 2:2000]

    bulk_matrix_genes_only = bulk_matrix.iloc[:, 2:-1]
    print(bulk_matrix_genes_only)

    pca = PCA(n_components=4)
    scores = pca.fit(bulk_matrix_genes_only).transform(bulk_matrix_genes_only)
    loadings = pca.components_.T

    print(np.shape(scores))
    print(np.shape(loadings))

    scores = scores[:, [0, 1, 2, 3]]

    loadings = loadings[:, [0, 1, 2, 3]]

    scores_df = pd.DataFrame(data = scores, columns = ["PC1", "PC2", "PC3", "PC4"])
    loadings_df = pd.DataFrame(data = loadings, columns = ["PC1", "PC2", "PC3", "PC4"])

    loadings_df["Gene"] = bulk_matrix_genes_only.columns.to_numpy()

    print(loadings_df)

    scores_df["Status"] = bulk_matrix["Status"].to_numpy()
    scores_df["Condition"] = bulk_matrix["Condition"].to_numpy()
    scores_df["Cell Type"] = bulk_matrix["Cell Type"].to_numpy()
    loadings_df["Gene"] = bulk_matrix_genes_only.columns.to_numpy()

    sns.scatterplot(data=scores_df, x="PC1", y="PC2", hue="Cell Type", ax=ax[0]).set_title("Scores")
    sns.scatterplot(data=scores_df, x="PC1", y="PC3", hue="Status", ax=ax[1]).set_title("Scores")
    sns.scatterplot(data=scores_df, x="PC1", y="PC4", ax=ax[2]).set_title("Scores")
    sns.scatterplot(data=scores_df, x="PC2", y="PC3", hue="Cell Type", ax=ax[3]).set_title("Scores")
    sns.scatterplot(data=scores_df, x="PC2", y="PC4", hue="Status", ax=ax[4]).set_title("Scores")
    sns.scatterplot(data=scores_df, x="PC3", y="PC4", ax=ax[5]).set_title("Scores")


    sns.scatterplot(data=loadings_df, x="PC1", y="PC2",ax=ax[6]).set_title("Loadings")
    sns.scatterplot(data=loadings_df, x="PC1", y="PC3",ax=ax[7]).set_title("Loadings")
    sns.scatterplot(data=loadings_df, x="PC1", y="PC4",ax=ax[8]).set_title("Loadings")
    sns.scatterplot(data=loadings_df, x="PC2", y="PC3",ax=ax[9]).set_title("Loadings")
    sns.scatterplot(data=loadings_df, x="PC2", y="PC4",ax=ax[10]).set_title("Loadings")
    sns.scatterplot(data=loadings_df, x="PC3", y="PC4",ax=ax[11]).set_title("Loadings")


    prop_var_df = pca.explained_variance_ratio_
    prop_var_df  = np.cumsum(pca.explained_variance_ratio_)
    print(prop_var_df)
    PC_Numbers = np.arange(pca.n_components_) + 1
    print(PC_Numbers)
    ax[12].scatter(x=PC_Numbers, y=prop_var_df)
    sns.lineplot(x=PC_Numbers, y=prop_var_df, ax=ax[12])
    plt.xlabel("PC Numbers")
    plt.ylabel("Variance")
    plt.title("R2X")
    plt.show
    return f
