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






def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((30, 8), (5, 5))

    # Add subplot labels
    subplotLabel(ax)
    X = read_h5ad("/opt/andrew/lupus/lupus_fitted_ann.h5ad")
    bulk_matrix, bulk_tensor = pseudobulk_lupus(X)
    print(bulk_matrix )
   
    print(np.shape(bulk_matrix))

    bulk_matrix_genes_only = bulk_matrix.iloc[0:, 2:-1]
    #print(Bulk_matrix_genes_only)

    pca = PCA(n_components=4)
    scores = pca.fit(bulk_matrix_genes_only).transform(bulk_matrix_genes_only)
    loadings = pca.components_.T

    print(np.shape(scores))
    print(np.shape(loadings))

    scores = scores[:, [0, 1, 2, 3,]]
    laodings = loadings[:, [0, 1, 2, 3,]]

    scores_df = pd.DataFrame(date = scores, columns = ["PC1", "PC2", "PC3", "PC4"])
    loadings_df = pd.DataFrame(date = laodings, columns = ["PC1", "PC2", "PC3", "PC4"])

    loadings_df["Gene"] = bulk_matrix_genes_only.columns.to_numpy()
    print(loadings_df)

    scores_df["Status"] = bulk_matrix["Status"].to_numpy()
    scores_df["Conditon"] = bulk_matrix["Condition"].to_numpy()
    scores_df["Cell Type"] = bulk_matrix["Cell Type"].to_numpy()
    loadings_df["Gene"] = bulk_matrix_genes_only.columns.to_numpy()

    sns.scatterplot(data=scores_df, x="PC1", y="PC2", ax=ax[0])
    sns.scatterplot(data=scores_df, x="PC1", y="PC3", hue="Status", ax=ax[1])
    sns.scatterplot(data=scores_df, x="PC1", y="PC4", hue="Condition", ax=ax[2])
    sns.scatterplot(data=scores_df, x="PC2", y="PC1", hue="Cell Type", ax=ax[3])
    sns.scatterplot(data=scores_df, x="PC2", y="PC3", ax=ax[4])
    sns.scatterplot(data=scores_df, x="PC2", y="PC4", hue="Status", ax=ax[5])
    sns.scatterplot(data=scores_df, x="PC3", y="PC1", hue="Condition", ax=ax[6])
    sns.scatterplot(data=scores_df, x="PC3", y="PC2", hue="Cell Type", ax=ax[7])
    sns.scatterplot(data=scores_df, x="PC3", y="PC4", ax=ax[8])
    sns.scatterplot(data=scores_df, x="PC4", y="PC1", hue="Status", ax=ax[9])
    sns.scatterplot(data=scores_df, x="PC4", y="PC2", hue="Condition", ax=ax[10])
    sns.scatterplot(data=scores_df, x="PC4", y="PC3", hue="Cell Type", ax=ax[11])

    sns.scatterplot(date=loadings_df, x="PC1", y="PC2", ax=ax[12])
    sns.scatterplot(date=loadings_df, x="PC1", y="PC3", ax=ax[13])
    sns.scatterplot(date=loadings_df, x="PC1", y="PC4", ax=ax[14])
    sns.scatterplot(date=loadings_df, x="PC2", y="PC1", ax=ax[15])
    sns.scatterplot(date=loadings_df, x="PC2", y="PC3", ax=ax[16])
    sns.scatterplot(date=loadings_df, x="PC2", y="PC4", ax=ax[17])
    sns.scatterplot(date=loadings_df, x="PC3", y="PC1", ax=ax[18])
    sns.scatterplot(date=loadings_df, x="PC3", y="PC2", ax=ax[19])
    sns.scatterplot(date=loadings_df, x="PC3", y="PC4", ax=ax[20])
    sns.scatterplot(date=loadings_df, x="PC4", y="PC1", ax=ax[21])
    sns.scatterplot(date=loadings_df, x="PC4", y="PC2", ax=ax[22])
    sns.scatterplot(date=loadings_df, x="PC4", y="PC3", ax=ax[23])

    prop_var_df = pca.explained_variance_ratio_
    prop_var_df = np.cumsum(pca.explained_variance_ratio_)
    print(prop_var_df)
    PC_Numbers = np.arrange(pca.n_components_) + 1
    print(PC_Numbers)
    ax[24].scatter(x=PC_Numbers, y=prop_var_df)
    sns.lineplot(x=PC_Numbers, y=prop_var_df, ax=ax[24])
    plt.xlabel("PC Number")
    plt.ylabel("Variance")
    plt.show
    return f









