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
    ax, f = getSetup((30, 8), (2, 4))

    # Add subplot labels
    subplotLabel(ax)


    X = read_h5ad("/opt/andrew/lupus/lupus_fitted_ann.h5ad")
    bulk_matrix, bulk_tensor = pseudobulk_lupus(X)
    print(bulk_matrix )
    # bulk_matrix_genes_only = bulk_matrix.iloc[:, 2:-1]



    # pca = PCA(n_components=2)
    # principalComponents = pca.fit_transform(bulk_matrix_genes_only)
    # print(np.shape(principalComponents))

    # sns.scatterplot(data=data, x="PC1", y="PC2", hue="Status", ax=ax[0])


    # url = (bulk_matrix)

    # df = pd.read_csv(url, names=['Condtion length','Condition width','Cell length','Cell width', 'Type length', 'Cell Width', 'Target'])

    # features = ["Condition", "Cell", "Type"]

    # x = df.loc[:, features].values

    # y = df.loc[:, ['target']].values

    # x = StandardScaler().fit_transform(bulk_matrix)


    # pca = PCA(n_components=2)

    # principalComponents = pca.fit_transform(bulk_matrix)

    # principalDf  = pd.DataFrame(data = principalComponents, columns = ["principal component 1", "princpal component 2"])

    # finalDF = pd.concat([principalDf, df[["target"]]], axis = 1)

    # #fig - plt.figure(figsize = (8, 8))

    # ax = eig.add_subplot(1,1,1)

    # ax.set_xlable('Principal Component 1', frontsize = 15)

    # ax.set_ylable("Principal Component 2", fontsize = 15)

    # ax.set_title("2 component PCA", fontsize = 20)

    # targets = ["Condition", "Cell", "Type"]

    # colors = ["r", "g", "b"]
    # for target, color in zip(targets, colors):
    #     indiciesToKeep = finalDF["taregt"] == target 
    #     ff = makeFigure()
    #     ax.scatter (finalDF.loc[indiciesToKeep, "principal component 1"], finalDF.loc[indiciesToKeep, "principal component 2"], c = color, s = 50)

    # pca.explained_variance_ratio_
    # ax.legen(targets)
    # ax.grid()

    #print(bulk_matrix)
    return f








