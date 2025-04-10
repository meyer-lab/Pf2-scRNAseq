"""
Cell type percentage per Leiden cluster per condition
"""

from anndata import read_h5ad
from ..imports import import_cytokine,import_pf2Cytokine30 
from .common import subplotLabel, getSetup
import seaborn as sns
from .commonFuncs.plotGeneral import cell_count_perc_df


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((20, 3), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    #X = import_pf2Cytokine30 ()
    #X=read_h5ad("/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Cytokine_Pf2_annotated_NB_new.h5ad")
    X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")

    
    df = cell_count_perc_df(X, celltype="CellType2")

    sns.barplot(
        data=df,
        x="Cell Type",
        y="Cell Type Percentage",
        hue="Condition",
        ax=ax[0],
        errorbar=None,
    )

    return f
