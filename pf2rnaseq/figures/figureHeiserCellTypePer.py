"""
Cell type percentage per Leiden cluster per condition
"""

from anndata import read_h5ad
from ..imports import import_Heiser
from .common import subplotLabel, getSetup
import seaborn as sns
from .commonFuncs.plotGeneral import cell_count_perc_df


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    # Get list of axis objects
    ax, f = getSetup((20, 3), (1, 1))

    # Add subplot labels
    subplotLabel(ax)

    
   
    X = import_Heiser()

    
    df = cell_count_perc_df(X, celltype="cell_type_1", grouping="expBatch")

    sns.barplot(
        data=df,
        x="Cell Type",
        y="Cell Type Percentage",
        hue="expBatch",
        ax=ax[0],
        errorbar=None,
    )

    return f
