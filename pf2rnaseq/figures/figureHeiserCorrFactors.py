"""
Figure A8:
"""

import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.axes import Axes
import anndata
from .common import subplotLabel, getSetup

from .figureHeiserCorrelation import partial_correlation_matrix
# from ..data_import import condition_factors_meta
# from .commonFuncs.plotGeneral import bal_combine_bo_covid


def convert_to_patients(data: anndata.AnnData, sample: bool = False) -> pd.Series:
    """
    Converts unique IDs to patient IDs.

    Parameters:
        data (anndata.AnnData): single-cell measurements
        sample (bool): return sample ID conversion

    Returns:
        conversions (pd.Series): maps unique IDs to patient IDs.
    """

    conversions = data.obs.loc[:, ["sample_id", "condition_unique_idxs"]]

    conversions.set_index("condition_unique_idxs", inplace=True, drop=True)
    conversions = conversions.loc[~conversions.index.duplicated()]
    conversions.sort_index(ascending=True, inplace=True)
    conversions = conversions.squeeze()

    return conversions


def condition_factors_meta(X):
    """Extract condition factors and match with metadata already in X"""
    # Get condition factors from PARAFAC2 result
    condition_factors = X.uns["Pf2_A"]

    # Get unique conditions with metadata (one row per condition)
    unique_conditions = X.obs.drop_duplicates(subset="condition_unique_idxs")
    sample_ids = unique_conditions["sample_id"].values

    # Create a conversion mapping
    sample_conversions = convert_to_patients(X)

    # Create DataFrame with condition factors
    condition_factors_df = pd.DataFrame(
        index=sample_conversions.values,
        data=condition_factors,
        columns=[f"Cmp. {i}" for i in np.arange(1, condition_factors.shape[1] + 1)],
    )

    # Extract relevant metadata columns from X.obs
    metadata_cols = ["treatment", "expBatch"]  # Replace with your actual columns
    metadata = unique_conditions.set_index("sample_id")[metadata_cols].copy()

    # Ensure alignment
    metadata = metadata.loc[metadata.index.isin(sample_conversions.values)]
    metadata = metadata.reindex(sample_conversions.values).dropna(axis=0, how="all")

    # Combine factors with metadata
    merged_df = pd.concat([condition_factors_df, metadata], axis=1)

    return merged_df


def makeFigure():
    """Get a list of the axis objects and create a figure."""
    ax, f = getSetup((6, 12), (6, 3))
    subplotLabel(ax)

    X = anndata.read_h5ad("/home/nicoleb/C3TAg_Pf2_30.h5ad")

    X.uns["Pf2_A"] -= np.min(X.uns["Pf2_A"], axis=0)
    X.uns["Pf2_A"] += np.median(X.uns["Pf2_A"], axis=0)
    X.uns["Pf2_A"] = np.log(X.uns["Pf2_A"])

    cmp_columns = [f"Cmp. {i}" for i in np.arange(1, X.uns["Pf2_A"].shape[1] + 1)]

    cond_fact_meta_df = condition_factors_meta(X)

    pc_df = partial_correlation_matrix(cond_fact_meta_df[cmp_columns])
    pc_df = remove_low_pc_cmp(pc_df, abs_threshold=0.4)

    pc_df["Var1"] = pc_df["Var1"].map(lambda x: x.lstrip("Cmp. ")).astype(int)
    pc_df["Var2"] = pc_df["Var2"].map(lambda x: x.lstrip("Cmp. ")).astype(int)

    pc_abs_df = pc_df.copy()
    pc_abs_df["Weight"] = np.abs(pc_df["Weight"])
    pc_abs_df = pc_abs_df.sort_values("Weight")

    for i in range(6):
        # Get the component pair and correlation coefficient
        cmp1 = pc_abs_df.iloc[-(i + 1), 0]
        cmp2 = pc_abs_df.iloc[-(i + 1), 1]
        

        # Pass to plotting functions
        plot_pair_gene_factors(X, cmp1, cmp2, ax[(3 * i)])
        plot_pair_cond_factors(
            cond_fact_meta_df,
            cmp1,
            cmp2,
            ax[(3 * i) + 1],
            label="treatment",
           
        )
        plot_pair_wp(X, cmp1, cmp2, ax[(3 * i) + 2], frac=0.001)

    return f


def plot_pair_gene_factors(
    X: anndata.AnnData, cmp1: int, cmp2: int, ax: Axes
):
    """Plots two gene components weights"""
    cmpWeights = np.concatenate(
        ([X.varm["Pf2_C"][:, cmp1 - 1]], [X.varm["Pf2_C"][:, cmp2 - 1]])
    )
    df = pd.DataFrame(
        data=cmpWeights.transpose(), columns=[f"Cmp. {cmp1}", f"Cmp. {cmp2}"]
    )
    sns.scatterplot(data=df, x=f"Cmp. {cmp1}", y=f"Cmp. {cmp2}", ax=ax, color="k")
    ax.set(title=f"Gene Factors")


def plot_pair_cond_factors(
    df: pd.DataFrame, cmp1: int, cmp2: int, ax: Axes, label: str
):
    """Plots two condition components weights"""
    sns.scatterplot(data=df, x=f"Cmp. {cmp1}", y=f"Cmp. {cmp2}", hue=label, ax=ax)
    ax.set(title=f"Condition Factors")


def plot_pair_wp(
    X: anndata.AnnData,
    cmp1: int,
    cmp2: int,
    ax: Axes,
    frac: float = 0.001
    
):
    """Plots two weighted projections components weights"""
    cmpWeights = np.concatenate(
        (
            [X.obsm["weighted_projections"][:, cmp1 - 1]],
            [X.obsm["weighted_projections"][:, cmp2 - 1]],
        )
    )
    df = pd.DataFrame(
        data=cmpWeights.transpose(), columns=[f"Cmp. {cmp1}", f"Cmp. {cmp2}"]
    )
    df = df.sample(frac=frac)

    sns.scatterplot(data=df, x=f"Cmp. {cmp1}", y=f"Cmp. {cmp2}", ax=ax, color="k")
    ax.set(title=f"WP {frac*100}% of Cells")


def remove_low_pc_cmp(pc_df: pd.DataFrame, abs_threshold=0.4):
    """Removes partial correlation values below absolute threshold"""
    pc_df = pc_df.where(np.triu(np.ones(pc_df.shape)).astype(bool))
    pc_df = pc_df.where(np.identity(pc_df.shape[0]) != 1, np.nan)
    pc_df = pc_df.stack().reset_index()

    pc_df.columns = ["Var1", "Var2", "Weight"]

    pc_df_top = pc_df.loc[pc_df["Weight"] > abs_threshold]
    pc_df_bot = pc_df.loc[pc_df["Weight"] < -abs_threshold]
    pc_df = pd.concat([pc_df_bot, pc_df_top])

    return pc_df
