import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Patch

from ..common import highly_weighted_cytokines

cmap = sns.diverging_palette(240, 10, as_cmap=True)


def plot_condition_factors(
    data: AnnData,
    ax: Axes,
    cond_group_labels: pd.Series | None = None,
    groupConditions=False,
    cond="Condition",
    log_scale=True
):
    """Plots Pf2 condition factors"""
    pd.set_option("display.max_rows", None)
    yt = pd.Series(np.unique(data.obs[cond]))
    X = np.array(data.uns["Pf2_A"])
    if log_scale:
        X = np.log10(X)
    
    X -= np.median(X, axis=0)
    X /= np.std(X, axis=0)+1e-3
    ind = reorder_table(X+1e-3)

    X = X[ind]
    yt = yt.iloc[ind]

    if cond_group_labels is not None:
        cond_group_labels = cond_group_labels.iloc[ind]
        if groupConditions is True:
            ind = cond_group_labels.argsort()
            cond_group_labels = cond_group_labels.iloc[ind]
            X = X[ind]
            yt = yt.iloc[ind]
        ax.tick_params(axis="y", which="major", pad=20, length=0)

        # extra padding to leave room for the row colors
        # get list of colors for each label:
        colors = sns.color_palette(
            n_colors=pd.Series(cond_group_labels).nunique()
        ).as_hex()
        lut = {}
        legend_elements = []
        for index, group in enumerate(pd.Series(cond_group_labels).unique()):
            lut[group] = colors[index]
            legend_elements.append(Patch(color=colors[index], label=group))
        row_colors = pd.Series(cond_group_labels).map(lut)
        for iii, color in enumerate(row_colors):
            ax.add_patch(
                plt.Rectangle(
                    xy=(-0.05, iii),
                    width=0.05,
                    height=1,
                    color=color,
                    lw=0,
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )
            )
        # add a little legend
        # ax.legend(handles=legend_elements, bbox_to_anchor=(0, 1.3))

    xticks = np.arange(1, X.shape[1] + 1)

    sns.heatmap(
        data=X,
        xticklabels=xticks,
        yticklabels=yt,
        ax=ax,
        center=0,
        cmap=cmap,
    )
    ax.tick_params(axis="y", rotation=0)
    ax.set(xlabel="Component")


def plot_condition_factors_groups(
    data: AnnData,
    ax: Axes,
    cond_group_labels: pd.Series | None = None,
    subgroup_labels: pd.Series | None = None,
    groupConditions=False,
    cond="Condition",
    main_group_title="Treatment",
    subgroup_title="Tumor Type",
    log_scale=True,
    sub_leg=True,
):
    """
    Plots Pf2 condition factors with two-level grouping capability.

    Parameters:
    -----------
    data: AnnData object containing the Pf2 results
    ax: Matplotlib axes to plot on
    cond_group_labels: Primary grouping labels (ie. treatments)
    subgroup_labels: Secondary grouping labels (ie. tumor types)
    groupConditions: Whether to sort conditions by groups
    cond: Column name in obs containing condition information
    main_group_title: Title for the main group legend
    subgroup_title: Title for the subgroup legend
    log_scale: Whether to apply log10 transformation to the data
    sub_leg: Whether to show the subgroup legend
    """
    pd.set_option("display.max_rows", None)
    yt = pd.Series(np.unique(data.obs[cond]))
    X = np.array(data.uns["Pf2_A"])
    if log_scale:
        # Apply log10 transformation for better visualization
        assert np.all(X >= 0)
        X = np.log10(X)

    X -= np.median(X, axis=0)
    X /= np.std(X, axis=0)

    if log_scale is False:
        X -= np.min(X, axis=0)

    # Hierarchically cluster conditions
    ind = reorder_table(X)
    X = X[ind]
    yt = yt.iloc[ind]

    # Sort by group labels if provided
    if cond_group_labels is not None:
        cond_group_labels = cond_group_labels.iloc[ind]

        # Initialize subgroup labels if provided
        if subgroup_labels is not None:
            subgroup_labels = subgroup_labels.iloc[ind]

        if groupConditions is True:
            # First create a composite index for sorting
            # We'll convert each group to a categorical and get codes
            if subgroup_labels is not None:
                # Create a DataFrame for sorting
                sort_df = pd.DataFrame(
                    {
                        "main_group": cond_group_labels,
                        "subgroup": subgroup_labels,
                        "original_idx": np.arange(len(cond_group_labels)),
                    }
                )

                # Sort first by main group, then by subgroup
                sort_df = sort_df.sort_values(["main_group", "subgroup"])

                # Get the reordering index
                composite_idx = sort_df["original_idx"].values

                # Apply the sort
                cond_group_labels = cond_group_labels.iloc[composite_idx]
                subgroup_labels = subgroup_labels.iloc[composite_idx]
                X = X[composite_idx]
                yt = yt.iloc[composite_idx]
            else:
                # If no subgroups, just sort by main groups
                ind = cond_group_labels.argsort()
                cond_group_labels = cond_group_labels.iloc[ind]
                X = X[ind]
                yt = yt.iloc[ind]

        # Add padding for row colors
        ax.tick_params(
            axis="y",
            which="major",
            pad=30 if subgroup_labels is not None else 30,
            length=0,
        )

        # Main group colors - use a distinct palette
        main_colors = sns.color_palette(
            "tab10", n_colors=pd.Series(cond_group_labels).nunique()
        ).as_hex()
        main_lut = {}
        main_legend_elements = []
        for index, group in enumerate(pd.Series(cond_group_labels).unique()):
            main_lut[group] = main_colors[index]
            main_legend_elements.append(Patch(color=main_colors[index], label=group))

        main_row_colors = pd.Series(cond_group_labels).map(main_lut)

        # Add color patches for main groups
        for i, color in enumerate(main_row_colors):
            ax.add_patch(
                plt.Rectangle(
                    xy=(-0.02, i),
                    width=0.02,
                    height=1,
                    color=color,
                    lw=0,
                    transform=ax.get_yaxis_transform(),
                    clip_on=False,
                )
            )

        # Add subgroup colors if provided
        if subgroup_labels is not None:
            # Use a different color palette for subgroups
            sub_colors = sns.color_palette(
                "Set3", n_colors=pd.Series(subgroup_labels).nunique()
            ).as_hex()
            sub_lut = {}
            sub_legend_elements = []
            for index, group in enumerate(pd.Series(subgroup_labels).unique()):
                sub_lut[group] = sub_colors[index]
                sub_legend_elements.append(Patch(color=sub_colors[index], label=group))

            sub_row_colors = pd.Series(subgroup_labels).map(sub_lut)

            # Add color patches for subgroups
            for i, color in enumerate(sub_row_colors):
                ax.add_patch(
                    plt.Rectangle(
                        xy=(-0.04, i),  # Position to left of main group colors
                        width=0.02,
                        height=1,
                        color=color,
                        lw=0,
                        transform=ax.get_yaxis_transform(),
                        clip_on=False,
                    )
                )

            # Add legends for both groupings with proper titles
            main_legend = ax.legend(
                handles=main_legend_elements,
                bbox_to_anchor=(0, 1.3),
                title=main_group_title,
                loc="upper left",
            )
            ax.add_artist(main_legend)  # Add first legend
            if sub_leg:
                # Add second legend for subgroups
                ax.legend(
                    handles=sub_legend_elements,
                    bbox_to_anchor=(0.5, 1.3),
                    title=subgroup_title,
                    loc="upper left",
                )

        else:
            # Add only main group legend if no subgroups
            ax.legend(
                handles=main_legend_elements,
                bbox_to_anchor=(0, 1.3),
                title=main_group_title,
            )

    # Create the heatmap
    xticks = np.arange(1, X.shape[1] + 1)
    sns.heatmap(
        data=X,
        xticklabels=xticks,
        yticklabels=yt,
        ax=ax,
        center=0,
        cmap=cmap,
    )
    ax.tick_params(axis="y", rotation=0)
    ax.set(xlabel="Component")


def plot_eigenstate_factors(data: AnnData, ax: Axes):
    """Plots Pf2 eigenstate factors"""
    rank = data.uns["Pf2_B"].shape[1]
    xticks = np.arange(1, rank + 1)
    X = data.uns["Pf2_B"]
    X = X / np.max(np.abs(np.array(X)))
    yt = np.arange(1, rank + 1)

    sns.heatmap(
        data=X,
        xticklabels=xticks,
        yticklabels=yt,
        ax=ax,
        center=0,
        cmap=cmap,
        vmin=-1,
        vmax=1,
    )
    ax.set(xlabel="Component")


def plot_gene_factors(
    data: AnnData, ax: Axes, trim=True
):  # yt gene names- input that will
    """Plots Pf2 gene factors"""
    rank = data.varm["Pf2_C"].shape[1]
    X = np.array(data.varm["Pf2_C"])
    yt = data.var.index.values
    if trim is True:
        max_weight = np.max(np.abs(X), axis=1)
        kept_idxs = max_weight > 0.08  # adjust this to sdjust amount of genes included
        X = X[kept_idxs]
        yt = yt[kept_idxs]
    # index for genes
    ind = reorder_table(X)
    X = X[ind]
    X = X / np.max(np.abs(X))
    yt = [yt[ii] for ii in ind]
    xticks = np.arange(1, rank + 1)

    sns.heatmap(
        data=X,
        xticklabels=xticks,
        yticklabels=yt,
        ax=ax,
        center=0,
        cmap=cmap,
        vmin=-1,
        vmax=1,
    )

    ax.set(xlabel="Component")


def plot_geneSet_factors(
    data: AnnData, ax: Axes, genes: np.array, trim=True
):  # yt gene names- input that will
    """Plots Pf2 gene factors for a set of genes"""
    rank = data.varm["Pf2_C"].shape[1]
    X = np.array(data.varm["Pf2_C"])
    yt = data.var.index.values

    kept_idxs = np.where(np.in1d(yt, genes))
    X = X[kept_idxs]
    yt = yt[kept_idxs]

    X = X / np.max(np.abs(X))

    xticks = np.arange(1, rank + 1)

    sns.heatmap(
        data=X,
        xticklabels=xticks,
        yticklabels=yt,
        ax=ax,
        center=0,
        cmap=cmap,
        vmin=-1,
        vmax=1,
    )
    ax.set_xlabel("Component", fontsize=12)
    ax.set_ylabel("Gene", fontsize=12)


def plot_gene_factors_partial(
    cmp: int, dataIn: AnnData, ax: Axes, geneAmount: int = 5, top=True
):
    """Plotting weights for gene factors for both most negatively/positively weighted terms"""
    cmpName = f"Cmp. {cmp}"

    df = pd.DataFrame(
        data=dataIn.varm["Pf2_C"][:, cmp - 1], index=dataIn.var_names, columns=[cmpName]
    )
    df = df.reset_index(names="Gene")
    df = df.sort_values(by=cmpName)

    if top:
        sns.barplot(
            data=df.iloc[-geneAmount:, :], x="Gene", y=cmpName, color="k", ax=ax
        )
    else:
        sns.barplot(data=df.iloc[:geneAmount, :], x="Gene", y=cmpName, color="k", ax=ax)

    ax.tick_params(axis="x", rotation=90)


def plot_factor_weight(X: AnnData, ax: Axes):
    """Plots weights from Pf2 model"""
    df = pd.DataFrame(data=np.transpose(X.uns["Pf2_weights"]), columns=["Value"])
    df["Value"] = df["Value"] / np.max(df["Value"])
    df["Component"] = np.arange(1, len(X.uns["Pf2_weights"]) + 1)
    sns.barplot(data=df, x="Component", y="Value", ax=ax)
    ax.tick_params(axis="x", rotation=90)


def reorder_table(projs: np.ndarray) -> np.ndarray:
    """Reorder a table's rows using heirarchical clustering"""
    assert projs.ndim == 2
    Z = sch.linkage(projs, method="complete", metric="cosine", optimal_ordering=True)
    return sch.leaves_list(Z)


def bot_top_genes(X, cmp, geneAmount=5):
    """Saves most pos/negatively genes"""
    df = pd.DataFrame(
        data=X.varm["Pf2_C"][:, cmp - 1], index=X.var_names, columns=["Component"]
    )
    df = df.reset_index(names="Gene")
    df = df.sort_values(by="Component")

    top = df.iloc[-geneAmount:, 0].values
    bot = df.iloc[:geneAmount, 0].values
    all_genes = np.concatenate([bot, top])

    return all_genes


def plot_geneSetScore(
    data: AnnData, ax: Axes, genes: np.array, trim=True
):  # yt gene names- input that will
    """Plots Pf2 gene set score: sum of gene set weights for each component"""
    rank = data.varm["Pf2_C"].shape[1]
    X = np.array(data.varm["Pf2_C"])
    yt = data.var.index.values

    # Filter the genes
    kept_idxs = np.where(np.in1d(yt, genes))[0]
    X = X[kept_idxs]
    yt = yt[kept_idxs]

    # Calculate the sum of X values for each component
    component_sums = np.sum(X, axis=0)

    # Find the top 3 components with highest absolute scores
    top_3_indices = np.argsort(np.abs(component_sums))[-3:]

    # Create colors array - highlight the top 3 components
    colors = []
    for i in range(len(component_sums)):
        if i in top_3_indices:
            colors.append("darkred")
        else:
            colors.append("steelblue")

    # Create the bar plot
    xticks = np.arange(1, rank + 1)
    sns.barplot(x=xticks, y=component_sums, ax=ax)

    ax.set_xlabel("Component", fontsize=12)
    ax.set_ylabel("Sum of Weights", fontsize=12)
    ax.set_title("Sum of Gene Factors per Component", fontsize=12)
    # Add labels to the top 3 components
    for idx in top_3_indices:
        component_num = idx + 1
        y_pos = component_sums[idx] + 0.01 * np.sign(component_sums[idx]) * np.max(
            np.abs(component_sums)
        )
        ax.text(
            component_num,
            y_pos,
            f"Comp. {component_num}",
            ha="center",
            va="bottom" if component_sums[idx] > 0 else "top",
            fontsize=12,
            fontweight="bold",
            color="darkred",
        )

    ax.set_xlabel("Component", fontsize=20)
    ax.set_ylabel("Sum of Weights", fontsize=20)
    ax.set_title("Signature Score", fontsize=25)
    ax.tick_params(axis="x", rotation=90, labelsize=16)
    ax.tick_params(axis="y", labelsize=16)


def plot_ttest(X: AnnData, ax: Axes):
    """Create heatmap figure showing dominant cytokines across components."""

    # Get total number of components
    total_components = X.uns["Pf2_A"].shape[1]

    # Get all cytokines directly without separate function
    all_cytokines = X.obs["cyt"].unique()

    # Get highly weighted cytokines per component
    results_df = highly_weighted_cytokines(X)

    # Create pivot table for all components
    if len(results_df) > 0:
        pivot_pvalues = results_df.pivot_table(
            index="Cytokine",
            columns="Component",
            values="PostHoc_pvalue_corrected",
            fill_value=1.0,
        )
    else:
        # Create empty pivot table with all cytokines and components
        pivot_pvalues = pd.DataFrame(
            1.0, index=all_cytokines, columns=range(1, total_components + 1)
        )

    # Ensure all components are represented
    all_components = list(range(1, total_components + 1))
    pivot_pvalues = pivot_pvalues.reindex(columns=all_components, fill_value=1.0)

    # Convert to -log10 for visualization
    log_p_matrix = -np.log10(pivot_pvalues + 1e-10)

    # Create heatmap
    sns.heatmap(
        log_p_matrix,
        cmap="YlOrRd",
        ax=ax,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "-log10(p-value)"},
        linewidths=0.5,
    )

    ax.set_title("Highly weighted Cytokines Across Components (ANOVA + Post-hoc)")
    ax.set_xlabel("Component")
    ax.set_ylabel("Cytokine")


def plot_comp_weights(
    data: AnnData,
    ax: Axes,
    comp: int,
    cond="Condition",
    sort_bars=True,
    top_n=3,
    include_lowest=True,
):
    """Plots component weights for each condition as a bar chart"""

    # Get condition names and factor matrix
    cond_df = (
        data.obs[[cond, "condition_unique_idxs"]]
        .drop_duplicates()
        .sort_values("condition_unique_idxs")
    )
    yt = cond_df[cond].to_numpy()
    X = np.array(data.uns["Pf2_A"])
    cond_mapping = data.obs.groupby("condition_unique_idxs", sort=True)[cond].first()

    # Extract condition names and indices
    condition_indices = cond_mapping.index.to_numpy()
    yt = cond_mapping.values

    # Extract weights
    component_weights = X[condition_indices, comp - 1]
    # Create DataFrame for plotting
    df = pd.DataFrame({"Condition": yt, "Weight": component_weights})

    # Get top N highest weighted conditions
    top_n_highest = df.nlargest(top_n, "Weight")

    # Conditionally get lowest weighted conditions
    if include_lowest:
        top_n_lowest = df.nsmallest(top_n, "Weight")
        # Combine and keep only top conditions
        df_filtered = pd.concat([top_n_highest, top_n_lowest]).drop_duplicates()
    else:
        df_filtered = top_n_highest

    # Sort by weight if requested
    if sort_bars:
        df_filtered = df_filtered.sort_values("Weight", ascending=False)

    # Create color mapping - highest in red, lowest in blue (if included)
    colors = []
    for condition in df_filtered["Condition"]:
        if condition in top_n_highest["Condition"].values:
            colors.append("darkred")  # Highest weights
        else:
            colors.append("darkblue")  # Lowest weights

    # Create bar plot with custom colors
    bars = ax.bar(df_filtered["Condition"], df_filtered["Weight"], color=colors)

    # Customize the plot title based on whether lowest are included
    if include_lowest:
        title = f"Component {comp} Weights by Condition (Top {top_n} Highest/Lowest)"
    else:
        title = f"Component {comp} Weights by Condition (Top {top_n} Highest)"

    ax.set_title(title, fontsize=20)
    ax.set_xlabel("Condition", fontsize=18)
    ax.set_ylabel("Weight", fontsize=18)
    ax.tick_params(axis="x", rotation=90, labelsize=18)
    ax.tick_params(axis="y", labelsize=18)

    # Add legend for color coding (only if lowest are included)
    if include_lowest:
        legend_elements = [
            Patch(facecolor="darkred", label=f"Top {top_n} Highest"),
            Patch(facecolor="darkblue", label=f"Top {top_n} Lowest"),
        ]
        ax.legend(handles=legend_elements, loc="upper right")

    # Add horizontal line at y=0 for reference
    ax.axhline(y=0, color="gray", linestyle="-", alpha=0.3)
