"""
This file contains functions that are used in multiple figures.
"""

import sys
import time
from string import ascii_letters

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from scipy.stats import f_oneway, ttest_ind
from statsmodels.stats.multitest import multipletests

matplotlib.use("AGG")

matplotlib.rcParams["legend.labelspacing"] = 0.2
matplotlib.rcParams["legend.fontsize"] = 8
matplotlib.rcParams["xtick.major.pad"] = 1.0
matplotlib.rcParams["ytick.major.pad"] = 1.0
matplotlib.rcParams["xtick.minor.pad"] = 0.9
matplotlib.rcParams["ytick.minor.pad"] = 0.9
matplotlib.rcParams["legend.handletextpad"] = 0.5
matplotlib.rcParams["legend.handlelength"] = 0.5
matplotlib.rcParams["legend.framealpha"] = 0.5
matplotlib.rcParams["legend.markerscale"] = 0.7
matplotlib.rcParams["legend.borderpad"] = 0.35
matplotlib.rcParams["svg.fonttype"] = "none"


def getSetup(
    figsize: tuple[float, float], gridd: tuple[int, int]
) -> tuple[list[plt.Axes], Figure]:
    """Establish figure set-up with subplots."""
    sns.set_theme(
        style="whitegrid",
        font_scale=0.7,
        color_codes=True,
        palette="colorblind",
        rc={"grid.linestyle": "dotted", "axes.linewidth": 0.6},
    )

    # Setup plotting space and grid
    f = plt.figure(figsize=figsize, layout="constrained")
    gs1 = gridspec.GridSpec(gridd[0], gridd[1], figure=f)

    # Get list of axis objects
    ax = [f.add_subplot(gs1[x]) for x in range(gridd[0] * gridd[1])]

    return ax, f


def subplotLabel(axs: list[plt.Axes]):
    """Place subplot labels on figure."""
    for ii, ax in enumerate(axs):
        ax.text(
            -0.2,
            1.2,
            ascii_letters[0],  # put ii back in later
            transform=ax.transAxes,
            fontweight="bold",
            va="top",
        )


def genFigure():
    """Main figure generation function."""
    start = time.time()
    nameOut = "figure" + sys.argv[1]

    exec(f"from pf2rnaseq.figures.{nameOut} import makeFigure", globals())
    ff = makeFigure()  # type: ignore # noqa: F821

    if ff is not None:
        ff.savefig(
            f"./output/{nameOut}.svg", dpi=300, bbox_inches="tight", pad_inches=0
        )

    print(f"Figure {sys.argv[1]} is done after {time.time() - start} seconds.\n")


def get_condition_data(X):
    """Extract condition data with cytokine labels."""
    condition_data = []
    unique_condition_indices = X.obs["condition_unique_idxs"].unique()
    for cond_idx in unique_condition_indices:
        cells = X.obs[X.obs["condition_unique_idxs"] == cond_idx]
        if not cells.empty:
            condition_data.append(
                {"condition_idx": cond_idx, "cytokine": cells["cyt"].iloc[0]}
            )
    return pd.DataFrame(condition_data)


def find_dominant_cytokine_per_component(X, alpha=0.05):
    """Determine dominant cytokines using ANOVA followed by post-hoc pairwise comparisons."""
    condition_df = get_condition_data(X)
    pf2_a = X.uns["Pf2_A"]
    all_cytokines = condition_df["cytokine"].unique()

    results = []

    for comp_idx in range(pf2_a.shape[1]):
        cytokine_values = {}
        cytokine_groups = []

        # Collect values for each cytokine group in this component
        for cytokine in all_cytokines:
            indices = condition_df[condition_df["cytokine"] == cytokine][
                "condition_idx"
            ].values
            if len(indices) > 0:
                values = pf2_a[indices, comp_idx]
                if len(values) > 0:
                    cytokine_values[cytokine] = values
                    cytokine_groups.append(values)

        # First perform ANOVA to test for overall significance
        if len(cytokine_groups) >= 2 and all(
            len(group) > 0 for group in cytokine_groups
        ):
            f_stat, anova_p_value = f_oneway(*cytokine_groups)

            # If ANOVA is significant, perform post-hoc tests to find dominant cytokines
            if anova_p_value < alpha:
                for target_cytokine in all_cytokines:
                    if target_cytokine in cytokine_values:
                        target_values = cytokine_values[target_cytokine]
                        other_values = [
                            val
                            for cyt, vals in cytokine_values.items()
                            if cyt != target_cytokine
                            for val in vals
                        ]

                        if len(target_values) > 0 and len(other_values) > 0:
                            # Test if this cytokine is significantly higher than others
                            _, posthoc_p_value = ttest_ind(
                                target_values, other_values, alternative="greater"
                            )

                            results.append(
                                {
                                    "Component": comp_idx + 1,
                                    "Cytokine": target_cytokine,
                                    "ANOVA_pvalue": anova_p_value,
                                    "PostHoc_pvalue": posthoc_p_value,
                                }
                            )

    results_df = pd.DataFrame(results)

    # Apply multiple testing correction to post-hoc p-values
    if len(results_df) > 0:
        _, corrected_pvals, _, _ = multipletests(
            results_df["PostHoc_pvalue"], alpha=alpha, method="fdr_bh"
        )
        results_df["PostHoc_pvalue_corrected"] = corrected_pvals

    return results_df
