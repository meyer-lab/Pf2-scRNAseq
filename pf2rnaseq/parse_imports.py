"""Provides functions to import various information from the Parse Biosciences dataset."""

import anndata
import pandas as pd


def get_labels(path: str, obs_column: str, unique: bool = True) -> pd.Series:
    """Get labels from the AnnData object using backed mode to avoid loading full dataset.

    Parameters:
    -----------
    path : str
        Path to the .h5ad file
    obs_column : str
        Name of the observation column to extract
    unique : bool
        If True, return unique values from the specified obs column

    Returns:
    --------
    pd.Series
    Union[pd.Series, np.ndarray]
        If unique is False, returns a pandas Series containing the values from the specified obs column.
        If unique is True, returns a numpy array of unique values from the specified obs column.
    """
    # Open in backed mode - only loads metadata, not the expression matrix
    adata = anndata.read_h5ad(path, backed="r")

    # Extract the specific column from obs
    if obs_column not in adata.obs.columns:
        raise KeyError(
            f"Column '{obs_column}' not found in obs. Available columns: {list(adata.obs.columns)}"
        )

    labels = adata.obs[obs_column].copy()
    if unique:
        labels = labels.unique()

    return labels


def get_cells(
    path: str,
    donor: str,
    cytokine: str,
    donor_column: str = "donor",
    cytokine_column: str = "cytokine",
) -> anndata.AnnData:
    """Get cells matching specific donor and cytokine using backed mode.

    Parameters:
    -----------
    path : str
        Path to the .h5ad file
    donor : str
        Donor identifier to filter by
    cytokine : str
        Cytokine identifier to filter by
    donor_column : str
        Name of the observation column containing donor info (default: "donor")
    cytokine_column : str
        Name of the observation column containing cytokine info (default: "cytokine")

    Returns:
    --------
    anndata.AnnData
        Subset of AnnData containing only cells matching the donor and cytokine criteria
    """
    # Open in backed mode - only loads metadata initially
    adata = anndata.read_h5ad(path, backed="r")

    # Check if required columns exist
    if donor_column not in adata.obs.columns:
        raise KeyError(
            f"Column '{donor_column}' not found in obs. Available columns: {list(adata.obs.columns)}"
        )

    if cytokine_column not in adata.obs.columns:
        raise KeyError(
            f"Column '{cytokine_column}' not found in obs. Available columns: {list(adata.obs.columns)}"
        )

    # Create boolean mask for filtering
    donor_mask = adata.obs[donor_column] == donor
    cytokine_mask = adata.obs[cytokine_column] == cytokine
    combined_mask = donor_mask & cytokine_mask

    # Check if any cells match the criteria
    if not combined_mask.any():
        print(
            f"Warning: No cells found matching donor='{donor}' and cytokine='{cytokine}'"
        )

    # Subset the data - this will load only the required portion
    filtered_adata = adata[combined_mask].to_memory()
    return filtered_adata
