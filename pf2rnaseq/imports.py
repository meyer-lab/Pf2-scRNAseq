from concurrent.futures import ProcessPoolExecutor

import anndata
import pandas as pd
import scanpy as sc
from parafac2.normalize import prepare_dataset


def import_citeseq() -> anndata.AnnData:
    """Imports 5 datasets from Hamad CITEseq."""
    files = ["control", "ic_pod1", "ic_pod7", "sc_pod1", "sc_pod7"]

    with ProcessPoolExecutor(max_workers=5) as executor:
        futures = [
            executor.submit(
                sc.read_10x_mtx,
                "/opt/andrew/HamadCITEseq/" + k,
                gex_only=False,
                make_unique=True,
            )
            for k in files
        ]

        data = {k: futures[i].result() for i, k in enumerate(files)}

    X = anndata.concat(data, merge="same", label="Condition")

    return prepare_dataset(X, "Condition", geneThreshold=0.1)


def import_cytokine() -> anndata.AnnData:
    """Import Meyer Cytokine PBMC dataset.
    -- columns from observation data:
    {'Stimulation': Cytokine and Dose}
    """
    X = anndata.read_h5ad("/opt/extra-storage/Treg_h5ads/Treg_raw.h5ad")

    # Remove multiplexing identifiers
    X = X[:, ~X.var_names.str.match("^CMO3[0-9]{2}$")]  # type: ignore

    return prepare_dataset(X, "Condition", geneThreshold=0.002)  # 0.1


def import_pf2Cytokine30() -> anndata.AnnData:
    """Import Meyer Cytokine PBMC dataset after pf2 run with 30 components.
    -- columns from observation data:
    {'Stimulation': Cytokine and Dose}
    """
    X = anndata.read_h5ad("/opt/extra-storage/pf2_results/cytok_pf2_30.h5ad")

    return X


def import_Heiser(deviance=False) -> anndata.AnnData:
    """Import Heiser C3TAg dataset.
    anndata.X is the raw counts

    """
    data = anndata.read_h5ad("/home/nicoleb/C3TAg.h5ad")
    if deviance:
        # Apply deviance transformation
        return prepare_dataset(data, "sample_id", geneThreshold=0.1, deviance=True)
    else:
        # Apply standard normalization and scaling
        return prepare_dataset(data, "sample_id", geneThreshold=0.1)


def import_MouseImmune() -> anndata.AnnData:
    """Import Mouse Immune Dictionary cytokine data.
     -- columns from observation data:
    {'biosample_id': cytokine and replicate info,
    'rep': replicate,
    'species': mouse species,
    'cytokine_family': cytokine family label,
    'cyt': cytokine mouse was treated with,
    'sex': sex of mouse,
    'celltype': cell type label,
    'organ__ontology_label': organ label,
    ...}"""
    X = anndata.read_h5ad("/home/nicoleb/MouseCytok.h5ad")
    # Filter out doublets
    X = X[X.obs["celltype"] != "doublet", :]

    return prepare_dataset(X, "biosample_id", geneThreshold=0.1)  # 0.01


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
        Series containing the values from the specified obs column
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
