from concurrent.futures import ProcessPoolExecutor
from parafac2.normalize import prepare_dataset
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_array, spmatrix
from sklearn.preprocessing import scale
from sklearn.utils.sparsefuncs import inplace_column_scale, mean_variance_axis


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


def pseudobulk_lupus(X, cellType="Cell Type"):
    """Average gene expression for each condition and cell type;
    creates matrix and tensor version"""
    X_df = X.to_df()
    X_df = X_df.subtract(X.var["means"].values)
    X_df["Condition"] = X.obs["Condition"].values
    X_df["Cell Type"] = X.obs[cellType].values
    X_df["Status"] = X.obs["SLE_status"].values
    X_matrix = (
        X_df.groupby(["Condition", "Cell Type"], observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    conds = pd.unique(X_matrix["Condition"])
    celltypes = pd.unique(X_matrix["Cell Type"])
    genes = X.var_names.values

    status = []
    for i, cond in enumerate(conds):
        all_status = X_df.loc[X_df["Condition"] == cond]["Status"]
        status = np.append(status, np.unique(all_status))

    X_matrix["Status"] = np.repeat(status, len(celltypes))

    X_tensor = np.empty((len(conds), len(celltypes), len(genes)))
    X_tensor[:] = np.nan

    for i, cond in enumerate(conds):
        for j, celltype in enumerate(celltypes):
            specific_df = X_matrix.loc[
                (X_matrix["Condition"] == cond) & (X_matrix["Cell Type"] == celltype)
            ]
            X_tensor[i, j, :] = specific_df.iloc[0, 2:-1].to_numpy()

    return X_matrix, X_tensor
