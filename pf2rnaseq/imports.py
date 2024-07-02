import glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import anndata
import scanpy as sc
from scipy.sparse import spmatrix, csr_matrix
from sklearn.utils.sparsefuncs import inplace_column_scale, mean_variance_axis
import pandas as pd


def prepare_dataset(
    X: anndata.AnnData, condition_name: str, geneThreshold: float
) -> anndata.AnnData:
    assert isinstance(X.X, spmatrix)
    assert np.amin(X.X.data) >= 0.0  # type: ignore

    # Filter out genes with too few reads
    readmean, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    X = X[:, readmean > geneThreshold]

    # Normalize read depth
    sc.pp.normalize_total(X, exclude_highly_expressed=False, inplace=True)

    # Scale genes by sum
    readmean, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    readsum = X.shape[0] * readmean
    inplace_column_scale(X.X, 1.0 / readsum)

    # Transform values
    X.X.data = np.log10((1000.0 * X.X.data) + 1.0)  # type: ignore

    # Get the indices for subsetting the data
    _, sgIndex = np.unique(X.obs_vector(condition_name), return_inverse=True)
    X.obs["condition_unique_idxs"] = sgIndex

    # Pre-calculate gene means
    means, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    X.var["means"] = means

    return X


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


def import_HTAN() -> anndata.AnnData:
    """Imports Vanderbilt's HTAN 10X data."""
    files = glob.glob("/opt/extra-storage/HTAN/*.mtx.gz")
    futures = []
    data = {}

    with ProcessPoolExecutor(max_workers=10) as executor:
        for filename in files:
            future = executor.submit(
                sc.read_10x_mtx,
                "/opt/extra-storage/HTAN/",
                gex_only=False,
                make_unique=True,
                prefix=filename.split("/")[-1].split("matrix.")[0],
            )
            futures.append(future)

        for i, k in enumerate(files):
            result = futures[i].result()
            data[k.split("/")[-1].split("_matrix.")[0]] = result

    X = anndata.concat(data, merge="same", label="Condition")

    return prepare_dataset(X, "Condition", geneThreshold=0.1)


def import_CCLE() -> anndata.AnnData:
    """Imports barcoded cell data."""
    # TODO: Still need to add gene names and barcodes.
    folder = "/opt/extra-storage/asm/Heiser-barcode/CCLE/"

    adatas = {
        "HCT116_1": anndata.read_text(
            Path(folder + "HCT116_tracing_T1.count_mtx.tsv")
        ).T,
        "HCT116_2": anndata.read_text(
            Path(folder + "HCT116_tracing_T2.count_mtx.tsv")
        ).T,
        "MDA-MB-231_1": anndata.read_text(
            Path(folder + "MDA-MB-231_tracing_T1.count_mtx.tsv")
        ).T,
        "MDA-MB-231_2": anndata.read_text(
            Path(folder + "MDA-MB-231_tracing_T2.count_mtx.tsv")
        ).T,
    }

    X = anndata.concat(adatas, label="sample")
    X.X = csr_matrix(X.X)

    return prepare_dataset(X, "sample", geneThreshold=0.1)


def import_cytokine() -> anndata.AnnData:
    """Import Meyer Cytokine PBMC dataset.
    -- columns from observation data:
    {'Stimulation': Cytokine and Dose}
    """
    X = anndata.read_h5ad("/opt/extra-storage/Treg_h5ads/Treg_raw.h5ad")

    # Remove multiplexing identifiers
    X = X[:, ~X.var_names.str.match("^CMO3[0-9]{2}$")]  # type: ignore

    return prepare_dataset(X, "Condition", geneThreshold=0.05)


def pseudobulk_lupus(X, cellType="Cell Type"):
    """Average gene expression for each condition and cell type; 
    creates matrix and tensor version"""
    X_df = X.to_df()
    X_df = X_df.subtract(X.var["means"].values)
    X_df["Condition"] = X.obs["Condition"].values
    X_df["Cell Type"] = X.obs[cellType].values
    X_df["Status"] = X.obs["SLE_status"].values
    X_matrix = X_df.groupby(["Condition", "Cell Type"], observed=False).mean(numeric_only=True).reset_index()

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
            specific_df = X_matrix.loc[(X_matrix["Condition"] == cond) & (X_matrix["Cell Type"] == celltype)] 
            X_tensor[i, j, :] = specific_df.iloc[0, 2:-1].to_numpy()
          
    return X_matrix, X_tensor

