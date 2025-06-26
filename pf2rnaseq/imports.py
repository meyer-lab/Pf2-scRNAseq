import glob
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import anndata
from pathlib import Path

import numpy as np
import pandas as pd
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_array, csr_matrix, spmatrix
from sklearn.preprocessing import scale
from sklearn.utils.sparsefuncs import inplace_column_scale, mean_variance_axis


def prepare_dataset_deviance(
    X: anndata.AnnData, condition_name, geneThreshold
) -> anndata.AnnData:
    X.X = csr_array(X.X)  # type: ignore
    assert np.amin(X.X.data) >= 0.0

    # Remove cells and genes with fewer than 10 reads
    # X = X[X.X.sum(axis=1) > 10, X.X.sum(axis=0) > 10]
    readmean, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    X = X[:, readmean > geneThreshold]

    # Copy so that the subsetting is preserved
    X._init_as_actual(X.copy())

    # deviance transform
    y_ij = X.X.toarray()  # type: ignore

    # counts per cell
    n_i = y_ij.sum(axis=1)

    # MLE of gene expression
    pi_j = y_ij.sum(axis=0) / np.sum(n_i)

    non_y_ij = n_i[:, None] - y_ij
    mu_ij = n_i[:, None] * pi_j[None, :]
    signs = np.sign(y_ij - mu_ij)

    first_term = 2 * y_ij * np.log(np.maximum(y_ij, 1.0) / mu_ij)
    second_term = 2 * non_y_ij * np.log(non_y_ij / (n_i[:, None] - mu_ij))

    X.X = signs * np.sqrt(np.maximum(first_term + second_term, 0.0))

    X.X = scale(X.X)

    _, sgIndex = np.unique(X.obs_vector(condition_name), return_inverse=True)
    X.obs["condition_unique_idxs"] = sgIndex
    X.obs["condition_unique_idxs"] = X.obs["condition_unique_idxs"].astype("category")

    # Pre-calculate gene means
    # X.var["means"] = np.zeros(X.shape[1])

    assert np.all(np.isfinite(X.X))  # type: ignore
    return X


def prepare_dataset(
    X: anndata.AnnData, condition_name: str, geneThreshold: float
) -> anndata.AnnData:
    assert isinstance(X.X, spmatrix)
    assert np.amin(X.X.data) >= 0.0  # type: ignore

    # Filter out genes with too few reads
    readmean, _ = mean_variance_axis(X.X, axis=0)  # type: ignore
    X = X[:, readmean > geneThreshold]
    X._init_as_actual(X.copy())
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
    X.obs["condition_unique_idxs"] = X.obs["condition_unique_idxs"].astype("category")

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
        return prepare_dataset_deviance(data, "sample_id", geneThreshold=0.1)
    else:
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
    ...}
    """
    X = anndata.read_h5ad("/home/nicoleb/MouseCytok.h5ad")

    return prepare_dataset(X, "biosample_id", geneThreshold=0.1)  # 0.01



