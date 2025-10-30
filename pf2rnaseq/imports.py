from concurrent.futures import ProcessPoolExecutor

import anndata
import pandas as pd
import scanpy as sc
from parafac2.normalize import prepare_dataset
from pathlib import Path

path_here = Path(__file__).parent.parent


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
    X = X[:, ~X.var_names.str.match("^CMO3[0-9]{2}$")].copy()  # type: ignore

    return prepare_dataset(X, "Condition", geneThreshold=0.002)  # 0.1


def import_pf2Cytokine30() -> anndata.AnnData:
    """Import Meyer Cytokine PBMC dataset after pf2 run with 30 components.
    -- columns from observation data:
    {'Stimulation': Cytokine and Dose}
    """
    X = anndata.read_h5ad("/opt/extra-storage/pf2_results/cytok_pf2_30.h5ad")

    return X


def import_Heiser() -> anndata.AnnData:
    """Import Heiser C3TAg dataset.
    anndata.X is the raw counts

    """
    data = anndata.read_h5ad("/home/nicoleb/C3TAg.h5ad")

    return prepare_dataset(data, "sample_id", geneThreshold=0.1)


def import_MouseImmune(geneThreshold=0.1) -> anndata.AnnData:
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
    X = X[X.obs["celltype"] != "doublet", :].copy()

    return prepare_dataset(X, "biosample_id", geneThreshold=geneThreshold)  # 0.01


def import_Parse(geneThreshold=0.1, doublet=False) -> anndata.AnnData:
    """Import Parse data .
    cytokine: cytokine treatment
    donor: donor identifier

    """
    X = anndata.read_h5ad("/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Parse_Donor11.h5ad")
    if doublet:
        doubletDF = pd.read_csv(
        path_here / "pf2rnaseq/Data/DN11Doublets.csv.gz",
        index_col=0  
    )
        X.obs = X.obs.join(doubletDF.reindex(X.obs.index))
        singlet_mask = X.obs["doublet"] == 0
        X = X[singlet_mask, :].copy()
        print(f"Kept {X.n_obs} singlet cells, removed {(~singlet_mask).sum()} doublets")

    return prepare_dataset(X, "cytokine", geneThreshold=geneThreshold)
