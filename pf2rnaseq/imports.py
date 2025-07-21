
from parafac2.normalize import prepare_dataset
import anndata




def import_cytokine(deviance=False) -> anndata.AnnData:
    """Import Meyer Cytokine PBMC dataset.
    -- columns from observation data:
    {'Stimulation': Cytokine and Dose}
    """
    X = anndata.read_h5ad("/opt/extra-storage/Treg_h5ads/Treg_raw.h5ad")

    # Remove multiplexing identifiers
    X = X[:, ~X.var_names.str.match("^CMO3[0-9]{2}$")]  # type: ignore

    return prepare_dataset(
        X, "Condition", geneThreshold=0.002, deviance=deviance
    )  # 0.1


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

    return prepare_dataset(
        data, "sample_id", geneThreshold=0.1, deviance=deviance
    )  


def import_MouseImmune(geneThreshold=0.1, deviance=False) -> anndata.AnnData:
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
    X = X.copy()

    return prepare_dataset(X, "biosample_id", geneThreshold, deviance=deviance)  # 0.01
