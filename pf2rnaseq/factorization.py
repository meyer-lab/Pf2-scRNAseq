import anndata
import cupy
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sps
from pacmap import PaCMAP
from parafac2.normalize import prepare_dataset
from parafac2.parafac2 import parafac2_nd, store_pf2
from scipy.stats import gmean
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from tensorly.cp_tensor import CPTensor
from tlviz.factor_tools import factor_match_score as fms
from tqdm import tqdm


def correct_conditions(X: anndata.AnnData):
    """Correct the conditions factors by overall read depth. Ensures that weighting is not affected by cell count difference"""
    # sgIndex = X.obs["condition_unique_idxs"]
    sgIndex = X.obs["condition_unique_idxs"].cat.codes
    counts = np.zeros((np.amax(sgIndex) + 1, 1))

    cond_mean = gmean(X.uns["Pf2_A"], axis=1)

    x_count = X.X.sum(axis=1)

    for ii in range(counts.size):
        counts[ii] = np.sum(x_count[X.obs["condition_unique_idxs"] == ii])

    lr = LinearRegression()
    lr.fit(counts, cond_mean.reshape(-1, 1))

    counts_correct = lr.predict(counts)

    return X.uns["Pf2_A"] / counts_correct


def pf2(
    X: anndata.AnnData,
    rank: int,
    random_state=1,
    doEmbedding: bool = True,
    tolerance=1e-9,
    r2x=False,
):
    cupy.cuda.Device(1).use()
    pf_out, R2X = parafac2_nd(
        X,
        rank=rank,
        random_state=random_state,
        tol=tolerance,
        n_iter_max=500,
    )

    X = store_pf2(X, pf_out)

    if doEmbedding:
        pcm = PaCMAP(random_state=random_state)
        X.obsm["X_pf2_PaCMAP"] = pcm.fit_transform(X.obsm["projections"])  # type: ignore
    if r2x:
        return X, R2X
    else:
        return X


def pf2_pca_r2x(X: anndata.AnnData, ranks):
    X = X.to_memory()
    XX = sps.csr_array(X.X)

    r2x_pf2 = np.zeros(len(ranks))

    for i in tqdm(range(len(r2x_pf2)), total=len(r2x_pf2)):
        _, R2X = parafac2_nd(X, rank=i + 1)
        r2x_pf2[i] = R2X

    pca = PCA(n_components=ranks[-1], svd_solver="arpack")
    pca.fit(XX)
    r2x_pca = np.cumsum(pca.explained_variance_ratio_)

    return r2x_pf2, r2x_pca[np.array(ranks) - 1]


def calculateFMS(A: anndata.AnnData, B: anndata.AnnData):
    """Calculates FMS between 2 factors"""
    factors = [A.uns["Pf2_A"], A.uns["Pf2_B"], A.varm["Pf2_C"]]
    A_CP = CPTensor(
        (
            A.uns["Pf2_weights"],
            factors,
        )
    )

    factors = [B.uns["Pf2_A"], B.uns["Pf2_B"], B.varm["Pf2_C"]]
    B_CP = CPTensor(
        (
            B.uns["Pf2_weights"],
            factors,
        )
    )

    return fms(A_CP, B_CP, consider_weights=False, skip_mode=1)  # type: ignore


def fms_percent_drop(
    X: anndata.AnnData,
    percentList: np.ndarray,
    runs: int,
    rank: int = 30,
):
    # Plots FMS score when percentage is removed from data
    dataX = pf2(X, rank, doEmbedding=False)

    fmsLists = []

    for j in range(0, runs, 1):
        scores = [1.0]

        for i in percentList[1:]:
            sampled_data: anndata.AnnData = sc.pp.subsample(
                X, fraction=1 - (i / 100), random_state=j, copy=True
            )  # type: ignore
            sampledX = pf2(sampled_data, rank, random_state=j + 2, doEmbedding=False)

            fmsScore = calculateFMS(dataX, sampledX)
            scores.append(fmsScore)

        fmsLists.append(scores)

    runsList_df = []
    for i in range(0, runs):
        for _j in range(0, len(percentList)):
            runsList_df.append(i)
    percentList_df = []
    for _i in range(0, runs):
        for j in range(0, len(percentList)):
            percentList_df.append(percentList[j])
    fmsList_df = []
    for sublist in fmsLists:
        fmsList_df += sublist
    df = pd.DataFrame(
        {
            "Run": runsList_df,
            "Percentage of Data Dropped": percentList_df,
            "FMS": fmsList_df,
        }
    )

    return df


def resample(data: anndata.AnnData) -> anndata.AnnData:
    """Bootstrapping dataset"""
    indices = np.random.randint(0, data.shape[0], size=(data.shape[0],))
    data = data[indices].copy()
    return data


def fms_diff_ranks(
    X: anndata.AnnData,
    ranksList: list[int],
    runs: int,
):
    # Plots FMS when using different Pf2 components
    fmsLists = []

    for j in range(0, runs, 1):
        scores = []
        for i in ranksList:
            dataX = pf2(X, rank=i, random_state=j, doEmbedding=False)

            sampledX = pf2(resample(X), rank=i, random_state=j, doEmbedding=False)

            fmsScore = calculateFMS(dataX, sampledX)
            scores.append(fmsScore)
        fmsLists.append(scores)

    runsList_df = []
    for i in range(0, runs):
        for _j in range(0, len(ranksList)):
            runsList_df.append(i)
    ranksList_df = []
    for _i in range(0, runs):
        for j in range(0, len(ranksList)):
            ranksList_df.append(ranksList[j])
    fmsList_df = []
    for sublist in fmsLists:
        fmsList_df += sublist
    df = pd.DataFrame(
        {"Run": runsList_df, "Component": ranksList_df, "FMS": fmsList_df}
    )

    return df


def downsample_counts_multinomial(
    X: anndata.AnnData,
    percent_drop: float,
    random_state: int = 0,
) -> anndata.AnnData:
    """
    Create a downsampled counts copy of AnnData using multinomial sampling.

    Parameters:
    -----------
    X : anndata.AnnData
        Input dataset
    percent_drop : float
        Percentage of counts to drop (0-100)
    random_state : int
        Random seed for reproducibility

    Returns:
    --------
    anndata.AnnData
        Downsampled copy of the input data
    """
    import scipy.sparse as sp

    # Handle 0% drop case
    if percent_drop == 0:
        return X.copy()

    # Set random seed
    np.random.seed(random_state)

    # Convert to CSR and extract structure
    original_csr = X.X.tocsr()
    data = original_csr.data.copy()
    indices = original_csr.indices
    indptr = original_csr.indptr

    # Process each cell
    for cell_idx in range(X.n_obs):
        start_idx = indptr[cell_idx]
        end_idx = indptr[cell_idx + 1]

        if start_idx == end_idx:
            continue

        cell_data = data[start_idx:end_idx]
        total_counts = int(np.sum(cell_data))

        if total_counts == 0:
            continue

        new_total = int(total_counts * (1 - percent_drop / 100))
        if new_total == 0:
            data[start_idx:end_idx] = 0
            continue

        # Convert to probabilities and normalize
        probs = cell_data / total_counts
        probs = probs / np.sum(probs)  # Ensure sum = 1.0

        # Multinomial sampling
        new_counts = np.random.multinomial(new_total, probs)
        data[start_idx:end_idx] = new_counts.astype(cell_data.dtype)

    # Create new sparse matrix
    sampled_csr = sp.csr_matrix((data, indices, indptr), shape=original_csr.shape)

    # Create new AnnData object
    sampled_data = X.copy()
    sampled_data.X = sampled_csr

    return sampled_data


def calculate_fms_downsample(
    X: anndata.AnnData,
    X_pf2: anndata.AnnData,
    percent_drop: float,
    rank: int = 30,
    deviance: bool = False,
    condition: str = "Condition",
    random_state: int = 0,
) -> float:
    """
    Calculate FMS for a single downsampling scenario.

    Parameters:
    -----------
    X : anndata.AnnData
        Original dataset for reference
    X_pf2 : anndata.AnnData
        Full factorized dataset
    percent_drop : float
        Percentage of counts to drop (0-100)
    rank : int
        Factorization rank
    deviance : bool
        Whether to use deviance normalization
    condition : str
        Condition column name
    random_state : int
        Random seed

    Returns:
    --------
    float
        FMS score
    """

    # Handle 0% drop case
    if percent_drop == 0:
        return 1.0

    # Create downsampled data
    sampled_data = downsample_counts_multinomial(
        X, percent_drop, random_state=random_state
    )

    # Apply same processing as reference
    sampled_data = prepare_dataset(
        sampled_data, condition, geneThreshold=0.0, deviance=deviance
    )

    # Factorization
    sampledX = pf2(sampled_data, rank, random_state=random_state + 2, doEmbedding=False)

    return calculateFMS(X_pf2, sampledX)


def fms_percent_drop_counts(
    X: anndata.AnnData,
    percentList: np.ndarray,
    rank: int = 30,
    deviance: bool = False,
    condition: str = "Condition",
    geneThreshold: float = 0.0,
    random_state: int = 0,
) -> pd.DataFrame:
    """
    Calculate FMS for multiple downsampling percentages (single run).

    Parameters:
    -----------
    X : anndata.AnnData
        Input dataset
    percentList : np.ndarray
        Array of dropout percentages to test
    rank : int
        Factorization rank
    deviance : bool
        Whether to use deviance normalization
    condition : str
        Condition column name
    geneThreshold : float
        Gene threshold for preparation
    random_state : int
        Random seed

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: Percentage of Counts Dropped, FMS
    """
    results = []
    X_prepared = prepare_dataset(
        X, condition, geneThreshold=geneThreshold, deviance=deviance
    )
    X_pf2 = pf2(X_prepared, rank, doEmbedding=False)

    for percent_drop in percentList:
        fms_score = calculate_fms_downsample(
            X=X,
            X_pf2=X_pf2,
            percent_drop=percent_drop,
            rank=rank,
            deviance=deviance,
            condition=condition,
            random_state=random_state,
        )

        results.append({"Percentage of Counts Dropped": percent_drop, "FMS": fms_score})

    return pd.DataFrame(results)
