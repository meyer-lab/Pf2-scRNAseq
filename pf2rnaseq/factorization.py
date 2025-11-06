import anndata
import cupy
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sps
from pacmap import PaCMAP
from parafac2.parafac2 import parafac2_nd, store_pf2
from scipy.optimize import minimize
from scipy.stats import gmean
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from tensorly.cp_tensor import CPTensor
from tlviz.factor_tools import factor_match_score as fms
from tqdm import tqdm


def correct_conditions(X: anndata.AnnData):
    """Correct the conditions factors by overall read depth. Ensures that weighting is not affected by cell count difference"""
    sgIndex = X.obs["condition_unique_idxs"]
    # sgIndex = X.obs["condition_unique_idxs"].cat.codes
    counts = np.zeros((np.amax(sgIndex) + 1, 1))
    min_val = np.min(X.uns["Pf2_A"])
    if min_val < 0:
        # Add the absolute value of the minimum (plus a small epsilon) to make all values positive
        X.uns["Pf2_A"] = X.uns["Pf2_A"] + abs(min_val) + 1e-10
        print(
            f"Warning: Found negative values in Pf2_A (min: {min_val:.6f}). Added {abs(min_val) + 1e-10:.6f} to all values."
        )

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
    cupy.cuda.Device(0).use()
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
        _, R2X = parafac2_nd(X, rank=ranks[i])
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


def deconvolution_cytokine(
    A: np.ndarray,
    alpha: float = 0.1,
    max_iter: int = 5000,
    random_state: int = 1,
    beta: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Decompose cytokine factor matrix:  A ≈ W @ H

    This decomposes observed cytokine effects into:
    1. Direct primary effects (H)
    2. Induced effects via other cytokines (W)

    Parameters
    ----------
    A : np.ndarray
        Input matrix (n_cytokines, n_components)
        Example: (91 cytokines, 100 Parafac2 components)
    alpha : float
        Regularization strength
    max_iter : int
        Maximum optimization iterations
    random_state : int
        Random seed

    Returns
    -------
    W : np.ndarray
        Cytokine interaction matrix (n_cytokines, n_cytokines)
        W[i, j] = total contribution of cytokine j to observed effect of i
        Diagonal W[i,i] = direct effect of cytokine i
    H : np.ndarray
        Effect basis matrix (n_cytokines, n_components)
        H[:, j] = cytokine effects for component j without indirect contributions
    """
    n_cytokines, n_components = A.shape

    np.random.seed(random_state)

    # W initialized as identity, H is original A
    W_init = np.eye(n_cytokines)
    H_init = A.copy()

    x0 = np.concatenate([W_init.ravel(), H_init.ravel()])

    print("Cytokine deconvolution:")
    print(f"  A shape: {A.shape} (cytokines × components)")
    print(f"  W shape: ({n_cytokines}, {n_cytokines}) (cytokine interactions)")
    print(f"  H shape: ({n_cytokines}, {n_components}) (effect basis)")

    w_size = n_cytokines * n_cytokines
    iteration_counter = [0]
    best_loss = [np.inf]

    def objective(x):
        W = x[:w_size].reshape(n_cytokines, n_cytokines)
        H = x[w_size:].reshape(n_cytokines, n_components)

        # Reconstruction:A ≈ W @ H

        reconstruction = W @ H
        mse = np.sum((A - reconstruction) ** 2)

        # Regularization: L1 penalty on both W and H
        l1_W = alpha * np.sum(np.abs(W))
        l1_H = alpha * np.sum(np.abs(H))

        total_loss = mse + l1_W + l1_H

        iteration_counter[0] += 1
        if total_loss < best_loss[0]:
            best_loss[0] = total_loss

        if iteration_counter[0] % 100 == 0:
            print(
                f"  Iter {iteration_counter[0]}: Loss={total_loss:.4f} "
                f"(MSE={mse:.4f}, L1_W={l1_W:.4f}, L1_H={l1_H:.4f})"
            )

        return total_loss

    def gradient(x):
        W = x[:w_size].reshape(n_cytokines, n_cytokines)
        H = x[w_size:].reshape(n_cytokines, n_components)

        # ===== Gradient w.r.t. W =====
        # 1. Reconstruction term: ∂/∂W [||A - WH||²] = 2(error @ H^T), L1 penalty: ∂/∂W [α||W||₁] = α * sign(W)
        grad_W = 2 * ((W @ H - A) @ H.T) + alpha

        # ===== Gradient w.r.t. H =====
        # 1. Reconstruction term: ∂/∂H [||A - WH||²] = 2(W^T @ error),  L1 penalty: ∂/∂H [α||H||₁] = α * sign(H)
        grad_H = 2 * (W.T @ (W @ H - A)) + alpha

        return np.concatenate([grad_W.ravel(), grad_H.ravel()])

    # Enforce non-negativity
    bounds = [(0, None)] * len(x0)

    print("\nStarting optimization...")

    result = minimize(
        fun=objective,
        x0=x0,
        method="L-BFGS-B",
        bounds=bounds,
        jac=gradient,
        options={"maxiter": max_iter, "disp": True},
    )

    W = result.x[:w_size].reshape(n_cytokines, n_cytokines)
    H = result.x[w_size:].reshape(n_cytokines, n_components)

    # Evaluate

    A_recon = W @ H

    recon_error = np.linalg.norm(A - A_recon, "fro")
    rel_error = recon_error / np.linalg.norm(A, "fro")

    # Statistics for W
    w_sparsity = np.sum(np.abs(W) < 1e-3) / W.size
    w_mean = np.abs(W).mean()
    w_max = np.abs(W).max()

    # Statistics for H
    h_sparsity = np.sum(np.abs(H) < 1e-3) / H.size
    h_mean = np.abs(H).mean()
    h_max = np.abs(H).max()

    print("\nOptimization complete:")
    print(f"  Success: {result.success}")
    print(f"  Iterations: {result.nit}")
    print(f"  Relative reconstruction error: {rel_error:.4%}")

    print("\n  W (cytokine interactions):")
    print(f"    Shape: {W.shape}")
    print(f"    Sparsity: {w_sparsity:.2%} (near-zero elements)")
    print(f"    Mean |W|: {w_mean:.4f}")
    print(f"    Max |W|: {w_max:.4f}")
    print(f"    Non-zeros: {np.sum(np.abs(W) > 1e-3)}/{W.size}")

    print("\n  H (effect patterns):")
    print(f"    Shape: {H.shape}")
    print(f"    Sparsity: {h_sparsity:.2%} (near-zero elements)")
    print(f"    Mean |H|: {h_mean:.4f}")
    print(f"    Max |H|: {h_max:.4f}")
    print(f"    Non-zeros: {np.sum(np.abs(H) > 1e-3)}/{H.size}")

    return W, H
