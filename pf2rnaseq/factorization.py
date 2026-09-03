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


def filter_components_by_cytokine_quantile(
    threshold: float,
    csv_path: str = "/home/nicoleb/Pf2-scRNAseq-1/pf2rnaseq/Data/donor_vs_cytokine_variability_sd.csv",
) -> np.ndarray:
    """Return 0-indexed components whose cytokine_quantile is at/above a threshold.

    Parameters
    ----------
    threshold : float
        Minimum cytokine_quantile value a component must have to be kept.
    csv_path : str
        Path to the CSV with columns component, donor_sd, cytokine_sd,
        donor_quantile, cytokine_quantile.

    Returns
    -------
    np.ndarray
        0-indexed component numbers with cytokine_quantile >= threshold.
    """
    df = pd.read_csv(csv_path)
    return df.loc[df["cytokine_quantile"] >= threshold, "component"].to_numpy()


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
        # Exclude diagonal of W from L1 penalty
        l1_W = alpha * np.sum(np.abs(W)) - alpha * np.diag(np.abs(W)).sum()
        l1_H = alpha * np.sum(np.abs(H))

        total_loss = mse + l1_W + l1_H

        iteration_counter[0] += 1
        if total_loss < best_loss[0]:
            best_loss[0] = total_loss

        if iteration_counter[0] % 10 == 0:
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
        grad_W = 2 * ((W @ H - A) @ H.T) + alpha * np.sign(W) - np.diag(alpha * np.sign(np.diag(W)))

        # ===== Gradient w.r.t. H =====
        # 1. Reconstruction term: ∂/∂H [||A - WH||²] = 2(W^T @ error),  L1 penalty: ∂/∂H [α||H||₁] = α * sign(H)
        grad_H = 2 * (W.T @ (W @ H - A)) + alpha * np.sign(H)

        return np.concatenate([grad_W.ravel(), grad_H.ravel()])

    print("\nStarting optimization...")

    result = minimize(
        fun=objective,
        x0=x0,
        method="L-BFGS-B",
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


def deconvolution_cytokine_admm(
    A: np.ndarray,
    alpha_h: float = 0.1,
    alpha_w: float = 0.01,
    rho_w_init: float | None = None,
    rho_h_init: float | None = None,
    max_iter: int = 10000,
    tol_abs: float = 1e-4,
    tol_rel: float = 1e-3,
    random_state: int = 1,
    adaptive_rho: bool = True,
    rho_bounds: tuple[float, float] = (1e-4, 1e4),
    non_negative_w: bool = True,
    non_negative_h: bool = True,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Decompose cytokine factor matrix using ADMM: A ≈ W @ H

    Parameters
    ----------
    A : np.ndarray
        Input matrix (n_cytokines, n_components)
    alpha_h : float
        L1 regularization for H
    alpha_w : float
        L1 regularization for W (off-diagonal only)
    rho_w_init, rho_h_init : float or None
        Initial ADMM penalty for the W- and H-subproblems. If None,
        computed via the spectral rule rho* = sqrt(lambda_min * lambda_max)
        of the relevant Gram matrix (H H^T for W-step, W^T W for H-step).
    max_iter : int
        Maximum iterations
    tol_abs, tol_rel : float
        Absolute and relative tolerance terms for the combined stopping
        criterion (Boyd et al. 2011, §3.3.1), applied per-block.
    random_state : int
        Random seed
    adaptive_rho : bool
        Whether to adaptively adjust rho_w and rho_h independently
    rho_bounds : tuple[float, float]
        (min, max) clip range for adaptive rho updates
    non_negative_w : bool
        If True, enforce W ≥ 0 (cytokines only activate, not inhibit)
    non_negative_h : bool
        If True, enforce H ≥ 0

    Returns
    -------
    Z_W : np.ndarray
        Cytokine interaction matrix (n_cytokines, n_cytokines)
    Z_H : np.ndarray
        Effect basis matrix (n_cytokines, n_components)
    history : dict
        Optimization history
    """
    n_cytokines, n_components = A.shape
    np.random.seed(random_state)

    # Initialize
    W = np.random.rand(n_cytokines, n_cytokines) * 0.1 + np.eye(n_cytokines)
    H = np.random.rand(n_cytokines, n_components) * np.mean(np.abs(A))
    Z_W = W.copy()
    Z_H = H.copy()
    U_W = np.zeros_like(W)
    U_H = np.zeros_like(H)

    off_diag_mask = ~np.eye(n_cytokines, dtype=bool)
    rho_min, rho_max = rho_bounds

    def soft_threshold(X, threshold):
        return np.sign(X) * np.maximum(np.abs(X) - threshold, 0)

    def spectral_rho(M):
        """Geometric mean of the extreme eigenvalues of M @ M.T."""
        eigs = np.linalg.eigvalsh(M @ M.T)
        eigs = eigs[eigs > 1e-12]
        if eigs.size == 0:
            return 1.0
        return float(np.sqrt(eigs.min() * eigs.max()))

    # --- Initialize rho_w, rho_h independently ---
    rho_w = rho_w_init if rho_w_init is not None else spectral_rho(H)
    rho_h = rho_h_init if rho_h_init is not None else spectral_rho(W)
    rho_w = np.clip(rho_w, rho_min, rho_max)
    rho_h = np.clip(rho_h, rho_min, rho_max)

    def update_W(H, Z_W, U_W, rho_w):
        lhs = H @ H.T + rho_w * np.eye(n_cytokines)
        rhs = A @ H.T + rho_w * (Z_W - U_W)
        return np.linalg.solve(lhs, rhs.T).T

    def update_H(W, Z_H, U_H, rho_h):
        W_TW = W.T @ W
        W_TA = W.T @ A
        lhs = W_TW + rho_h * np.eye(n_cytokines)
        rhs = W_TA + rho_h * (Z_H - U_H)
        return np.linalg.solve(lhs, rhs)

    def update_Z_W(W, U_W, alpha, rho_w):
        Z_W_new = soft_threshold(W + U_W, alpha / rho_w)
        if non_negative_w:
            np.maximum(Z_W_new, 0, out=Z_W_new)
        np.fill_diagonal(Z_W_new, 1.0)
        return Z_W_new

    def update_Z_H(H, U_H, alpha, rho_h):
        Z_H_new = soft_threshold(H + U_H, alpha / rho_h)
        if non_negative_h:
            np.maximum(Z_H_new, 0, out=Z_H_new)
        return Z_H_new

    history = {
        "objective": [],
        "primal_residual_w": [],
        "dual_residual_w": [],
        "primal_residual_h": [],
        "dual_residual_h": [],
        "rho_w": [],
        "rho_h": [],
        "w_sparsity": [],
        "h_sparsity": [],
    }

    print("Cytokine deconvolution with ADMM:")
    print(f"  A shape: {A.shape}")
    print(f"  Alpha_W: {alpha_w}, Alpha_H: {alpha_h}")
    print(f"  Initial rho_w: {rho_w:.4g}, rho_h: {rho_h:.4g}")
    print(f"  Tolerance (abs, rel): ({tol_abs:.1e}, {tol_rel:.1e})")
    print(f"  Non-negative W: {non_negative_w}")
    print(f"  Non-negative H: {non_negative_h}")
    print("\nStarting ADMM iterations...")

    for iteration in range(max_iter):
        Z_W_old = Z_W.copy()
        Z_H_old = Z_H.copy()

        # ADMM updates
        W = update_W(H, Z_W, U_W, rho_w)
        H = update_H(W, Z_H, U_H, rho_h)
        Z_W = update_Z_W(W, U_W, alpha_w, rho_w)
        Z_H = update_Z_H(H, U_H, alpha_h, rho_h)
        U_W = U_W + (W - Z_W)
        U_H = U_H + (H - Z_H)

        # Per-block residuals, normalized by entry count
        r_w = np.linalg.norm(W - Z_W) / np.sqrt(W.size)
        s_w = np.linalg.norm(rho_w * (Z_W - Z_W_old)) / np.sqrt(W.size)
        r_h = np.linalg.norm(H - Z_H) / np.sqrt(H.size)
        s_h = np.linalg.norm(rho_h * (Z_H - Z_H_old)) / np.sqrt(H.size)

        # Combined absolute + relative stopping criteria (Boyd et al. 2011, §3.3.1)
        eps_pri_w = tol_abs + tol_rel * max(
            np.linalg.norm(W) / np.sqrt(W.size), np.linalg.norm(Z_W) / np.sqrt(W.size)
        )
        eps_dual_w = tol_abs + tol_rel * np.linalg.norm(rho_w * U_W) / np.sqrt(W.size)
        eps_pri_h = tol_abs + tol_rel * max(
            np.linalg.norm(H) / np.sqrt(H.size), np.linalg.norm(Z_H) / np.sqrt(H.size)
        )
        eps_dual_h = tol_abs + tol_rel * np.linalg.norm(rho_h * U_H) / np.sqrt(H.size)

        # Compute objective
        recon_error = np.sum((A - W @ H) ** 2)
        l1_W = alpha_w * np.sum(np.abs(Z_W[off_diag_mask]))
        l1_H = alpha_h * np.sum(np.abs(Z_H))
        objective = recon_error + l1_W + l1_H

        # Track sparsity
        w_sparsity = np.sum(np.abs(Z_W[off_diag_mask]) < 1e-3) / np.sum(off_diag_mask)
        h_sparsity = np.sum(np.abs(Z_H) < 1e-3) / Z_H.size

        # Store history
        history["objective"].append(objective)
        history["primal_residual_w"].append(r_w)
        history["dual_residual_w"].append(s_w)
        history["primal_residual_h"].append(r_h)
        history["dual_residual_h"].append(s_h)
        history["rho_w"].append(rho_w)
        history["rho_h"].append(rho_h)
        history["w_sparsity"].append(w_sparsity)
        history["h_sparsity"].append(h_sparsity)

        if iteration % 100 == 0 or iteration < 10:
            print(
                f"  Iter {iteration:4d}: Obj={objective:.4e}, "
                f"r_w={r_w:.3e}, s_w={s_w:.3e}, ρ_w={rho_w:.3g} | "
                f"r_h={r_h:.3e}, s_h={s_h:.3e}, ρ_h={rho_h:.3g}"
            )

        # Adaptive rho update — W and H blocks handled independently
        if adaptive_rho and iteration > 0:
            if r_w > 10 * s_w:
                rho_w = np.clip(rho_w * 2, rho_min, rho_max)
                U_W = U_W / 2
            elif s_w > 10 * r_w:
                rho_w = np.clip(rho_w / 2, rho_min, rho_max)
                U_W = U_W * 2

            if r_h > 10 * s_h:
                rho_h = np.clip(rho_h * 2, rho_min, rho_max)
                U_H = U_H / 2
            elif s_h > 10 * r_h:
                rho_h = np.clip(rho_h / 2, rho_min, rho_max)
                U_H = U_H * 2

        # Convergence check — both blocks must satisfy their combined criteria
        if (
            r_w < eps_pri_w
            and s_w < eps_dual_w
            and r_h < eps_pri_h
            and s_h < eps_dual_h
        ):
            print(f"\n✓ Converged at iteration {iteration}")
            print(
                f"  W block: r={r_w:.4e} < {eps_pri_w:.4e}, s={s_w:.4e} < {eps_dual_w:.4e}"
            )
            print(
                f"  H block: r={r_h:.4e} < {eps_pri_h:.4e}, s={s_h:.4e} < {eps_dual_h:.4e}"
            )
            break

    # Final statistics
    A_recon = W @ H
    rel_error = np.linalg.norm(A - A_recon, "fro") / np.linalg.norm(A, "fro")

    w_sparsity = np.sum(np.abs(Z_W[off_diag_mask]) < 1e-3) / np.sum(off_diag_mask)
    h_sparsity = np.sum(np.abs(Z_H) < 1e-3) / Z_H.size

    print("\nOptimization complete:")
    print(f"  Iterations: {iteration + 1}/{max_iter}")
    print(f"  Relative reconstruction error: {rel_error:.4%}")
    print(f"  Final rho_w: {rho_w:.4g}, rho_h: {rho_h:.4g}")

    print("\n  W (cytokine interactions):")
    print(f"    Off-diagonal sparsity: {w_sparsity:.2%}")
    print(f"    Off-diagonal non-zeros: {np.sum(np.abs(Z_W[off_diag_mask]) > 1e-3)}")
    print(f"    Mean |W_offdiag|: {np.abs(Z_W[off_diag_mask]).mean():.4f}")
    print(f"    Min value: {W.min():.4f}")
    print(f"    Max value: {W.max():.4f}")
    print(f"    Diagonal: all 1.0 (constrained)")

    print("\n  H (effect patterns):")
    print(f"    Sparsity: {h_sparsity:.2%}")
    print(f"    Non-zeros: {np.sum(np.abs(Z_H) > 1e-3)}/{Z_H.size}")
    print(f"    Mean |H|: {np.abs(Z_H).mean():.4f}")
    print(f"    Min value: {H.min():.4f}")
    print(f"    Max value: {H.max():.4f}")
    print(f"    Negative values: {np.sum(H < 0)} ({100 * np.sum(H < 0) / H.size:.1f}%)")

    return Z_W, Z_H, history

def correct_conditions_from_npy(
    A: np.ndarray,
    cyt_order: np.ndarray,
    h5ad_path: str = "/opt/data/Parse_10M_PBMC_cytokines.h5ad",
    cytokine_col: str = "cytokine",
) -> np.ndarray:
    """Correct cytokine factor matrix by total read depth per condition.

    Mirrors correct_conditions() but works with a standalone numpy factor matrix
    and reads the h5ad in backed mode, using the obs cytokine column to group
    cells without loading X.X into memory.

    Parameters
    ----------
    A : np.ndarray
        Cytokine factor matrix, shape (n_cytokines, n_components).
        Rows must correspond to cytokines in cyt_order.
    cyt_order : np.ndarray
        Cytokine names in the same row order as A.
    h5ad_path : str
        Path to the h5ad file opened in backed mode.
    cytokine_col : str
        Column in obs that identifies each cell's cytokine condition.

    Returns
    -------
    np.ndarray
        Corrected factor matrix, same shape as A.
    """
    import anndata as an

    X = an.read_h5ad(h5ad_path, backed="r")
    obs = X.obs

    if "total_counts" in obs.columns:
        per_cell = obs["total_counts"]
    else:
        per_cell = pd.Series(np.ones(len(obs), dtype=np.float64), index=obs.index)

    # Sum read counts (or cell counts) per cytokine condition
    cyt_totals = per_cell.groupby(obs[cytokine_col]).sum()
    X.file.close()

    # Align to the row order of A using cyt_order
    counts = cyt_totals.reindex(cyt_order).to_numpy(dtype=np.float64).reshape(-1, 1)

    A = A.copy()
    min_val = np.min(A)
    if min_val < 0:
        A = A + abs(min_val) + 1e-10

    cond_mean = gmean(A, axis=1)

    lr = LinearRegression()
    lr.fit(counts, cond_mean.reshape(-1, 1))
    counts_correct = lr.predict(counts)

    return A / counts_correct

