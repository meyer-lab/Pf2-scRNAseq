"""
Test the cytokine deconvolution method.
"""

import numpy as np
import pytest

from ..factorization import deconvolution_cytokine_admm


def test_deconvolution_cytokine_admm_sparse():
    """
    Test deconvolution_cytokine_admm with sparse ground truth matrices.

    This test generates sparse W (cytokine interaction) and H (effect basis) matrices,
    computes A = W @ H, and verifies that the deconvolution recovers the structure.
    """
    np.random.seed(42)

    # Dimensions
    n_cytokines = 8
    n_components = 12

    # Generate sparse ground truth W (cytokine interaction matrix)
    # W should have 1s on diagonal and sparse off-diagonal elements
    W_true = np.eye(n_cytokines)

    # Add sparse off-diagonal interactions (only 20% of off-diagonal elements)
    off_diag_mask = ~np.eye(n_cytokines, dtype=bool)
    n_off_diag = np.sum(off_diag_mask)
    n_nonzero_w = int(0.2 * n_off_diag)

    # Randomly select positions for non-zero off-diagonal elements
    off_diag_positions = np.where(off_diag_mask)
    nonzero_indices = np.random.choice(n_off_diag, n_nonzero_w, replace=False)

    for idx in nonzero_indices:
        i, j = off_diag_positions[0][idx], off_diag_positions[1][idx]
        # Use small positive values for cytokine interactions
        W_true[i, j] = np.random.uniform(0.1, 0.5)

    # Generate sparse ground truth H (effect basis matrix)
    # H should have about 30% non-zero elements
    H_true = np.zeros((n_cytokines, n_components))
    n_nonzero_h = int(0.3 * n_cytokines * n_components)

    for _ in range(n_nonzero_h):
        i = np.random.randint(0, n_cytokines)
        j = np.random.randint(0, n_components)
        # H can have both positive and negative values
        H_true[i, j] = np.random.uniform(-2.0, 2.0)

    # Compute the observed matrix A
    A = W_true @ H_true

    # Add small noise
    noise_level = 0.01
    A_noisy = A + noise_level * np.random.randn(n_cytokines, n_components)

    # Run deconvolution
    W_recovered, H_recovered, history = deconvolution_cytokine_admm(
        A_noisy,
        alpha_h=0.1,
        alpha_w=0.05,
        rho=1.0,
        max_iter=1000,
        tol=1e-6,
        random_state=42,
        adaptive_rho=True,
        non_negative_w=True,
    )

    # Verify shapes
    assert W_recovered.shape == (n_cytokines, n_cytokines)
    assert H_recovered.shape == (n_cytokines, n_components)

    # Verify diagonal of W is constrained to 1
    np.testing.assert_allclose(np.diag(W_recovered), np.ones(n_cytokines), atol=1e-10)

    # Verify non-negativity of W
    assert np.all(W_recovered >= -1e-10), "W should be non-negative"

    # Verify reconstruction quality
    A_reconstructed = W_recovered @ H_recovered
    reconstruction_error = np.linalg.norm(
        A_noisy - A_reconstructed, "fro"
    ) / np.linalg.norm(A_noisy, "fro")
    assert reconstruction_error < 0.1, (
        f"Reconstruction error too high: {reconstruction_error}"
    )

    # Verify sparsity of W (off-diagonal should be sparse)
    w_sparsity = np.sum(np.abs(W_recovered[off_diag_mask]) < 1e-3) / np.sum(
        off_diag_mask
    )
    assert w_sparsity > 0.5, f"W should be sparse, but sparsity is only {w_sparsity}"

    # Verify sparsity of H
    h_sparsity = np.sum(np.abs(H_recovered) < 1e-3) / H_recovered.size
    assert h_sparsity > 0.3, f"H should be sparse, but sparsity is only {h_sparsity}"

    # Verify history contains expected keys
    assert "objective" in history
    assert "primal_residual" in history
    assert "dual_residual" in history
    assert "rho" in history
    assert "w_sparsity" in history
    assert "h_sparsity" in history

    # Verify objective decreases (generally)
    assert len(history["objective"]) > 0
    # Check that final objective is lower than initial (with some tolerance for fluctuations)
    initial_obj = history["objective"][0]
    final_obj = history["objective"][-1]
    assert final_obj < initial_obj * 1.1, "Objective should generally decrease"

    print("\nTest passed!")
    print(f"Reconstruction error: {reconstruction_error:.4f}")
    print(f"W off-diagonal sparsity: {w_sparsity:.2%}")
    print(f"H sparsity: {h_sparsity:.2%}")
    print(f"Converged in {len(history['objective'])} iterations")


def test_deconvolution_cytokine_admm_small():
    """
    Test with a small problem to ensure basic functionality.
    """
    np.random.seed(999)

    n_cytokines = 3
    n_components = 5

    # Simple test matrix
    A = np.random.randn(n_cytokines, n_components)

    # Run with default parameters
    W, H, history = deconvolution_cytokine_admm(
        A, max_iter=100, tol=1e-6, random_state=999
    )

    # Basic checks
    assert W.shape == (n_cytokines, n_cytokines)
    assert H.shape == (n_cytokines, n_components)
    assert len(history["objective"]) > 0

    # Verify diagonal constraint
    np.testing.assert_allclose(np.diag(W), np.ones(n_cytokines), atol=1e-10)

    print("\nSmall problem test passed!")
