"""
Hyperparameter sweep for Pf2 using Weights & Biases
Optimizing rank and regularization parameter
"""

import numpy as np
import wandb
from factorization import pf2
from imports import import_cytokine
from tensorly.cp_tensor import CPTensor
from tlviz.factor_tools import factor_match_score as fms

ranks = np.arange(1, 31)
# Define the sweep configuration
sweep_config = {
    "method": "grid",  # grid search for thorough exploration
    "metric": {
        "name": "fms",  # optimize for factor match score
        "goal": "maximize",  # we want to maximize factor stability
    },
    "parameters": {
        "rank": {
            "values": ranks  # Different component numbers to test
        },
        "regParam": {
            "values": [
                0.0,
                1e-6,
                1e-5,
                5e-5,
                1e-4,
            ]  # Different L1 regularization strengths
        },
    },
}


def resample(data):
    """Bootstrapping dataset"""
    indices = np.random.randint(0, data.shape[0], size=(data.shape[0],))
    return data[indices].copy()


def calculateFMS(A, B):
    """Calculates FMS between 2 factorizations"""
    A_factors = [A.uns["Pf2_A"], A.uns["Pf2_B"], A.varm["Pf2_C"]]
    A_CP = CPTensor((A.uns["Pf2_weights"], A_factors))

    B_factors = [B.uns["Pf2_A"], B.uns["Pf2_B"], B.varm["Pf2_C"]]
    B_CP = CPTensor((B.uns["Pf2_weights"], B_factors))

    return fms(A_CP, B_CP, consider_weights=False, skip_mode=1)


def calculate_sparsity(matrix, threshold=1e-6):
    """Calculate sparsity (proportion of near-zero elements)"""
    total_elements = matrix.size
    near_zero_elements = np.sum(np.abs(matrix) < threshold)
    return near_zero_elements / total_elements


def train():
    """Main training function for wandb sweep"""
    # Initialize a new wandb run
    with wandb.init() as run:
        # Get parameters from wandb
        config = wandb.config

        # Load data (do this once per run to save time)
        X = import_cytokine()
        print(f"Running with rank={config.rank}, regParam={config.regParam}")

        # Set number of bootstrap samples
        n_bootstrap = 3

        # Run base factorization with current parameters
        base_model, r2x = pf2(
            X,
            rank=config.rank,
            random_state=42,
            doEmbedding=False,
            regParam=config.regParam,
            r2x=True,
        )

        sparsity_C = calculate_sparsity(base_model.varm["Pf2_C"])

        # Log R2X and sparsity metrics
        wandb.log({"r2x": r2x, "sparsity_C": sparsity_C})

        # Calculate FMS across bootstrap samples
        fms_scores = []
        for i in range(n_bootstrap):
            # Create bootstrap sample
            bootstrap_data = resample(X)

            # Run factorization on bootstrap sample
            bootstrap_model = pf2(
                bootstrap_data,
                rank=config.rank,
                random_state=i,
                doEmbedding=False,
                regParam=config.regParam,
            )

            # Calculate FMS between base model and bootstrap model
            fms_score = calculateFMS(base_model, bootstrap_model)
            fms_scores.append(fms_score)

            # Log individual bootstrap FMS
            wandb.log({f"fms_bootstrap_{i}": fms_score})

        # Calculate and log average FMS
        avg_fms = np.mean(fms_scores)
        wandb.log({"fms": avg_fms})

        print(
            f"Completed run: rank={config.rank}, regParam={config.regParam}, R2X={r2x:.4f}, FMS={avg_fms:.4f}"
        )


if __name__ == "__main__":
    # Initialize wandb
    wandb.login()

    # Create the sweep
    sweep_id = wandb.sweep(sweep_config, project="Pf2_parameter_optimization2")

    # Run the sweep
    wandb.agent(
        sweep_id, function=train, count=None
    )  # Set count if you want to limit runs
