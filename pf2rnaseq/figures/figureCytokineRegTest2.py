
"""FMS score for different regularization parameters"""
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.axes import Axes
from tensorly.cp_tensor import CPTensor
from tlviz.factor_tools import factor_match_score as fms
import wandb
from ..factorization import pf2, pf2_pca_r2x
from .common import getSetup, subplotLabel
from ..imports import import_cytokine

run = wandb.init(
    # Set the wandb entity where your project will be logged (generally your team name).
    entity="nbedanova-ucla",
    # Set the wandb project where this run will be logged.
    project="L1_params",
    # Track hyperparameters and run metadata.
    config={
        "Regularize": 'A and C',
        "Runs": 3,
    },
)



def makeFigure():
    ax, f = getSetup((6, 3), (1, 1))
    subplotLabel(ax)

    X = import_cytokine()
    lambdaList = [5e-6, 1e-5, 5e-5, 1e-4]
    rank=np.arange(1, 31, 5)
    plot_fms_diff_reg(X, ax[0], rank, lambdaList, runs=3)

    

    return f


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


def resample(data: anndata.AnnData) -> anndata.AnnData:
    """Bootstrapping dataset"""
    indices = np.random.randint(0, data.shape[0], size=(data.shape[0],))
    data = data[indices].copy()
    return data


def plot_fms_diff_reg(
    X: anndata.AnnData,
    ax: Axes,
    rank: int,
    regsList: list[float],
    runs: int,
):
    #Plots FMS when using different regularization params
    fmsLists = []

    
    
    for r in rank:
       
        for i in regsList:
            scores = []
            dataX,r2X = pf2(X, rank=r, random_state=0, doEmbedding=False,regParam=i, regularize_A=True, r2x=True)
            for j in range(0, runs, 1):
                

                sampledX = pf2(resample(X), rank=r, random_state=j, doEmbedding=False, regParam=i,regularize_A=True)
                
                fmsScore = calculateFMS(dataX, sampledX)
                scores.append(fmsScore)
            
            run.log({"fms": np.mean(scores),"regParam": i,"rank": r, "R2X": r2X})
    for r in rank:
        scores = []
        dataX,r2X = pf2(X, rank=r, random_state=0, doEmbedding=False, r2x=True)
        for j in range(0, runs, 1):
            sampledX = pf2(resample(X), rank=r, random_state=j, doEmbedding=False)
            
            fmsScore = calculateFMS(dataX, sampledX)
            scores.append(fmsScore)
        
        run.log({"fms": np.mean(scores),"regParam": 0,"rank": r, "R2X": r2X})
        
        

   
run.finish()

