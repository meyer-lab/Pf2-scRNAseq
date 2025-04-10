"""
exports csv of top 30 and bottom 30 genes per component
"""

from anndata import read_h5ad
from imports import import_cytokine,import_pf2Cytokine30 
import pandas as pd
import numpy as np
import os
#from common import getSetup, subplotLabel






def bot_top_genes(X, cmp, geneAmount=5):
    """Saves most pos/negatively genes"""
    df = pd.DataFrame(
        data=X.varm["Pf2_C"][:, cmp - 1], index=X.var_names, columns=["Component"]
    )
    df = df.reset_index(names="Gene")
    df = df.sort_values(by="Component")

    top = df.iloc[-geneAmount:, 0].values
    bot = df.iloc[:geneAmount, 0].values
    all_genes = np.concatenate([bot, top])

    return all_genes
def export_genes_to_csv(X, output_path="component_genes.csv", components=range(1, 31), geneAmount=30):
    """
    Exports top and bottom genes for each component to a single CSV file.
    For each component:
    - One row for bottom genes
    - One row for bottom weights
    - One row for top genes
    - One row for top weights
    
    Parameters:
    X: AnnData object
    output_path: Path to save the CSV file
    components: Range of components to process
    geneAmount: Number of top/bottom genes to extract
    """
    # Create a DataFrame to store results
    results = []
    
    # Process each component
    for comp in components:
        # Get the component data
        df = pd.DataFrame(
            data=X.varm["Pf2_C"][:, comp - 1], index=X.var_names, columns=["Weight"]
        )
        df = df.reset_index(names="Gene")
        df = df.sort_values(by="Weight")
        
        # Get bottom genes with their weights
        bottom_genes = df.iloc[:geneAmount, 0].values
        bottom_weights = df.iloc[:geneAmount, 1].values
        
        # Add bottom genes as a row
        results.append({
            'Component': comp,
            'Type': 'Bottom',
            'Content': 'Genes',
            'Values': ','.join(bottom_genes),
        })
        
        # Add bottom weights as a separate row
        results.append({
            'Component': comp,
            'Type': 'Bottom',
            'Content': 'Weights',
            'Values': ','.join(map(str, bottom_weights)),
        })
        
        # Get top genes with their weights
        top_genes = df.iloc[-geneAmount:, 0].values
        top_weights = df.iloc[-geneAmount:, 1].values
        
        # Add top genes as a row
        results.append({
            'Component': comp,
            'Type': 'Top',
            'Content': 'Genes',
            'Values': ','.join(top_genes),
        })
        
        # Add top weights as a separate row
        results.append({
            'Component': comp,
            'Type': 'Top',
            'Content': 'Weights',
            'Values': ','.join(map(str, top_weights)),
        })
    
    # Convert to DataFrame and save to CSV
    results_df = pd.DataFrame(results)
    
    # Note: Change file extension to .csv since we're saving a CSV file
    if output_path.endswith('.xlsx'):
        output_path = output_path.replace('.xlsx', '.csv')
        
    results_df.to_csv(output_path, index=False)
    print(f"Exported to {output_path}")
    
    return output_path


X = read_h5ad("/home/nicoleb/Cytokine_Pf2_annotated_NB_031725.h5ad")
export_genes_to_csv(X, output_path="component_genes.xlsx", components=range(1, 31), geneAmount=30)    


