#!/usr/bin/env python3
"""
Differential expression analysis for HSC-MPPs
"""

import scanpy as sc
import pandas as pd
import numpy as np

# Load the data
print("Loading adata_hem.h5ad...")
adata = sc.read_h5ad('examples/adata_hem.h5ad')

print(f"\nDataset: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"\nUsing 'Cluster' column for cell type annotations")
print(f"HSC-MPPs: {(adata.obs['Cluster'] == 'HSC-MPPs').sum()} cells")

# Check if DE results already exist
if 'rank_genes' in adata.uns:
    print("\nFound existing differential expression results in adata.uns['rank_genes']")
    print("Extracting HSC-MPPs markers...")

    # Get the results
    de_results = adata.uns['rank_genes']

    # Check what groups are available
    if 'HSC-MPPs' in de_results['names'].dtype.names:
        print("\nTop 20 marker genes for HSC-MPPs:")
        print("="*80)

        # Extract results for HSC-MPPs
        genes = de_results['names']['HSC-MPPs'][:20]
        scores = de_results['scores']['HSC-MPPs'][:20]
        pvals = de_results['pvals']['HSC-MPPs'][:20]
        pvals_adj = de_results['pvals_adj']['HSC-MPPs'][:20]
        logfoldchanges = de_results['logfoldchanges']['HSC-MPPs'][:20]

        # Create a nice table
        results_df = pd.DataFrame({
            'Gene': genes,
            'Score': scores,
            'Log2FC': logfoldchanges,
            'P-value': pvals,
            'Adj P-value': pvals_adj
        })

        print(results_df.to_string(index=False))

        # Save to file
        results_df.to_csv('examples/HSC_MPPs_markers.csv', index=False)
        print("\n\nFull results saved to examples/HSC_MPPs_markers.csv")

    else:
        print(f"\nAvailable groups in DE results: {de_results['names'].dtype.names}")
        print("HSC-MPPs not found in existing results. Running new DE analysis...")
        run_new_de = True
else:
    print("\nNo existing DE results found. Running differential expression analysis...")
    run_new_de = True

# If we need to run new DE analysis
if 'run_new_de' in locals() and run_new_de:
    print("\nRunning Wilcoxon rank-sum test for HSC-MPPs vs rest...")

    # Use the counts layer for DE analysis
    if 'counts' in adata.layers:
        print("Using 'counts' layer for analysis")
        adata_de = adata.copy()
        adata_de.X = adata_de.layers['counts']
    else:
        print("Using main matrix (X) for analysis")
        adata_de = adata

    # Run differential expression
    sc.tl.rank_genes_groups(
        adata_de,
        groupby='Cluster',
        groups=['HSC-MPPs'],
        reference='rest',
        method='wilcoxon',
        key_added='rank_genes_hsc_mpps'
    )

    # Get results
    result = sc.get.rank_genes_groups_df(adata_de, group='HSC-MPPs', key='rank_genes_hsc_mpps')

    print("\nTop 20 marker genes for HSC-MPPs:")
    print("="*80)
    print(result.head(20).to_string(index=False))

    # Save results
    result.to_csv('examples/HSC_MPPs_markers.csv', index=False)
    print("\n\nFull results saved to examples/HSC_MPPs_markers.csv")

    # Also save top 50 for reference
    print(f"\nTotal significant genes (adj p-val < 0.05): {(result['pvals_adj'] < 0.05).sum()}")

print("\nAnalysis complete!")
