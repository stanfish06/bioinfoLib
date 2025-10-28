#!/usr/bin/env python3
"""
Explore single-cell RNA-seq data in h5ad format
"""

import scanpy as sc
import numpy as np
import pandas as pd

# Load the data
print("Loading adata_hem.h5ad...")
adata = sc.read_h5ad('examples/adata_hem.h5ad')

print("\n" + "="*60)
print("DATASET OVERVIEW")
print("="*60)
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"Shape: {adata.shape}")

print("\n" + "="*60)
print("OBSERVATION (CELL) METADATA")
print("="*60)
print(f"Available columns: {list(adata.obs.columns)}")
print(f"\nFirst few rows:")
print(adata.obs.head())

print("\n" + "="*60)
print("VARIABLE (GENE) METADATA")
print("="*60)
print(f"Available columns: {list(adata.var.columns)}")
print(f"\nFirst few rows:")
print(adata.var.head())

print("\n" + "="*60)
print("LAYERS")
print("="*60)
print(f"Available layers: {list(adata.layers.keys())}")

print("\n" + "="*60)
print("UNSTRUCTURED ANNOTATIONS")
print("="*60)
print(f"Available uns keys: {list(adata.uns.keys())}")

print("\n" + "="*60)
print("OBSM (Cell embeddings/dimensionality reductions)")
print("="*60)
print(f"Available obsm keys: {list(adata.obsm.keys())}")

print("\n" + "="*60)
print("VARM (Gene embeddings)")
print("="*60)
print(f"Available varm keys: {list(adata.varm.keys())}")

# Calculate basic QC metrics if not already present
print("\n" + "="*60)
print("BASIC QC METRICS")
print("="*60)

if 'n_genes_by_counts' not in adata.obs.columns:
    print("Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

print(f"\nCounts per cell:")
print(f"  Mean: {adata.obs['total_counts'].mean():.2f}")
print(f"  Median: {adata.obs['total_counts'].median():.2f}")
print(f"  Min: {adata.obs['total_counts'].min():.2f}")
print(f"  Max: {adata.obs['total_counts'].max():.2f}")

print(f"\nGenes detected per cell:")
print(f"  Mean: {adata.obs['n_genes_by_counts'].mean():.2f}")
print(f"  Median: {adata.obs['n_genes_by_counts'].median():.2f}")
print(f"  Min: {adata.obs['n_genes_by_counts'].min():.2f}")
print(f"  Max: {adata.obs['n_genes_by_counts'].max():.2f}")

# Check for mitochondrial genes
mt_genes = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
print(f"\nMitochondrial genes detected: {mt_genes.sum()}")

if mt_genes.sum() > 0 and 'pct_counts_mt' not in adata.obs.columns:
    adata.var['mt'] = mt_genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    print(f"Mitochondrial % per cell:")
    print(f"  Mean: {adata.obs['pct_counts_mt'].mean():.2f}%")
    print(f"  Median: {adata.obs['pct_counts_mt'].median():.2f}%")
elif 'pct_counts_mt' in adata.obs.columns:
    print(f"Mitochondrial % per cell:")
    print(f"  Mean: {adata.obs['pct_counts_mt'].mean():.2f}%")
    print(f"  Median: {adata.obs['pct_counts_mt'].median():.2f}%")

print("\n" + "="*60)
print("ANALYSIS COMPLETE")
print("="*60)
