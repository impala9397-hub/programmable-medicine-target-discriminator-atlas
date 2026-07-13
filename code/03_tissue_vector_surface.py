"""
tissue vector surface — produces tissue_vec_surface.parquet
Artifact version: e5b3ae2f-7dfc-4482-8150-b9f1725a4183
Environment: python
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import pandas as pd
import numpy as np

gtex_prior = pd.read_csv("inputs/STEP2_gtex_normal_expression.csv")

surf_narrow = pd.read_parquet("surface_narrowed.parquet")

surf_genes_n = surf_narrow.index.tolist() if surf_narrow.index.name == 'gene' else surf_narrow['gene'].tolist()

gtex_surf_full = gtex_prior.set_index('geneSymbol')
tissue_cols = list(gtex_surf_full.columns)

surf_tissue_vec = gtex_surf_full.reindex([g for g in surf_genes_n if g in gtex_surf_full.index])

surf_tissue_vec.to_parquet('tissue_vec_surface.parquet')