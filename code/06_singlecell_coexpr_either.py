"""
singlecell coexpr either — produces coexpr_either_malignant.parquet
Artifact version: bb3a699b-e8b4-4888-8127-9fae7c1b4501
Environment: census
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import os
os.environ['NUMBA_CACHE_DIR']=os.path.join(os.getcwd(),'.numba_cache')
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp

adata = ad.read_h5ad("inputs/tnbc_sc_norm.h5ad")
genes = adata.var_names.tolist()
X = adata.X.tocsc()
counts = adata.layers['counts'].tocsc()
coarse = adata.obs['coarse'].values

mal_mask = (coarse=='malignant')
n_mal = mal_mask.sum()
gi = {g:i for i,g in enumerate(genes)}

counts_csr2 = adata.layers['counts'].tocsr()
mal_idx = np.where(mal_mask)[0]

surf_genes_n = pd.read_parquet("inputs/surface_narrowed.parquet")
sgn = surf_genes_n['gene'].tolist() if 'gene' in surf_genes_n.columns else surf_genes_n.index.tolist()

mal_counts = counts_csr2[mal_idx]  # malignant x genes (counts)
pos = (mal_counts>0)  # boolean sparse
sgn_present = [g for g in sgn if g in gi]
pos_dense = np.asarray(pos[:,[gi[g] for g in sgn_present]].todense())  # n_mal x n_surf
posA = pos_dense.astype(np.float32)
n = posA.shape[0]
# co-occurrence matrix
both = (posA.T @ posA)/n           # frac both positive
pA = posA.mean(axis=0)             # frac each positive
# either = pA + pB - both
either = pA[:,None]+pA[None,:]-both
coexp_either = pd.DataFrame(either, index=sgn_present, columns=sgn_present)
coexp_either.to_parquet('coexpr_either_malignant.parquet')