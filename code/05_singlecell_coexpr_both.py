"""
singlecell coexpr both — produces coexpr_both_malignant.parquet
Artifact version: 8357ae32-6ad5-47d0-b4f3-24921c8b75e5
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
cats = pd.Index(sorted(pd.unique(coarse)))

gi = {g:i for i,g in enumerate(genes)}

surf_genes_n = pd.read_parquet("inputs/surface_narrowed.parquet")
sgn = surf_genes_n['gene'].tolist() if 'gene' in surf_genes_n.columns else surf_genes_n.index.tolist()

mal_mask = (coarse=='malignant')
n_mal = mal_mask.sum()
counts_csr2 = adata.layers['counts'].tocsr()
mal_idx = np.where(mal_mask)[0]

sgn_present = [g for g in sgn if g in gi]
mal_counts = counts_csr2[mal_idx]
pos = (mal_counts>0)
pos_dense = np.asarray(pos[:,[gi[g] for g in sgn_present]].todense())
posA = pos_dense.astype(np.float32)
n = posA.shape[0]
both = (posA.T @ posA)/n
coexp_both = pd.DataFrame(both, index=sgn_present, columns=sgn_present)
coexp_both.to_parquet('coexpr_both_malignant.parquet')