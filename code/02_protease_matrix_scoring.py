"""
protease matrix scoring — produces protease_matrix_TNBC.csv
Artifact version: d86c0f9c-6c31-4452-870a-c788b9853081
Environment: python
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import pandas as pd
import numpy as np
import pickle
import json

# Load protease matrix (initial scored version)
prot = pd.read_csv("inputs/protease_matrix_TNBC_scored.csv").set_index('gene')

# Load TCGA Basal RSEM matrix (genes x samples)
rsem = pd.read_parquet("inputs/tcga_basal_rsem.parquet")
rsem = rsem[~rsem.index.duplicated()]

# Load GTEx protease normal expression
gp = pd.read_parquet("inputs/gtex_protease_normal.parquet")

# Inhibitor families mapped to their inhibitor genes
inh_map = {
    'MMP': ['TIMP1', 'TIMP2', 'TIMP3', 'TIMP4'],
    'ADAM': ['TIMP1', 'TIMP2', 'TIMP3', 'TIMP4'],
    'ADAMTS': ['TIMP1', 'TIMP2', 'TIMP3', 'TIMP4'],
    'serine': ['SERPINE1', 'SERPINB2', 'SERPINA1', 'SERPINA3', 'SERPINB5'],
    'cysteine': ['CST3', 'CST4', 'CST6', 'CSTA', 'CSTB'],
    'aspartic': [], 'glycosidase': [],
}

# Median TME expression per gene
med = rsem.median(axis=1)

# Inhibitor pool TME level per family (sum of median expression of its inhibitors)
fam_inh = {f: float(med.reindex(genes).fillna(0).sum()) for f, genes in inh_map.items()}

# Per-protease: ratio = inhibitor_pool / (protease_expr + inhibitor_pool)
prot['inh_tme'] = prot['family'].map(fam_inh)
prot['ratio'] = prot['inh_tme'] / (prot['spec_tme_tpm'] + prot['inh_tme'])

# Rebuild spec_score_ACTIVITY_PROXY with ratio-based penalty
from scipy.stats import rankdata

base = (0.55 * (rankdata(prot['spec_ratio_log2']) / len(prot)) +
        0.45 * (1 - rankdata(np.log1p(prot['tox_load'])) / len(prot)))

prot['inhibitor_penalty_minmax_OLD'] = prot['inhibitor_penalty']
prot['inhibitor_penalty'] = prot['ratio'].round(3)
prot['spec_score_ACTIVITY_PROXY'] = (base * (1 - 0.25 * prot['inhibitor_penalty'])).round(4)
prot['inhibitor_penalty_method'] = 'inhibitor/(protease+inhibitor) TME median RSEM ratio'

prot = prot.drop(columns=['ratio', 'inh_tme'])

# Recompute inter_cov_pct_patients at 250 TPM
PROT_THRESH = 250

def cov_at(g, t):
    return float((rsem.loc[g] > t).mean() * 100) if g in rsem.index else np.nan

prot['inter_cov_pct_patients_OLD_50'] = prot['inter_cov_pct_patients']
prot['inter_cov_pct_patients'] = prot.index.map(lambda g: round(cov_at(g, PROT_THRESH), 1))
prot['inter_cov_threshold_tpm'] = PROT_THRESH

# Add TOIL-based columns (tumor/normal ratio for ADC suitability)
# Load gene_organ_matrix (TOIL comparable TPM, genes x tissues+TNBC)
gene_organ = pd.read_csv('gene_organ_matrix.csv', index_col=0)

norm_cols = [c for c in gene_organ.columns if c != 'TNBC']

prot['toil_tnbc_tpm'] = np.nan
prot['toil_normal_max_tpm'] = np.nan
prot['toil_tumor_norm_ratio'] = np.nan
prot['toil_n_normal_active_gt10'] = np.nan
prot['toil_gini_normal'] = np.nan

for p in prot.index:
    if p in gene_organ.index:
        row = gene_organ.loc[p]
        v = row[norm_cols].astype(float)
        nmax = v.max()
        tnbc = row['TNBC']
        prot.loc[p, 'toil_tnbc_tpm'] = round(tnbc, 2)
        prot.loc[p, 'toil_normal_max_tpm'] = round(nmax, 2)
        prot.loc[p, 'toil_tumor_norm_ratio'] = round(tnbc / (nmax + 1e-9), 3)
        prot.loc[p, 'toil_n_normal_active_gt10'] = int((v > 10).sum())
        s = np.sort(v.values)
        n = len(s)
        gini = (2 * np.sum((np.arange(1, n + 1)) * s) / (n * s.sum() + 1e-9)) - (n + 1) / n
        prot.loc[p, 'toil_gini_normal'] = round(gini, 3)

# Housekeeping flag: broad (n_gt10 > 8) OR uniform (gini < 0.4)
def is_housekeeping(p):
    if pd.isna(prot.loc[p, 'toil_n_normal_active_gt10']):
        return None
    n10 = prot.loc[p, 'toil_n_normal_active_gt10']
    g = prot.loc[p, 'toil_gini_normal']
    return bool((n10 > 8) or (g < 0.4))

prot['toil_housekeeping_flag'] = [is_housekeeping(p) for p in prot.index]

prot.to_csv('protease_matrix_TNBC.csv')
print("saved protease_matrix_TNBC.csv, shape:", prot.shape)