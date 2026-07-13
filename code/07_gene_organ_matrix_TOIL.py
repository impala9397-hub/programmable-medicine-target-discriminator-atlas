"""
gene organ matrix TOIL — produces gene_organ_matrix.csv
Artifact version: 0fb4ee3c-b1ad-49f2-9398-9100adc31e03
Environment: python
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import pickle
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False

# Load dependencies
prot_scored = pd.read_csv("inputs/protease_matrix_TNBC.csv").set_index('gene')
sm = pd.read_parquet("inputs/surface_matrix_TNBC.parquet").set_index('gene')
sm = sm[~sm.index.duplicated()]
narrowed = pd.read_parquet("inputs/surface_narrowed.parquet").set_index('gene')

# Load TOIL expression data from pickle (pre-fetched from Xena)
d = pickle.load(open('/tmp/toil_expr.pkl', 'rb'))
g2f = d['g2f']; want = d['want']; gtex = set(d['gtex']); tnbc = set(d['tnbc'])
site_map = dict(zip(d['samples'], d['site']))
inv = {v: k for k, v in g2f.items()}
gene_order = [inv[f] for f in d['fields']]

X = pd.DataFrame(d['mat'], index=gene_order, columns=want).astype(float)
lin = np.power(2, X) - 0.001
lin[lin < 0] = 0

gtex_cols = [c for c in want if c in gtex]
sites = sorted(set(site_map[c] for c in gtex_cols if site_map[c]))
FULL = pd.DataFrame({s: lin[[c for c in gtex_cols if site_map[c] == s]].median(axis=1) for s in sites})
FULL['TNBC'] = lin[[c for c in want if c in tnbc]].median(axis=1)

# gene classification
prot_set = set(prot_scored.index)
all_prot = [g for g in gene_order if g in prot_set]
all_surf = [g for g in gene_order if g not in prot_set]

# apply activity proxy to proteases
FULLp = FULL.copy()
for g in all_prot:
    FULLp.loc[g] = FULL.loc[g] * (1 - 0.25 * prot_scored.loc[g, 'inhibitor_penalty'])

# order columns: normals alpha, Breast last, TNBC
norm_sites = [c for c in FULLp.columns if c not in ('Breast', 'TNBC')] + ['Breast']
Mh = FULLp[norm_sites + ['TNBC']].loc[all_surf + all_prot]

Mh.to_csv('gene_organ_matrix.csv')