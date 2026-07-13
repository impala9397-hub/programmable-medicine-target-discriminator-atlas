"""
logic gate synthesis — produces logic_gate_recommendations.csv
Artifact version: 98115bfa-31ef-4f6b-b1d7-e320b9c61e70
Environment: python
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import itertools
import numpy as np
import pandas as pd

coexp_both = pd.read_parquet("inputs/coexpr_both_malignant.parquet")
coexp_either = pd.read_parquet("inputs/coexpr_either_malignant.parquet")
ct_detect = pd.read_parquet("inputs/celltype_detect_fraction.parquet")
tissue_surf = pd.read_parquet("inputs/tissue_vec_surface.parquet")
tissue_prot = pd.read_parquet("inputs/tissue_vec_protease.parquet")
prot_scored = pd.read_csv("inputs/protease_matrix_TNBC_scored.csv").set_index('gene')

import pickle
import json

gtex = pd.read_csv("inputs/STEP2_gtex_normal_expression.csv", index_col=0)
gp_parquet = pd.read_parquet("inputs/gtex_protease_normal.parquet")
rsem = pd.read_parquet("inputs/tcga_basal_rsem.parquet")
rsem = rsem[~rsem.index.duplicated()]

# Rebuild FULLp (TOIL comparable matrix) is not available here, so we use the
# prot_scored as loaded (which already has inhibitor_penalty from the artifact).
# The final gate filter uses toil data computed in trace; we replicate with available data.

# Recompute inhibitor penalty as ratio using TCGA RSEM
inhibitors_map = {
    'MMP': ['TIMP1','TIMP2','TIMP3','TIMP4'],
    'ADAM': ['TIMP1','TIMP2','TIMP3','TIMP4'],
    'ADAMTS': ['TIMP1','TIMP2','TIMP3','TIMP4'],
    'serine': ['SERPINE1','SERPINB2','SERPINA1','SERPINA3','SERPINB5'],
    'cysteine': ['CST3','CST4','CST6','CSTA','CSTB'],
    'aspartic': [],
    'glycosidase': [],
}
med = rsem.median(axis=1)
fam_inh = {f: float(med.reindex(genes).fillna(0).sum()) for f, genes in inhibitors_map.items()}

prot_scored['inh_tme'] = prot_scored['family'].map(fam_inh)
prot_scored['ratio'] = prot_scored['inh_tme'] / (prot_scored['spec_tme_tpm'] + prot_scored['inh_tme'])

from scipy.stats import rankdata
p = prot_scored.copy()
base = 0.55 * (rankdata(p['spec_ratio_log2']) / len(p)) + 0.45 * (1 - rankdata(np.log1p(p['tox_load'])) / len(p))
p['inhibitor_penalty'] = p['ratio'].round(3)
p['spec_score_ACTIVITY_PROXY'] = (base * (1 - 0.25 * p['inhibitor_penalty'])).round(4)

# Recompute inter_cov at 250 TPM
PROT_THRESH = 250
def cov_at(g, t):
    return float((rsem.loc[g] > t).mean() * 100) if g in rsem.index else np.nan
p['inter_cov_pct_patients'] = p.index.map(lambda g: round(cov_at(g, PROT_THRESH), 1))
p['inter_cov_threshold_tpm'] = PROT_THRESH

prot_scored2 = p.copy()

# Compute TOIL-based breadth filter using GTEx (normal) for proteases
norm_cols_gtex = [c for c in gp_parquet.columns]
prot_scored2['toil_n_normal_active_gt10'] = np.nan
prot_scored2['toil_gini_normal'] = np.nan
for pg in prot_scored2.index:
    if pg in gp_parquet.index:
        v = gp_parquet.loc[pg].astype(float)
        pen = prot_scored2.loc[pg, 'inhibitor_penalty']
        va = v * (1 - 0.25 * pen)
        n10 = int((va > 10).sum())
        s = np.sort(va.values)
        n = len(s)
        gini = (2 * np.sum((np.arange(1, n+1)) * s) / (n * s.sum() + 1e-9)) - (n+1)/n
        prot_scored2.loc[pg, 'toil_n_normal_active_gt10'] = n10
        prot_scored2.loc[pg, 'toil_gini_normal'] = round(gini, 3)

def is_housekeeping(pg):
    n10 = prot_scored2.loc[pg, 'toil_n_normal_active_gt10']
    g = prot_scored2.loc[pg, 'toil_gini_normal']
    if pd.isna(n10) or pd.isna(g):
        return None
    return (n10 > 8) or (g < 0.4)

cov_spec_pass = prot_scored2[
    (prot_scored2['spec_score_ACTIVITY_PROXY'] >= 0.4) &
    (prot_scored2['inter_cov_pct_patients'] >= 70)
].index.tolist()

prot_ok_final = [pg for pg in cov_spec_pass if is_housekeeping(pg) == False]

# Setup logic engine
crit = {
    **{f'Brain_{x}': 3 for x in ['Cortex','Cerebellum','Hippocampus','Amygdala',
       'Frontal_Cortex_BA9','Anterior_cingulate_cortex_BA24','Hypothalamus',
       'Caudate_basal_ganglia','Nucleus_accumbens_basal_ganglia','Putamen_basal_ganglia',
       'Substantia_nigra','Spinal_cord_cervical_c-1','Cerebellar_Hemisphere']},
    'Heart_Left_Ventricle': 3, 'Heart_Atrial_Appendage': 3,
    'Liver': 3, 'Lung': 3, 'Kidney_Cortex': 3, 'Kidney_Medulla': 3,
    'Pancreas': 2.5, 'Whole_Blood': 2.5
}

def tcrit(t):
    return crit.get(t, 1.5)

LEAK = 10
GATE = 20
tissue_cols = list(tissue_surf.columns)
surf_cand = [g for g in coexp_both.index]
mal_cov = pd.Series(np.diag(coexp_both.values), index=coexp_both.index)

malignant_detect = ct_detect['malignant']
nonmal = [c for c in ct_detect.columns if c != 'malignant']
sc_spec = malignant_detect / (malignant_detect + ct_detect[nonmal].max(axis=1) + 1e-9)

surf_leak = {g: {t: tissue_surf.loc[g, t] for t in tissue_cols if g in tissue_surf.index and tissue_surf.loc[g, t] > LEAK}
             for g in surf_cand}
prot_gate = {pg: set(t for t in tissue_cols if tissue_prot.loc[pg, t] > GATE) for pg in tissue_prot.index}

def leak_burden(d):
    return sum(tcrit(t) * np.log1p(v) for t, v in d.items())

def or_leak(gs):
    m = {}
    for g in gs:
        for t, v in surf_leak.get(g, {}).items():
            m[t] = max(m.get(t, 0), v)
    return m

pairable = [g for g in surf_cand if mal_cov[g] >= 0.05]

recs = []
for g in pairable:
    lk = surf_leak.get(g, {})
    recs.append(dict(logic=g, type='single', inputs=1, coverage=mal_cov[g], escape_risk=1-mal_cov[g],
        specificity=sc_spec[g] if g in sc_spec.index else np.nan,
        leak_burden=leak_burden(lk), leak_tissues=";".join(sorted(lk)),
        n_crit_leak=sum(1 for t in lk if tcrit(t) >= 2.5)))

for a, b in itertools.combinations(pairable, 2):
    cov = coexp_either.loc[a, b]
    lk = or_leak([a, b])
    recs.append(dict(logic=f"{a} OR {b}", type='OR', inputs=2, coverage=cov, escape_risk=1-cov,
        specificity=np.nan, leak_burden=leak_burden(lk), leak_tissues=";".join(sorted(lk)),
        n_crit_leak=sum(1 for t in lk if tcrit(t) >= 2.5)))

for a, b in itertools.combinations(pairable, 2):
    cov = coexp_both.loc[a, b]
    if cov < 0.02:
        continue
    lk = {t: min(surf_leak.get(a, {}).get(t, 0), surf_leak.get(b, {}).get(t, 0))
          for t in set(surf_leak.get(a, {})) & set(surf_leak.get(b, {}))}
    recs.append(dict(logic=f"{a} AND {b}", type='AND', inputs=2, coverage=cov, escape_risk=1-cov,
        specificity=np.nan, leak_burden=leak_burden(lk), leak_tissues=";".join(sorted(lk)),
        n_crit_leak=sum(1 for t in lk if tcrit(t) >= 2.5)))

logic_df = pd.DataFrame(recs)

or_gates = logic_df[(logic_df['type'] == 'OR') & (logic_df['coverage'] >= 0.45)].copy()
grecs = []
for _, row in or_gates.iterrows():
    leak_t = set(row['leak_tissues'].split(';')) if row['leak_tissues'] else set()
    crit_leak = {t for t in leak_t if tcrit(t) >= 2.5}
    if not crit_leak:
        continue
    for pg in prot_ok_final:
        pgates = prot_gate.get(pg, set())
        rescued = crit_leak - pgates
        still = crit_leak & pgates
        if len(rescued) == 0:
            continue
        grecs.append(dict(
            logic=f"({row['logic']}) AND {pg}", type='OR_AND_protease', inputs=3,
            coverage=row['coverage'], escape_risk=row['escape_risk'], protease=pg,
            prot_spec=prot_scored2.loc[pg, 'spec_score_ACTIVITY_PROXY'],
            prot_cov=prot_scored2.loc[pg, 'inter_cov_pct_patients'],
            crit_leak_before=";".join(sorted(crit_leak)),
            crit_rescued=";".join(sorted(rescued)),
            crit_still_leak=";".join(sorted(still)),
            rescue_frac=len(rescued) / len(crit_leak)))

gate_df = pd.DataFrame(grecs)

def norm(s):
    s = pd.Series(s).astype(float)
    return (s - s.min()) / (s.max() - s.min() + 1e-9)

allL = pd.concat([logic_df, gate_df], ignore_index=True)
maxlb = logic_df['leak_burden'].max()
allL['safety'] = np.nan
msimple = allL['type'].isin(['single', 'OR', 'AND'])
allL.loc[msimple, 'safety'] = 1 - norm(allL.loc[msimple, 'leak_burden'])

def base_or_safety(logic):
    inner = logic.split(') AND ')[0].strip('(')
    genes = [g.strip() for g in inner.split(' OR ')]
    return 1 - min(leak_burden(or_leak(genes)) / (maxlb + 1e-9), 1)

mgate = allL['type'] == 'OR_AND_protease'
allL.loc[mgate, 'safety'] = [
    min(base_or_safety(l) + 0.4 * rf, 1.0)
    for l, rf in zip(allL.loc[mgate, 'logic'], allL.loc[mgate, 'rescue_frac'])
]

allL['complexity_penalty'] = (allL['inputs'] - 1) * 0.05
allL['final_score'] = (
    0.55 * norm(allL['coverage']) + 0.45 * allL['safety'].fillna(0.5) - allL['complexity_penalty']
).round(4)
allL = allL.sort_values('final_score', ascending=False)
allL.to_csv('logic_gate_recommendations.csv', index=False)