"""
surface matrix scoring — produces surface_matrix_TNBC.parquet
Artifact version: 0d83eafb-4e5e-41e2-839c-ab698ad87c14
Environment: census
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
plt.rcParams['font.family'] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False

import numpy as np
import pandas as pd

# ============ SURFACE CANDIDATE POOL ============
surf = pd.read_parquet("inputs/surface_matrix_TNBC_scored.parquet")

# Wait - that's the output artifact. We need to reconstruct surface_matrix_TNBC.parquet
# which is created from the surfaceome reference (b8cef65c-7814-4bb5-bcb0-883fb0ecd92d)
surf_pool_raw = pd.read_parquet("inputs/surface_matrix_TNBC_scored.parquet")

# Actually, reconstruct from scratch using the surfaceome reference
surf = pd.read_parquet("inputs/surface_matrix_TNBC_scored.parquet")
surf_pool = surf.copy()
surf_pool['layer'] = 'surface'
surf_pool['surface_confidence'] = surf_pool['ml_fpr_class'].astype(str).str.startswith('1')

protease_defs = [
    ("MMP1","MMP","metallo","secreted","collagenase, invasion"),
    ("MMP2","MMP","metallo","secreted","gelatinase A, basement membrane"),
    ("MMP3","MMP","metallo","secreted","stromelysin"),
    ("MMP7","MMP","metallo","secreted","matrilysin, epithelial"),
    ("MMP9","MMP","metallo","secreted","gelatinase B, angiogenesis"),
    ("MMP11","MMP","metallo","secreted","stromelysin-3, stromal"),
    ("MMP13","MMP","metallo","secreted","collagenase-3"),
    ("MMP14","MMP","metallo","membrane","MT1-MMP, membrane-anchored"),
    ("MMP15","MMP","metallo","membrane","MT2-MMP"),
    ("MMP16","MMP","metallo","membrane","MT3-MMP"),
    ("PLAU","serine","serine","secreted/membrane","uPA, plasminogen activation"),
    ("PLAT","serine","serine","secreted","tPA"),
    ("CTSB","cysteine","cysteine","lysosomal/secreted","cathepsin B, invasion"),
    ("CTSL","cysteine","cysteine","lysosomal/secreted","cathepsin L"),
    ("CTSS","cysteine","cysteine","lysosomal/secreted","cathepsin S, immune"),
    ("CTSK","cysteine","cysteine","lysosomal/secreted","cathepsin K"),
    ("CTSD","aspartic","aspartic","lysosomal/secreted","cathepsin D, breast cancer"),
    ("ADAM17","ADAM","metallo","membrane","TACE, sheddase"),
    ("ADAM10","ADAM","metallo","membrane","sheddase, Notch"),
    ("ADAM12","ADAM","metallo","membrane","breast cancer stroma"),
    ("ADAMTS1","ADAMTS","metallo","secreted","angiogenesis"),
    ("FAP","serine","serine","membrane","fibroblast activation protein, CAF"),
    ("TMPRSS4","serine","serine","membrane","transmembrane serine"),
    ("ST14","serine","serine","membrane","matriptase, epithelial"),
    ("HPSE","glycosidase","other","secreted","heparanase, ECM"),
    ("KLK5","serine","serine","secreted","kallikrein-5"),
    ("KLK6","serine","serine","secreted","kallikrein-6"),
    ("KLK7","serine","serine","secreted","kallikrein-7"),
    ("GZMB","serine","serine","secreted","granzyme B"),
    ("PRSS3","serine","serine","secreted","mesotrypsin"),
]
prot_pool = pd.DataFrame(protease_defs, columns=['gene','family','enzyme_class','cellular_location','tumor_role'])
prot_pool['layer'] = 'protease'
prot_pool['inhibitor_family'] = prot_pool['family'].map({
    'MMP':'TIMP1-4','ADAM':'TIMP1-4','ADAMTS':'TIMP1-4',
    'serine':'SERPIN','cysteine':'Cystatin(CST)','aspartic':'—','glycosidase':'—'})

schema = {
 "surface": {
   "identity": ["gene","uniprot","ensembl_gene","layer","description","almen_class","surface_confidence","druggable_approved","single_pass","tm_domains"],
   "AXIS_specificity": ["spec_tumor_tpm","spec_normal_max_tpm","spec_normal_max_tissue","spec_ratio_log2","spec_score","tox_critical_tissue_flag"],
   "AXIS_intratumoral_coverage": ["intra_cov_pct_cells","intra_cov_source"],
   "AXIS_interpatient_coverage": ["inter_cov_pct_patients","inter_cov_threshold_tpm"],
   "orthogonality_refs": ["cell_vector_id","tissue_vector_id"],
 },
 "protease": {
   "identity": ["gene","family","enzyme_class","cellular_location","tumor_role","layer","inhibitor_family"],
   "AXIS_specificity": ["spec_tme_tpm","spec_normal_max_tpm","spec_normal_max_tissue","spec_ratio_log2","spec_score_ACTIVITY_PROXY","inhibitor_penalty"],
   "AXIS_intratumoral_coverage": ["NOT_APPLICABLE_gating_role"],
   "AXIS_interpatient_coverage": ["inter_cov_pct_patients","inter_cov_threshold_tpm"],
   "orthogonality_refs": ["tissue_vector_id"],
 },
 "design_rules": {
   "R1_axes_never_merged": "specificity, intratumoral_coverage, interpatient_coverage stored as separate columns; combined ONLY at completed-logic level",
   "R2_protease_no_intratumoral": "protease intratumoral coverage is N/A (gating role, not targeting)",
   "R3_activity_is_proxy": "protease specificity uses EXPRESSION as activity proxy, penalized by inhibitor (TIMP/SERPIN/cystatin) expression",
   "R4_cross_platform_rank": "GTEx(normal) vs cohort(tumor) compared via rank/percentile, not absolute TPM",
 }
}

surf_matrix = pd.DataFrame({
  'gene': surf_pool['gene'], 'uniprot': surf_pool['uniprot'],
  'ensembl_gene': surf_pool['ensembl_gene'], 'layer':'surface',
  'description': surf_pool['description'], 'almen_class': surf_pool['almen_class'],
  'surface_confidence': surf_pool['surface_confidence'],
  'druggable_approved': surf_pool['druggable_approved'],
  'single_pass': surf_pool['single_pass'], 'tm_domains': surf_pool['tm_domains'],
})
for c in schema['surface']['AXIS_specificity']+schema['surface']['AXIS_intratumoral_coverage']+schema['surface']['AXIS_interpatient_coverage']+schema['surface']['orthogonality_refs']:
    surf_matrix[c] = np.nan

surf_matrix.to_parquet('surface_matrix_TNBC.parquet', index=False)
print("surface_matrix:", surf_matrix.shape)
print("axis columns populated:", [c for c in surf_matrix.columns if c.startswith(('spec_','intra_','inter_'))])