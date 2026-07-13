"""
tissue vector protease — produces tissue_vec_protease.parquet
Artifact version: 082f019c-2ee7-481d-82b5-da94d344999b
Environment: python
Extracted from the artifact's execution lineage (the actual code that produced it).
"""
import pandas as pd
import json

# Load GTEx protease normal data
gp = json.load(open("inputs/gtex_protease.json"))
gpdf = pd.DataFrame(gp['rows'])
gtex_prot = gpdf.pivot_table(index='geneSymbol', columns='tissueSiteDetailId', values='median', aggfunc='first')

# Load GTEx surface normal data to get tissue columns
gtex_prior = pd.read_csv("inputs/STEP2_gtex_normal_expression.csv")
gtex_surf_full = gtex_prior.set_index('geneSymbol')
tissue_cols = list(gtex_surf_full.columns)

prot_genes = ['MMP1','MMP2','MMP3','MMP7','MMP9','MMP11','MMP13','MMP14','MMP15','MMP16',
 'PLAU','PLAT','CTSB','CTSL','CTSS','CTSK','CTSD','ADAM17','ADAM10','ADAM12','ADAMTS1',
 'FAP','TMPRSS4','ST14','HPSE','KLK5','KLK6','KLK7','GZMB','PRSS3']

prot_tissue_vec = gtex_prot.reindex(columns=tissue_cols).reindex([g for g in prot_genes if g in gtex_prot.index])

prot_tissue_vec.to_parquet('tissue_vec_protease.parquet')