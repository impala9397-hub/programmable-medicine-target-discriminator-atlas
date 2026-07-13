# Analysis code

These are the **actual scripts that produced each data artifact**, extracted from the execution
lineage of the pipeline. The analysis was run interactively (cell by cell, across sessions), so
each file is the recorded code for one output rather than a single top-to-bottom pipeline. Read
them in numbered order; each header names the artifact version it produced and the environment it
ran in.

| File | Produces | What it does |
|---|---|---|
| `01_surface_matrix_scoring.py` | `data/surface_matrix_TNBC.*` | 3-axis scoring of 2,886 surface proteins — `spec_score` (tumor vs normal-tissue, criticality-weighted `tox_load`), interpatient coverage (% patients >100 TPM), intratumoral coverage (% malignant cells, single-cell), and single-cell tumor-cell specificity. |
| `02_protease_matrix_scoring.py` | `data/protease_matrix_TNBC.csv` | Protease activity-proxy scoring — expression as an activity proxy, penalized by same-family inhibitor expression (`inhibitor_penalty`); interpatient coverage at 250 TPM. |
| `03_tissue_vector_surface.py` | `tissue_vec_surface` | Per-surface-gene normal-tissue leak vector (GTEx). |
| `04_tissue_vector_protease.py` | `tissue_vec_protease` | Per-protease normal-tissue gate vector (GTEx). |
| `05_singlecell_coexpr_both.py` | co-expression (AND) matrix | % of malignant cells expressing **both** genes of a surface pair (single-cell). |
| `06_singlecell_coexpr_either.py` | co-expression (OR) matrix | % of malignant cells expressing **either** gene (OR coverage). |
| `07_gene_organ_matrix_TOIL.py` | `data/gene_organ_matrix.csv` | Gene × organ TPM matrix on the UCSC Xena TOIL scale (GTEx + TCGA reprocessed together), for the tissue-orthogonality heatmap. |
| `08_logic_gate_synthesis.py` | `data/logic_gate_recommendations.csv` | Enumerates and scores 2–3-input logic gates `(A OR B) AND protease`; drops housekeeping proteases; computes coverage, safety, and `final_score`; ranks the recommendations. |

## Reproducing

The scripts read the input matrices from the analysis workspace and write the artifacts listed
above. Data sources: SURFY surfaceome, cBioPortal TCGA-BRCA PanCancer Atlas (Basal-like, n=171),
GTEx, CELLxGENE Census (TNBC single-cell), and UCSC Xena TOIL. See `../docs/methods.md` for the
full formulas, thresholds, and rationale behind every score.

> **Note.** Input paths have been normalized to `inputs/<filename>` (the upstream matrices each
> script consumes). The scripts are provided for transparency and method review, not as a turnkey
> `python 01…py` run — the intermediate matrices are the artifacts in `../data/` and the earlier
> numbered steps.
