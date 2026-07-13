# Methods

## 1. Candidate pools

- **Surface pool** — SURFY surfaceome, 2,886 proteins (2,096 high-confidence class-1;
  408 druggable/approved).
- **Protease pool** — 30 literature-curated tumor-associated, extracellularly active proteases
  (MMP ×10, serine ×10, cysteine ×4, ADAM ×3, ADAMTS/aspartic/glycosidase ×1 each), out of the
  ~550–600 human degradome. A pilot subset, not the full degradome.

## 2. Scoring axes (kept separate; combined only at the completed-logic level)

### Surface — three independent axes
- **Bulk-organ specificity (`spec_score`)** = `0.45 × pctrank(log1p tumor) + 0.55 × (1 − pctrank(log1p tox_load))`.
  - `tumor` = median RSEM in TCGA Basal-like (n=171).
  - `tox_load` = Σ (GTEx median TPM × tissue-criticality weight). Criticality: brain/heart/liver/lung/kidney = 3.0,
    pancreas/blood = 2.5, … breast/skin/adipose = 1.0, cultured lines = 0.5, default = 1.5.
  - Rank/percentile transform absorbs the RSEM-vs-TPM platform difference.
- **Intratumoral coverage** = fraction of malignant cells expressing the gene (single-cell detection rate;
  a dropout-limited lower bound).
- **Interpatient coverage** = % of the 171 patients with tumor expression ≥ **100 TPM**
  (density threshold reflecting CAR-recognition-relevant surface density).

### Protease — activity proxy, no intratumoral axis
Proteases gate by environment, so intratumoral coverage is **not applicable** — the enzyme only needs to be
active in the tumor microenvironment.
- **Activity-proxy specificity** = `(0.55 × rank(spec_ratio) + 0.45 × (1 − rank(tox))) × (1 − 0.25 × inhibitor_penalty)`.
- **`inhibitor_penalty`** = endogenous-inhibitor expression ÷ (protease + inhibitor) at TME median RSEM
  (e.g. TIMPs for MMPs) — net activity, not raw abundance.
- **Interpatient coverage** = % patients with TME expression above threshold.

### z-score (complementary)
TCGA "z-score vs normal" compares tumor breast tissue against **adjacent normal breast** — orthogonal to
`spec_score`, which compares TNBC against all GTEx organs.

## 3. Narrowing to the single-cell set

Single-cell co-expression is expensive, so the surface pool was first reduced by a gate
`(spec_tumor_tpm ≥ 30) AND (inter_cov ≥ 40%) AND (spec_score ≥ 0.5)` → 699, then top-45 by a rank-blended
pick score + 11 forced validation anchors (EGFR, MET, FOLH1, MSLN, MUC1, TACSTD2, ERBB2, EPCAM, CD276,
MUC16, VTCN1) = **56** genes carried into single-cell analysis.

## 4. Orthogonality

- **Cell orthogonality (single-cell)** — CELLxGENE Census (2025-11-08), real "triple-negative breast carcinoma"
  cells (161,048 cells; 48,862 malignant). Pairwise malignant-cell co-expression → OR/AND coverage.
- **Tissue orthogonality (GTEx / TOIL)** — per-organ leak burden; surface leak threshold 10 TPM,
  protease gate threshold 20 TPM.

## 5. Logic synthesis

- **Coverage** — single-cell malignant co-expression (single = detection rate, OR = either, AND = both).
- **Safety** = `1 − normalized(leak_burden)`; a protease AND-gate raises safety by up to +0.4 in proportion
  to the leak organs it silences (`min(base + 0.4 × rescue_fraction, 1.0)`).
- **Final score** = `0.55 × norm(coverage) + 0.45 × safety − 0.05 × (inputs − 1)`.
  (An earlier `escape` term was dropped: `escape ≈ 1 − coverage` was redundant with coverage.)
- **Housekeeping filter** — proteases active in >8 tissues (>10 TPM) or with Gini < 0.4 are excluded as
  non-specific; final gate proteases = **MMP11, MMP7, ADAM12, FAP, PLAU**. MMP11 is the representative gate.

## 6. Comparable-scale heatmap

GTEx median TPM and TCGA RSEM differ ~30–100× in scale and cannot be compared directly. The organ heatmap
therefore uses **UCSC Xena TOIL** (`TcgaTargetGtex_rsem_gene_tpm`), where GTEx and TCGA are reprocessed on a
single log2(TPM+0.001) pipeline. 32 GTEx sites + TCGA BRCA basal (171) shown; e.g. EPCAM TNBC 323 vs Breast 16.

## Limitations

- Two grid cells (surface protein, extracellular protease), one tumor type (TNBC) — a pilot.
- Bulk tumor expression uses PAM50 Basal-like as a TNBC proxy (TCGA has no IHC receptor status).
- Protease "activity" is an expression-based proxy corrected for inhibitor expression, not a measured assay.
- Single-cell detection rate is a dropout-limited lower bound on true expression.
- Extends to the full surfaceome/degradome and other tumor types with the same framework.
