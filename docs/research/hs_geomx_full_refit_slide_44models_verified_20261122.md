# HS GSE294009 — final original-DCC full-nuisance LOO and slide-adjusted JCHAIN audit

**Updated 2026-11-22 KST · Technical sensitivity PASS, scientific disease-stage biomarker confirmation HOLD.**

## Provenance and actually executed CI

- Original HS spatial paper: [Spatial transcriptomics maps early B-cell and T-cell activation in hidradenitis suppurativa](https://pmc.ncbi.nlm.nih.gov/articles/PMC13552446/), DOI [10.1016/j.jid.2026.01.041](https://doi.org/10.1016/j.jid.2026.01.041).
- Original GSE294009 WTA data: **204 public DCC**, raw matrix reconstructed into **18,676 endogenous genes × 132 author-S2-linked HS ROIs** (15 independent HS donors); 14 public acne conglobata ROI are not independent controls.
- Original full-paper **156 HS ROI = 132 new HS ROIs + 24 older 2023 HS ROIs**, with the latter quantitative source still unavailable. There is **no** full 156 ROI author replication.
- Matrix SHA-256: `943dff9ffcd6a9b75b0356a1c388863b48ec95804cf1d004ca17571a4ce5c169`.
- [Original raw WTA reconstruction CI #44909141520](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909141520): SUCCESS.
- [Independent donor/slide design rank preflight #44909306514](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909306514): SUCCESS.
- [Real full-nuisance/slide model fitting #44925755817](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44925755817): **SUCCESS**, source-verified 44 model fits. Artifact `hs-132-full-nuisance-LOO-and-slide-sensitivity` ID **14715499837**, `HS_FULL_REFIT_SLIDE_DONOR_LOO_AUDIT.json`.
- [Independent result-artifact re-audit #44926880011](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44926880011): **SUCCESS**, 44/44 model records, 40/40 positive donor-fold directions, 10 unique withheld donors, exact matrix provenance, and all publication-significance guardrails verified. Artifact `hs-final-44-full-nuisance-models-JCHAIN-direction-QC` ID **14715707094**.
- Original script `projects/hs_spatial_gse294009/full_refit_slide_loo.R`; independent verifier `verify_full_refit_slide_artifact.py`; GitHub workflow `.github/workflows/hs-geomx-standr-ruv4-pilot.yml` and `hs-geomx-final-donor-slide-audit.yml`.
- Technical correction: raw source slide names such as `HS slide 2 plz` contain spaces, which are invalid in `limma::makeContrasts` coefficient names. `make.names` was applied to **design matrix coefficient names only**. Source slide names, grouping and original metadata remain unchanged. The initial failed slide fits are coding faults, not negative biological results.

## Restricted source population and designs

- **Same-site Yonsei University FFPE HS samples**, excluding S1/S2 clinical-stage discordant **HS12**: **82 ROI from 10 HS donors** (57 ROI Early/7 donors; 25 ROI Late/3 donors).
- Primary matched tissue morphology **deep dermal inflammatory aggregates**: **24 ROI** (15 Early ROI from 7 patients; 9 Late ROI from 3 patients).
- Per model, `standR::geomxNorm` with `method="TMM"` or `"upperquartile"`. The authors' actual normalization method parameter remains unverified.
- Negative control genes from `standR::findNCGs` per cohort/slide (`top_n=200`), with `standR::geomxBatchCorrection(method="RUV4",k=2,factors="Group")`.
- Deep dermis regression (i) `~0 + Stage + ruv_W1 + ruv_W2`, (ii) `~0 + Stage + SlideName + ruv_W1 + ruv_W2`, with `limma::voom` and `limma::duplicateCorrelation(block=Donor)` before coefficient fitting.
- Every one of 10 leave-one-donor-out folds **removes every ROI from that patient before upstream standR normalization, NCG estimation and RUV4**, then recomputes gene filtering, voom and donor-blocked correlation. This is stricter than the preceding conditional LOO, which retained nuisance estimates from all donors.
- Stage and donor fixed effects cannot be entered together due stage-within-donor alias; the correlation method addresses repeated ROI, but does not increase independent donor count.

## Actual model results from CI and independent audit

All effects are **exploratory model estimates** of log₂(Late − Early) differential expression. No nominal p-value, FDR or confidence interval was asserted. "PASS" means model fitting and predeclared sign robustness, **not statistical significance**.

| Normalization | Covariates | Full data JCHAIN log₂FC | Successful held-out-donor refits | Full reestimated LOO log₂FC range | Sign reversal |
|---|---|---:|---:|---:|---:|
| TMM | Stage + RUV4(k=2) | **+2.19444** | **10/10** | **+1.04553 to +2.93155** | 0 |
| TMM | Stage + slide + RUV4(k=2) | **+3.65453** | **10/10** | **+1.94354 to +4.87994** | 0 |
| Upper-quartile | Stage + RUV4(k=2) | **+2.04522** | **10/10** | **+1.03131 to +2.89449** | 0 |
| Upper-quartile | Stage + slide + RUV4(k=2) | **+3.13406** | **10/10** | **+1.54564 to +4.01065** | 0 |

**44/44 model fits successful** (4 full-donor models + 40 donor-withheld fits); **40/40 LOO JCHAIN effects are positive**. These counts were directly verified by the fast independent Actions audit from the original completed R artifact, not inferred from an Actions SUCCESS badge.

## Interpretation and scientific limits

1. **Exploratory directional stability is substantial.** With independent-donor unit removal and full nuisance reestimation, JCHAIN's Late > Early direction persists across both normalization choices, with and without explicit slide effects.
2. **Magnitude is model-sensitive.** Adding slide factors shifts TMM log₂FC from +2.19 to +3.65 and upper-quartile from +2.05 to +3.13. This is **not evidence that the biologic effect became stronger**; it raises concerns about limited sample support, model conditioning, nuisance-factor collinearity, and slide-dependent estimation. A confidence interval / sensitivity on donor-level counts and model condition numbers is needed before effect size is interpreted.
3. **Small independent cohort:** only **3 Late patients**. Omitting a Late patient leaves **2 Late patients**. The number 24 ROI is not n=24 independent biological replicates.
4. **This is not predictive cross-validation.** The withheld donor's data are not used in that fold's upstream preprocessing or coefficient fit, but no held-out donor prediction, prediction loss or accuracy was evaluated. It is **leave-one-donor-out model refit sensitivity**, not external replication.
5. **Not the author-equivalent paper model.** The authors reported a 156-ROI integration with their own normalization/NCG, RUV4, limma-voom choices, and 2023 historical data; their exact settings and all expression matrices are not replicated.
6. **Uncertain clinical novelty.** HS lesion-associated JCHAIN/IgA biology was published in an independent 2026 study [DOI 10.1038/s41467-026-70883-5](https://doi.org/10.1038/s41467-026-70883-5). Confirming JCHAIN in HS skin is not a new discovery. This study explores only a spatial Early-versus-Late compartment-restricted difference, without orthogonal stage-specific independent validation.
7. **Treatment outcome and clinical prediction remain outside the demonstrated evidence.** No confirmed stage-DE FDR, therapy prediction, causality or prognostic effect was measured.

## Decision

- **PASS technical evidence:** source provenance; 18,676 × 132 count reconstruction; donor/slide design-rank preflight; RUV4 + patient-blocked modeling; 4 full fits and 40 full-nuisance donor-refits, all positive JCHAIN direction; independent JSON audit.
- **GO limited hypothesis:** JCHAIN is an **exploratory spatial stage-associated candidate**, suitable for orthogonal, independent validation design.
- **HOLD claims:** validated biomarker, stage-specific p/FDR, novel HS-associated gene, 156-ROI manuscript reconstruction, or treatment-response causal inference.
- **Next priorities:** report patient-level expression distributions and confidence-interval uncertainty with only 7 Early/3 Late donors; inspect batch covariate correlation/condition numbers; identify an independent stage-specific spatial/protein cohort; do not infer from ROI-level n or mix inconsistent stage labels.

All files are in the research PR; no main merge, customer outreach, or Salesforce changes.
