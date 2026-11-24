# HS GeoMx JCHAIN: verified standR RUV4+voom donor leave-one-out sensitivity

**Checked 2026-11-22 (KST). Outcome: CONDITIONAL LOO 20/20 models PASS, all JCHAIN directions positive. Scientific conclusion remains exploratory stage-associated signal, not independently validated stage-specific HS biomarker.**

## Grounded execution and source

- [Full two-normalization standR pilot CI #44913688795](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44913688795): **SUCCESS**, original baseline fits.
- [Conditional leave-one-donor-out CI #44920062325](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44920062325): **SUCCESS**, actual finished run with auditable 10 donor-withheld refits per normalization.
- LOO run artifact: `hs-132-stanDR-RUV4-k2-donor-blocked-voom-pilot`, ID **14713199465**. Contains `HS_STANDR_RUV4_VOOM_FEASIBILITY.json` with 10 donor-specific results, excluded donor, sample support, estimated correlation and logFC.
- Source: [GSE294009](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294009), official 2026 article S2 patient-ROI crosswalk; [full WTA source matrix action #44909141520](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909141520): 18,676 genes × 132 S2-linked HS ROIs, 15 donors, checksum-pinned.
- Primary Yonsei FFPE, exclude discordant HS12: 82 ROI, 10 patients (Early 7, Late 3). Deep inflammatory dermis: 24 ROIs, Early 15/7 donors, Late 9/3 donors.
- R 4.5.2, standR 1.14.0, `standR::geomxNorm` on TMM and upperquartile separately, `geomxBatchCorrection(method="RUV4", k=2)`; downstream design `~0 + Stage + ruv_W1 + ruv_W2`; `limma::voom` and `duplicateCorrelation(block=Donor)`, robust eBayes.

## Direct observations from successful execution logs

| Normalization | All-donor JCHAIN exploratory log2 Late/Early | Baseline intra-patient rho | Filter-passing genes | Conditional LOO fits | Conditional LOO all positive | Conditional LOO logFC range |
|---|---:|---:|---:|---:|---|---:|
| TMM | **+2.19444** | 0.12359 | 2,360 | **10/10** | Yes | **+0.98545 to +2.73647** |
| Upper-quartile | **+2.04522** | 0.21681 | 2,937 | **10/10** | Yes | **+0.90985 to +2.86661** |

Each holdout removes one independent donor from the downstream count matrix and re-estimates feature filtering, voom weights, patient correlation and `limma::lmFit` StageLate−StageEarly coefficient. A Late donor holdout leaves only **2 independent Late patients**; therefore direction stability does not imply reliable magnitude, precision or statistical significance.

## Major qualifications

1. **Conditional nuisance estimation:** `geomxNorm`, NCG selection and RUV4 nuisance factors are fitted to the original **82-ROI primary set and held fixed** through the LOO fitting. The held-out donor contributed to the upstream nuisance estimates. This is **not** a fully held-out standR pipeline nor a model-independent validation.
2. **Stage/site and cell composition:** The primary model checks a same-site FFPE subset, but **does not include the slide fixed-effect dummies** in the final Stage+RUV design. The design rank audit shows a stage coefficient is estimable after slide adjustment, not that slide confounding has been eliminated by RUV4. Spatial compartment membership is a morphology-based definition, not single-cell compositional adjustment.
3. **Small independent patient count:** only seven Early and three Late donors; conditional LOO cannot replace independent cohorts.
4. **Different from original author pipeline:** absent source 2023 historical 24 ROIs; the original 156-ROI study and author normalization method argument, 200 NCG choices and complete inferential design are not replicated. No reliable paper-equivalent FDR or ROI inferential claim was produced.
5. **Novelty limitation:** HS lesion-associated JCHAIN is already shown in the literature, including qPCR and Western blot. This experiment tests only a restricted stage/compartment hypothesis.

## Research decision

- **PASS:** raw 132-ROI source provenance and two normalization model fits.
- **PASS:** 10/10 within-model donor exclusions in each method; **direction positive in every fit**.
- **GO for further testing only:** prioritize JCHAIN for exact slide-adjusted, nuisance-reestimated donor LOO and independent spatial/protein verification.
- **HOLD:** stage-specific biomarker validation, causal trajectory, therapy prediction, original paper reproduction and author-equivalent DE statistics.

Full log and artifact are preserved in GitHub Actions. No new patient identifiers created, main merge, external outreach or production write.
