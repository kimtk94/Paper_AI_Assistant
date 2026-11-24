# IPF plasma biomarker study — proposed design and launch gates (2026-11-24)

**Internal research design only. HOLD actual enrollment and clinical inference. No new participants, specimens or clinical records were obtained. All thresholds and analytic choices below are proposals requiring investigators, laboratory and biostatistician approval before outcome-linked analysis.**

## Evidence and novelty

- [Korean IPF cohort, 2026, DOI 10.1186/s12931-026-03504-w](https://doi.org/10.1186/s12931-026-03504-w): 107 people including 19 familial IPF; longitudinal lung function and transplant/survival. It does not establish that 107 matched baseline plasma specimens can be reused. Published cohort IRB/consent does not authorize this proposed study.
- [IPF-PRO mass-spectrometry proteomics, 2025](https://doi.org/10.1186/s12931-025-03377-5): 299 IPF, 99 controls after QC. [Korean multicenter observational biomarker study, 2026](https://pubmed.ncbi.nlm.nih.gov/43065028/): 118 pirfenidone-treated patients; 12-month outcome and nine blood markers. [ISABELA biomarker analysis, 2026](https://pubmed.ncbi.nlm.nih.gov/41704719/): 1,280 IPF trial participants, progression/mortality biomarkers. Generic biomarker/prognostic novelty is not established.
- [PR #23 recheck](https://github.com/kimtk94/Paper_AI_Assistant/blob/research/ipf-pxd010965-feasibility-20261121/docs/review/2026-11-24-ipf-fcn3-gpx3-evidence.md): raw vs median-centered plasma log2FC FCN3 -0.846/-0.337, GPX3 -2.102/-1.593; consistent negative LOO direction. Source-matched lung RNA and independent cohort provide limited endothelial/stromal context, not matched protein origin or clinical prediction. The new review uses 20 cells/donor-celltype and >=5 donors/group; the prior protocol quoted different thresholds. No direction of FVC association is assumed from case-control differences.

## Cohort, index and outcome estimand

- Adults with independently adjudicated IPF, one unique subject key and one index plasma specimen per patient. Resolve duplicate/discordant records before inclusion; proposed EDTA-plasma primary cohort, not pooled serum or other anticoagulants without assay-comparability review.
- Index baseline lung function is the closest eligible, quality-reviewed native-lung FVC to index plasma within +/-30 days (choose the nearest by prespecified rule, tie earlier; investigator adjudication before protein outcomes). Exposure measured at baseline. FVC values use a consistent reference equation/units; record FVC %pred and mL. Followup chosen closest to 365.25 days from baseline FVC within **270-450 days**, tie earlier, never after death/transplant. Capture actual dates. These are day-window operational choices, not exact calendar months.
- Primary outcome: annualized absolute FVC %pred change = (followup %pred - baseline %pred) * 365.25 / observed days, units **percentage points/year**, not relative percent change. Primary estimand is association among observed eligible native-lung FVC records, with explicit survivor/observation selection limitation.
- Record death, transplant, missed visits, withdrawal/loss, and last-known observation. Do not code death/transplant/missing as stable FVC or zero change; do not impute post-death/transplant FVC. If followup endpoint is ambiguous, hold adjudication. Secondary exploratory time-to-first documented >=10 percentage-point absolute FVC decline from baseline, all-cause death, or lung transplant in 365 days, with ascertainment and censoring rules specified before analysis. No assumption that people with missing FVC had no event. If adequate serial FVC, sensitivity mixed-effects trajectory (biomarker x time); assess informative truncation/attrition with statistician, and report limitations. No automatic clinical recommendation.
- Baseline covariates: baseline FVC %pred, age, sex, smoking, antifibrotic use/timing, batch/site; capture baseline DLCO, family history, height/reference equation, eGFR/liver disease, acute infection/exacerbation, and treatment changes where available. Freeze covariate coding and missing-data rules before analysis. Avoid adding every variable to a small model.

## Assay and specimen plan

| Stage | Proposed activity | Release gate |
|---|---|---|
| Inventory | Confirm IRB/data-use approval, sample IDs, EDTA tube, volume, processing delay, storage temperature/duration, freeze-thaw, hemolysis and FVC linkage | Clinical/lab owner reviews actual inventory; no guessed specifications |
| Pilot | 6-10 independent representative baseline specimens; blinded repeated pooled QC, blanks, suitable standards; DIA detectability and target peptide specificity | Lab predefines validated SOP, peptide identity, usable signal range, precision, interference and specimen requirements |
| Target verification | PRM with appropriate internal standards or validated immunoassay; confirm matrix effects, linearity, LLOQ and batch reproducibility | Candidate-specific method qualification; provisional planning targets >=80% quantifiable pilot specimens and pooled QC CV <=20%, subject to lab review; not universal standards |
| Exploratory association | 60-80 complete records if attainable; randomize/balance plates by baseline characteristics independently of outcomes; blinded analysis, QC decisions locked before unblinding | Independent cohort and SAP review; not confirmatory clinical validity |
| Independent replication | Separate nonoverlapping temporal/site cohort using frozen assay and analysis plan | Reproducible effect, precision, incremental value and clinical validation before claims |

Do not claim Novogene-specific specimen volume, tube, shipping or platform specifications from this draft; obtain actual laboratory manual and assay-team confirmation. Capture sex, age, smoking, treatment, batch and sample QC; record missingness and exclusions by stage. Do not upload any patient-level manifest, assay result, identifiers or small-cell clinical report to public GitHub/CI. Use authorized local systems only.

## Sample size — illustrative, not formal power

For two prespecified two-sided biomarkers, use Bonferroni family alpha 0.05 (per marker 0.025) for conservative planning, 80% power and Fisher-z approximation for a simple unadjusted correlation. Proposed reporting uses Holm for the two primary p-values. These are assumed correlations, not estimates from prior cross-sectional effect sizes.

| Assumed absolute r | Complete records | Screened at 80% baseline eligibility x 80% followup |
|---:|---:|---:|
| 0.25 | 149 | 233 |
| 0.30 | 103 | 161 |
| 0.40 | 56 | 88 |

Single marker alpha .05, r=.30: 85 complete. 80 complete at 64% coverage needs 125 screened. With 107 screened at 64% coverage, expected complete count ~68.5; actual yield is unknown. These do not account for adjusted partial effect, residual SD, model degrees of freedom, measurement error, informative loss or external validation. Formal power/simulation and a feasible model must be approved from aggregate pilot inputs; an 80-case pilot is not an adequately powered confirmatory trial.

## Statistical analysis plan (proposed)

1. Clinical-only reference model: annualized FVC percentage-point change, with a prespecified parsimonious core covariate set (baseline FVC, age, sex, baseline antifibrotic use, smoking; batch/site adjustment if identifiable). Final coding/parameter count decided before seeing associations; no univariate outcome-driven covariate selection.
2. Fit separate nested candidate models with log2 FCN3 or GPX3, report beta per doubling with 95% CI, two-sided Holm adjustment for the two primary tests. Assess residual behavior/influence, robust standard errors where appropriate, nonlinearities only as labeled sensitivity. A combined model and DIA-wide FDR are secondary exploratory analyses. Baseline DLCO, renal/hepatic variables and treatment changes as prespecified sensitivity covariates if justified. Do not infer causal drug effect.
3. Inspect missingness and eligible-vs-excluded distributions in secure aggregate review. Do not replace non-detects with zero or choose imputation based on favorable p-values. Use lab- and SAP-defined missingness rules; sensitivity for censored values and informative FVC availability. Secondary clinical event endpoint requires adequate events and ascertainment; no model if sparse.
4. Estimate incremental performance versus clinical-only model using optimism-corrected internal resampling (all preprocessing inside resampling), report calibration and out-of-sample error with uncertainty as exploratory only. Independent temporal/site cohort required for external validation. No arbitrary classifier threshold, AUC promise or clinical utility claim.

## Code and authorization

Legacy feasibility.py remains a preliminary complete-record screen. Duplicate subject rows are now all excluded until adjudication, and implausible adult ages and malformed CSV schema are rejected. Its GO labels do not authorize research. The stricter study_readiness.py v2 requires explicit source/diagnostic/preanalytic/endpoint review flags and separates observed native-lung FVC from terminal events, missing and unknown outcomes. It emits aggregate counts only; all readiness outputs explicitly disallow clinical-analysis authorization and prognosis validation. The proposed EDTA restriction, 18-120 age plausibility and 270-450 day window require investigator approval.

```sh
python proposals/ipf_clinical/study_readiness.py --out /secure/planning.json
python proposals/ipf_clinical/study_readiness.py --manifest /secure/authorized_v2_manifest.csv --out /secure/aggregate_readiness.json
python -m unittest discover -s proposals/ipf_clinical -p 'test_*.py' -v
```

The local aggregate v2 output is not a substitute for independent source or IRB review. Without actual authorized data, the status is HOLD_NO_CLINICAL_MANIFEST. Before research launch require institutional approval, data-use agreement, assay-team feasibility, exact specimen availability, endpoint adjudication, SAP/power signoff and an independent validation plan. No contact, quotation, SFDC update or clinical recruitment is authorized by this document.
