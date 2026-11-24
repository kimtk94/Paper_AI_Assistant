# HS GeoMx GSE294009: complete source WTA reconstruction and standR feasibility

Date: 2026-11-21. **Full raw WTA source count reconstruction PASS; independent donor/model design preflight PASS; author-equivalent 156 ROI inferential reproduction remains HOLD.**

## Source evidence

- Official GSE294009 GEO 204 raw DCC files and 18,677-target GeoMx WTA PKC.
- Published S2 resolves **132 HS GEO ROI from 15 independent HS patients**, plus **14 acne conglobata GEO ROI from just 1 donor**; 24 HS historical 2023 ROI do not have verified matching public gene expression.
- Manuscript's full 156 HS ROI count = 132 new + 24 historical; these are distinct expression sources.
- Original 8-page Methods: `standR::geomxNorm`, RUV4 via `geomxBatchCorrection(k=2)`, `limma-voom` with disease status and RUV factors plus robust `eBayes` and BH. The authors' explicit `geomxNorm(method=...)`, NCG list, precise LOQ, original 24-ROI input and full design code are not recovered. Thus the following is **an independently designed feasibility experiment**, not exact original reproduction.

## Executed original data recovery

[GitHub Actions #44909141520](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909141520): **SUCCESS**. Artifact `hs-2026-132-full-wta-source-counts`.

- Raw 204/204 GSE DCC GSM matched; source series, PKC and raw tar SHA-256 matched an earlier independently verified CI.
- **18,676 endogenous WTA gene rows × 132 GSM** integer count matrix in `HS_132_HS_GENE_WTA_COUNTS.csv.gz`.
- Unique S2 donor-ROI metadata in `HS_132_AUTHOR_S2_MANIFEST.csv` and probe-to-gene counts in a separate source audit.
- Matrix SHA-256: `943dff9ffcd6a9b75b0356a1c388863b48ec95804cf1d004ca17571a4ce5c169`.
- All 15 HS donors represented. No formal DE or normalization has been claimed.
- An arbitrary minimum of 12,000 observed RTS entries per DCC was found unjustified by the raw source and **removed**; actual DCC identifiers and full genomic feature/source checks are now the acceptance criterion.

## Independent-patient and slide-stage design audit

[GitHub Actions #44909306514](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909306514): **SUCCESS**. Artifact `hs-132-stanDR-donor-rank-preflight`.

- For the primary analysis, **Yonsei University FFPE only**, excluding S1/S2 discordant `HS12`: **82 ROI, 10 distinct HS donors, 57 Early ROI from 7 donors, 25 Late ROI from 3 donors**.
- 4 slide-label groups. The model matrix including slide dummy variables and an Early/Late contrast **has full additional stage rank**, so stage is not perfectly aliased to slide in this restricted cohort.
- In contrast, adding all donor IDs as fixed effects makes the stage coefficient **perfectly aliased**, since every patient belongs to exactly one stage. Do **not** fit `~donor + stage` expecting an estimable stage coefficient.
- Proper design requires patient-level aggregation, random/correlation modeling (for instance limma `duplicateCorrelation(block=Donor)`), or independently planned repeated-measures statistics. It cannot treat 82 ROI as 82 independent biological replicates.

| Compartment | Early ROI | Late ROI | Independent patients (E/L) | Identifiable after slide fixed-effect adjustment? |
|---|---:|---:|---:|---|
| Deep dermal inflammatory aggregates | 15 | 9 | 7 / 3 | YES, for model feasibility only |
| Epidermis | 19 | 5 | 7 / 2 | YES, but low Late n |
| Superficial dermal inflammatory aggregates | 3 | 5 | 2 / 3 | NO (slide confounding in this subset) |
| Sinus tract | 10 | 6 | 3 / 1 | NO |
| Hair follicle | 6 | 0 | 2 / 0 | NO |
| Gland | 4 | 0 | 2 / 0 | NO |

## standR feasibility pilot

An independent experimental R workflow is in the PR at `.github/workflows/hs-geomx-standr-ruv4-pilot.yml` using:

1. Verified 18,676 × 132 source matrix; isolate 82 same-site FFPE ROI.
2. `standR::readGeoMxFromDGE` on raw endogenous WTA counts, then `geomxNorm(method='TMM')` and `geomxNorm(method='upperquartile')` in separate branches. **Neither is asserted as the paper's exact unspecified method.**
3. `findNCGs(batch_name='SlideName', top_n=200)`; heuristic NCGs are not the authors' original list.
4. `geomxBatchCorrection(k=2, factors='Group', method='RUV4')` with Group combining S2 stage and ROI compartment.
5. Primary DeepDermis ROI (15 Early vs 9 Late) with `~0+Stage+ruv_W1+ruv_W2` model design, `voom` and `duplicateCorrelation(block=Donor)` plus `limma` exploratory model fit.
6. Reject structural alias, unavailable NCGs/W covariates, ill-conditioned design, unstable within-donor correlation, and missing JCHAIN.
7. Even if all steps run, **do not publish nominal or FDR findings as manuscript replication**; assay QC and only three Late independent donors remain the limiting factor.

Workflow run: [#44909510661](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909510661).

## Scientific decision

**GO:** reproducibly reconstructed raw complete WTA count input and independent patient design feasibility for an exploratory standR pilot.

**HOLD:** exact manuscript normalization/QC, full integrated 156-ROI expression, formally validated early/late DE, JCHAIN stage-specific biomarker novelty or prognostic/treatment-response claims. These holds persist even if the R pilot can fit a mathematical model.

No user-facing patient identifier beyond published author S2; no outreach, CRM write, branch merge, or external changes.
