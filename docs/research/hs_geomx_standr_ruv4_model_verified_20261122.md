# HS GSE294009: verified standR RUV4 k=2 donor-blocked voom pilot

**Status (2026-11-22 KST): source-validated R model feasibility PASS; original 156-ROI study reproduction and stage biomarker inference HOLD.**

## Execution evidence

- [Real CI run #44913688795](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44913688795) — **completed SUCCESS** at 2026-11-21 14:29:56 UTC.
- Artifact `hs-132-stanDR-RUV4-k2-donor-blocked-voom-pilot` (GitHub artifact ID **14710962934**) includes full machine-readable R model diagnostics.
- Original prerequisite full WTA matrix run [#44909141520](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909141520): 18,676 endogenous WTA genes × 132 HS GEO ROI, 15 HS patients. Raw matrix SHA-256 `943dff9ffcd6a9b75b0356a1c388863b48ec95804cf1d004ca17571a4ce5c169`.
- Donor and slide rank preflight [#44909306514](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44909306514): primary 82 Yonsei FFPE ROI, 10 HS patients after excluding discordant HS12 (Early 7, Late 3).
- **Deep dermal inflammatory aggregates: 24 ROIs from 10 patients, Early 15 ROI/7 donors, Late 9 ROI/3 donors**.
- Ubuntu 24.04, R **4.5.2**, standR **1.14.0**. Missing native `curl/curl.h` / ImageMagick C++ library fixed in workflow; earlier failed run #44909510661 never fitted a model.

## Actual CI-log results

| Model | RUV4 k | Deep dermis rows | Patient correlation (duplicateCorrelation) | Expressed genes after filtering | JCHAIN exploratory model log2 Late/Early |
|---|---:|---:|---:|---:|---:|
| `standR::geomxNorm(method="TMM")` | 2 | 24 ROI | 0.12359 | 2,360 | **+2.19444** |
| `standR::geomxNorm(method="upperquartile")` | 2 | 24 ROI | 0.21681 | 2,937 | **+2.04522** |

Both models passed the following actual controls: NCG selection; RUV4 factor output (`ruv_W1`, `ruv_W2`); `~0 + Stage + ruv_W1 + ruv_W2` rank; patient `duplicateCorrelation`; `voom` precision weights and `limma::lmFit` with `block=Donor`, followed by robust `eBayes`.

**Key caveats:** The RUV4/NCG fitting was performed on a cohort restricted to 82 Yonsei FFPE ROI. The manuscript instead integrated 156 HS ROI (including currently unavailable historical 24), and the original normalization method argument and precise NCG/QC choices remain unspecified. These logFCs cannot be described as the paper's original effects, nor do they establish p-value/FDR significance, effect precision, within-donor causal changes, or independent biomarker validation.

The gene's association with HS lesions has independent literature precedent. The untested novel claim is **within-compartment early-to-late trajectory**, requiring independent patient, QC, cell-composition, and ideally orthogonal assay validation. Site/slide effect assumptions, only three Late patients, and original 24 ROI source access remain serious limitations.

## Scientific interpretation

- **PASS:** complete public 132-ROI source matrix + 82-ROI same-site sample design + executable two-normalization RUV4/voom patient-blocked pilot.
- **GO:** independent-donor leave-one-out robustness and nuisance-factor leakage sensitivity testing.
- **HOLD:** inferential stage-DE, manuscript-equivalent normalization, 156-ROI replication, clinical performance/therapy prediction and novel JCHAIN biomarker claims.

No main merge. Original code committed to PR #25.
