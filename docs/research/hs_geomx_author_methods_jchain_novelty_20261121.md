# HS GeoMx: original Methods audit, JCHAIN validation, and novelty gate

Updated 2026-11-21. **Scientific status: GO for exploratory hypothesis generation, HOLD for author-equivalent quantitative replication and stage-specific biomarker claims.**

## Verified original sources and GitHub Actions evidence

- Source publication: [Spatial transcriptomics maps early B-cell and T-cell activation in hidradenitis suppurativa](https://pmc.ncbi.nlm.nih.gov/articles/PMC13552446/), DOI [10.1016/j.jid.2026.01.041](https://doi.org/10.1016/j.jid.2026.01.041).
- Supplementary Methods source: `NIHMS2203138-supplement-1.pdf`, 8 original source pages from Europe PMC, SHA256 `02546a7a94e71e6c1a7551d32f9a97a1e9adbc13ef2d9a38825810e29b8e113d`.
- [Original Methods provenance CI #44908096673](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44908096673) — SUCCESS. Its bounded report and excerpts are retained, but the licensed original PDF is **not uploaded to the repository**.
- [Original DCC/PKC negative-probe LOQ/Q3/Q90 CI #44907538365](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44907538365) — SUCCESS, **132/132** author-linked HS ROI exploratory LOQ computed.
- [Methods equivalence and JCHAIN novelty CI #44908294784](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44908294784) — SUCCESS, with **fail-closed assertions** preventing misrepresentation of original methods and novelty.
- Machine-readable verdict: workflow artifact `hs-geomx-author-methods-equivalence-and-jchain-novelty`.

## What the original Methods actually report

**Source PDF p. 1**, spatial integration and normalization:
- GeoMx NGS Pipeline DCC files; `standR` in R 4.2.2.
- Normalization using **`geomxNorm`**.
- Unwanted-variation correction using **`RUV4` via `geomxBatchCorrection`**, with **`k=2`**.
- The authors state that donor and slide effects can induce correlation between ROIs; RUV4 addresses latent nuisance factors.
- PCA on quantile-normalized values after log transformation; `limma::removeBatchEffect` also appears in their variance-assessment description. The exact use of the latter relative to the primary differential-expression model should not be assumed to be a second sequential adjustment without the analysis code.

**Source PDF p. 2**, differential expression:
- `limma-voom` on `DGEList`.
- Disease status plus RUV factors included in the design matrix.
- Precision weighting with `voom`, robust `eBayes`, and BH-adjusted p-values.
- The description discusses mitigating ROI pseudoreplication; **RUV4 alone is not proof of independent donor-level hypothesis testing**. The exact statistical blocking/correlation implementation has not been independently reconstructed.

**LOQ:** no explicit LOQ equation was detected in the searchable 8-page official Methods PDF. A standard GeoMx negative-probe threshold is not evidence of the precise preprocessing executed by these authors. This means the **original LOQ filtering regime remains unspecified**, not that no LOQ or QC was used.

## Method-equivalence check

| Analysis property | Original paper | Current 132-ROI work | Status |
|---|---|---|---|
| Raw DCC/PKC provenance | GeoMx DCC | Raw 204 DCC, S2 HS 132 | Partial source match |
| Normalization | `standR::geomxNorm`; quantile expression referenced | per-ROI Q3 ratio; Q90 stress test | **NOT reproduced** |
| Batch correction | `RUV4 k=2`, `geomxBatchCorrection` | Yonsei FFPE restriction; no RUV4 | **NOT reproduced** |
| Inferential model | `limma-voom`, RUV covariates, robust `eBayes`, BH | independent-donor medians, leave-one-out sign | **NOT reproduced** |
| LOQ | Not specified in retrieved Methods PDF | Exploratory clamped negative-probe geometric threshold | **Method unknown** |
| Patient/ROI population | 156 HS ROI including 24 old-study ROI | 132 matched new-study HS ROI | **Partial cohort** |

The exploratory source-data JCHAIN deep-dermis ratio of **+1.844 log₂(Late/Early)** after per-donor medians is **not** an author-equivalent effect size and is **not** formal DEG replication.

## Independent literature: JCHAIN is not a new HS lesion-associated discovery

A separate study, [IgA autoantibodies promote inflammation, Th17 polarization and fibrotic responses in hidradenitis suppurativa](https://www.nature.com/articles/s41467-026-70883-5), DOI [10.1038/s41467-026-70883-5](https://doi.org/10.1038/s41467-026-70883-5), already reported:
- JCHAIN lesional HS skin increased relative to healthy controls by qPCR (**15 HS lesions versus 3 controls**).
- **J chain protein** elevated in HS lesional lysates by Western blot; IgA localized to HS lesions by immunofluorescence.
- Independent public bulk RNA-seq **GSE155176** supports lesional versus non-lesional JCHAIN/IGHA expression differences.
- Single-cell RNA-seq **GSE249793** localized JCHAIN/IGHA transcripts in B cells and plasma cells.

**Interpretation:** This supports the presence of a JCHAIN-related plasma-cell/IgA program in HS lesions, but does **not** validate a **deep dermis, Early-vs-Late, donor-blocked** spatial effect. HS-associated JCHAIN elevation is an established observation, so claiming its basic discovery as novel is inappropriate.

## Current GO/HOLD decisions

| Claim or use | State |
|---|---|
| Provenance-checked GEO raw ROI donor mapping (132 HS ROI, 15 HS donors) | **PASS** |
| Exploratory negative-probe-QC sensitivity and patient-median trend | **GO for hypothesis only** |
| JCHAIN disease-level HS lesion association | **Published prior evidence; not novel** |
| JCHAIN Early-vs-Late deep dermis marker | **HOLD: no method-equivalent or independent stage-specific replication** |
| 156-ROI complete quantification (24 historic missing expression) | **HOLD** |
| Exact LOQ/standR/RUV/limma author pipeline and FDR reproduction | **HOLD** |
| Treatment response or causal claim | **HOLD** |

## Next executable scientific priority

Reconstruct a **132-ROI-only** `standR::geomxNorm` / `RUV4 k=2` / `limma-voom` proof-of-feasibility in GitHub Actions, with explicit donor ID and design-matrix diagnostics. Do **not** claim the 156-ROI paper was reproduced; compare donor-level direction with the current Q3 stress test only after model/normalization choices are verified. For novel clinical relevance, the key untested endpoint is **stage-specific within-compartment** JCHAIN variation, not its existence in HS lesions. Human donor blocking remains a mandatory inferential gate.

This report changes only the research PR; there was no production/main merge, outreach, or customer communication.
