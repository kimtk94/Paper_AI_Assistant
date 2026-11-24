# HS response / MuSK-MG / NMOSD source and cohort review — 2026-11-24

**Decisions: HS previous exploratory paired analysis available; HOLD clinical response prediction and drug transfer. MG HOLD MuSK-specific cohort. NMOSD antigen-enriched cell-level mechanism only; HOLD donor-level or relapse prediction.**

## Audited archive bytes

Five original artifact ZIP SHA-256 values were recalculated and matched the GitHub API metadata. This review inspected outputs; the old analysis inputs were not freshly downloaded, and no new DESeq2 fit was performed.

| Run | Artifact | SHA-256 | Expiration UTC |
|---|---|---|---|
| [45375886569](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886569) HS failed mapping | 14912065536 | 9cc945ebdba6c5c47ed0c68a641115a02fdecf7375ab3713bd81d2b8dd33cc50 | 2027-01-08 10:37:09 |
| Same run NMOSD | 14911865416 | 13dfef6a8b21031d96bcb521986596ae96e9d143b3f15567f4e8006cd90b63e9 | 2027-01-08 10:34:13 |
| Same run MG | 14911765307 | 8ee775fe401e37da9122a51f49b06eb1ba68c9c2dbe2b9c2c63323fb2406c891 | 2027-01-08 10:34:16 |
| [45375886987](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886987) HS response | 14911391744 | 42c33d4dd895a68a70cd7c249ccb909c81297918f8e910cccaa3730c1d095336 | 2027-01-08 10:37:31 |
| [45375886069](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886069) HS covariates | 14910204503 | b34d47449cae1853f5eeb2cfaa98580735e351a245381ddaada2eaa7dc25f9e6 | 2027-01-08 10:34:06 |

## HS: anti-TNF, not anti-IL17A

Official [GSE213761](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213761) has 59 tissue samples, baseline/week-12 adalimumab, PMID 36689500. A failed run lacked series metadata and archive access (HTTP 403) despite obtaining the combined matrix. The successful response workflow used that same combined matrix SHA-256 65d8ab8b95946d556d323267ca2a1542c9e7240a8f4a9fa10a7606adfe5049d6 and raw archive SHA-256 43abe109dd41fe2ca5902906eeb58d6b334a8911aec4a335a34c311156f73bce. Prior mapping reported 59/59 and 14 lesional paired patients, but the old implementation allowed filename/partial signature matches. It is an archived result, not a newly certified full-vector mapping.

- 14 independent patients, 7 GEO responders and 7 nonresponders, 28 lesional biopsies.
- Response x time model: 20,970 tested rows, zero genome-wide BH FDR <0.05.
- Baseline YES vs NO: 20,611 tested rows, zero genome-wide BH FDR <0.05.
- 12 observed candidate genes in 3,432 exact label partitions: no change association survives targeted BH. Baseline CCL20 difference +1.50073 log2 normalized expression, exact p=0.0017483, targeted BH=0.020979; corresponding DESeq2 baseline genome-wide padj=0.985735. Retain only as a hypothesis. The normalization uses the same retrospective cohort, not a validated predictor.
- Clinical covariate output documents age, BMI, severity and smoking/gender summaries in a small cohort. Unadjusted balance p-values do not demonstrate equivalence or randomization. No adjudicated responder endpoint or independent cohort is verified by this audit.
- Require full individual count-vector verification before future fits, explicit patient/tissue/visit uniqueness, 14 lesions and consistent YES/NO 7/7. Do not interpret adalimumab as secukinumab response. New source guard remains **live-unverified** until a fresh research workflow runs.

## MG: serotype and biological unit

[GSE227835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227835), PMID 38711503: archived 40 explicit PBMC sample records, 10 AChR-positive, 10 HC, 10 seronegative pre, 10 seronegative post; 10 matching pre/post title IDs. This is not 40 independent subjects or a MuSK-positive cohort. CD22 is a reported general/seronegative MG mechanistic candidate, not validated here as a MuSK predictor.

[GSE297547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297547), BioProject PRJNA1265170: official series has 25 registrations, 19 PBMC titles (C01–C10, MG01–MG09) plus 6 V(D)J-related registrations, enriched T/B cells, AbSeq and some multiplexed immune repertoire data. Earlier CI could not obtain the series matrix. Do not count V(D)J records as additional donors; no MuSK-positive donor-level serology, treatment and response crosswalk established. No 5.2 GB archive downloaded in this review. Need actual author serology, longitudinal clinical outcome, treatment and paired sample mapping; otherwise HOLD.

## NMOSD: mechanism context, not patient-level differential expression

[GSE242232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242232), PMID 39226901, DOI 10.1016/j.immuni.2024.08.005: archived metadata has 928 cell records, 928 unique Track labels; matrix has four annotation fields plus 924 mapped cell columns. Four missing records are in the T-AQP4 title group (262 metadata vs 258 counted). T-HD 266, T-Flu 192, T-CAND 192, Treg 16 counted. Title suffixes alone are not verified biological donors. All six targeted gene rows in prior CSV are descriptive, no normalization/donor DE. Observed detection fractions for T-AQP4 vs T-HD: FOXP3 0.4729/0.0414, CTLA4 0.6783/0.4586, PDCD1 0.3062/0.4211, LAG3 0.2171/0.2782. Marker patterns are not uniformly higher; no pooled patient inference or relapse validation. Gene IDs are Ensembl (not Entrez). Require exact Track uniqueness, explicit missing-track accounting and finite nonnegative integer target counts in future runs. Donor IDs, diagnosis/serology, treatment timing and clinical endpoints must be validated independently before any donor-level test.

## Code and operational gates

- Replace partial HS source mapping with full per-library exact count and gene-universe comparison; reject duplicate identities, incomplete 59-source coverage or incompatible labels; block downstream models before R installation unless eligibility is verified.
- Strengthen NMOSD source header and Track uniqueness, preserve the known 924/928 scope, reject ragged/invalid/duplicate target rows; retain cell-only status.
- Remove stale hard-coded GSE297547 fetch result in favor of a current metadata audit; do not infer MuSK status from generic MG labels.
- Synthetic tests run in lightweight PR CI; no external research downloads, patient data publication, outreach or merge. Previously archived numerical results are unchanged. New full-vector guard has not yet passed a live GSE213761 archive fit; it intentionally blocks unsupported source formats rather than inferring linkage.
- Archive digests do not preserve original bytes; retain authorized files outside public source control before expiry.
