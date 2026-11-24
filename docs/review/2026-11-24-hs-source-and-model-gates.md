# HS GeoMx evidence-gate review — 2026-11-24

## Decision

**PASS: source S2 identity and the archived producer-run technical model gate.**
**HOLD: complete 156-HS-ROI paper replication, historical 24 ROI expression counts, clinical biomarkers, and inferential claims.**

This is an independent *review of actual GitHub Actions artifact bytes*, not a claim that static CI reproduces a public scientific study. Every link below points to its originating public repository run/artifact. The materials concern published/public deidentified metadata and exploratory model outputs only.

## Evidence 1 — the official workbook IS available (corrected from earlier review)

- **Official source:** Europe PMC `PMC13552446/supplementaryFiles`, retrieved in [run 45375886636](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886636).
- **Archived artifact:** [14910769107](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886636/artifacts/14910769107) (`hs-spatial-publication-supplement-linkage`), expiry **2026-12-24 10:34 UTC**.
- **Exact XLSX member:** `NIHMS2203138-supplement-Supplementary_Table_1-3.xlsx`, **48,726 bytes**, SHA-256 `fb65bdc9fd578eb4746c428faf9c7d4d4c37d6ef248205fcf36ae0df5317de9a`.
- The workbook contains **S1, S2, S3**. S2 contains **146 explicitly same-row GSM-to-patient matches** (132 HS ROI and 14 acne/control ROI), plus **24 historical HS ROI records without GEO GSM**. The checked-in crosswalk is `author_s2_geo_crosswalk.csv`.
- The independent audit `HS_SOURCE_STAGE_SITE_QC.json` reports `source_GSM_found=146`, `source_same_row_patient_explicit=146`, `patient_source_verdict=EXPLICIT_SAME_ROW`. This is **PASS for metadata source linkage**, not expression quantification or statistical significance.
- The `SUPPLEMENT_DISCOVERY.json` file still has `decision=HOLD_NO_VERIFIED_DONOR_MAPPING` as its *initial discovery-stage gate*. It is **superseded for S2 crosswalk identity only** by the separately executed verified `HS_SOURCE_STAGE_SITE_QC.json` gate in the same artifact. Do not use the discovery-stage label as the final integrated outcome.
- Code now fails closed if any S2 GSM is absent, repeated, on a wrong worksheet, or not explicitly paired with the expected patient on that same row. Offline unit tests cover these cases.

## Evidence 2 — independent sample counts and unresolved source gap

| Requirement | Verified evidence | Gate |
| --- | --- | --- |
| Original S2 GEO identities | 146 unique GSM+patient pairs in the official workbook | **PASS_METADATA** |
| Public 2026 HS source counts | 132 HS ROI, 15 donors in GEO GSE294009, within 146 GSM/ROI entries | **PASS_SOURCE_SCOPE** |
| Historical ROIs documented by 2026 paper S2 | HS41: 6 ROI; HS42: 6 ROI; HS43: 12 ROI (3 donors, total 24 ROI) | **PASS_METADATA_ONLY** |
| Historical 24 expression values | [2023 BJD source audit 45375886768](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886768), [artifact 14910169962](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886768/artifacts/14910169962) found **0 verified 2023 paper-linked ROI quantitative files** | **HOLD_NO_ORIGINAL_QUANT** |
| Full paper scope | 132 public HS GEO ROI + 24 historical metadata ROI = 156 paper HS ROI; quantitative evidence available only for the 132 subset | **HOLD_156_ROI_REPLICATION** |
| S1 versus S2 HS12 clinical stage | Four HS12 GSM: `GSM8896388`, `GSM8896393`, `GSM8896394`, `GSM8896410` disagree between S1 and S2 staging; exact basis needs author clarification | **HOLD_STAGE_DISCORDANCE** |
| Donor-aware primary subset | 82 Yonsei FFPE ROI excluding HS12: Early 57 ROI / 7 donors and Late 25 ROI / 3 donors | **PASS_COHORT_DESCRIPTION_ONLY** |

The historical 2023 DOI is `10.1093/bjd/ljad274`. The run returned HTTP 403/429/503 for some endpoints, while other metadata endpoints succeeded. "No verified files found" **does not mean the raw files definitively do not exist**.

## Evidence 3 — official methods PDF and transient retrieval

- A previous official-supplement retrieval **did succeed**: [run 44908096673](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44908096673), [artifact 14708172624](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44908096673/artifacts/14708172624), expiry **2026-12-21 13:18 UTC**.
- It indexed `NIHMS2203138-supplement-1.pdf`, **8 pages**, **16,447,790 bytes**, member SHA-256 `02546a7a94e71e6c1a7551d32f9a97a1e9adbc13ef2d9a38825810e29b8e113d`.
- [Later run 45375886638](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886638) failed with an official Europe PMC **HTTP 503**, not a proven invalid source. Subsequent code implements bounded transient retries and offline tests; a post-change official-PDF download has **not** been demonstrated by that lightweight CI alone.
- The ZIP archive's SHA-256 differed between retrievals despite equal observed length. Treat **member content SHA-256 (the PDF and XLSX)** as the durable identity, not the archive ZIP checksum alone.
- The 8-page PDF extraction established source availability but **did not** establish exact original normalization, LOQ settings, batch correction, DE model or FDR equivalence. Preserve the methods-equivalence HOLD.

## Evidence 4 — latest producer-side 44-fit QC is now tied to its source bytes

- **Producer workflow:** [run 45375886224](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886224) at commit `1a520e0dd952efeff26070788d0247e0b17c5502`, **completed successfully**. Its in-job independent verifier step also succeeded.
- **Actual producer JSON:** [artifact 14911168655](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886224/artifacts/14911168655) containing `HS_FULL_REFIT_SLIDE_DONOR_LOO_AUDIT.json`, SHA-256 `23924f48b44b16dc9a91247285094bcb1d1e847cbcbcfb3ad67c08357de1efd6`.
- **Independent QC JSON:** [artifact 14910914934](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886224/artifacts/14910914934) containing `HS_FINAL_QC_GATE.json`. Its `verified_input_sha256` is **identical** to the actual producer JSON SHA-256 above.
- QC reports 4/4 full model fits, 40/40 donor-LOO nuisance-refitted models, 10 unique donors, **0 direction reversals**; source matrix SHA-256 `943dff9ffcd6a9b75b0356a1c388863b48ec95804cf1d004ca17571a4ce5c169`. Both artifacts expire **2027-01-08 10:44 UTC**.
- The code was independently verified *inside that producer run*; this is stronger than the earlier final verifier's incorrectly pinned historical artifact. The producer run predates this documentation/test-only branch head, so do not mislabel it as a full re-run at current PR HEAD.
- The model is for the **82-ROI / 10-donor Yonsei FFPE subset**, not all 156 paper ROI. It is **exploratory direction robustness**, not a differential-expression significance result, independent prediction test, or validated biomarker.

## Practical release checklist

- [x] Official supplemental XLSX actually retrieved with a preserved member checksum.
- [x] 146 original S2 GSM/patient identities verified in the same source row.
- [x] 24 historical ROI source metadata identified, with donor IDs and ROI labels.
- [x] Producer 44-fit exploratory artifact SHA matches its producer-side independent QC gate.
- [x] Independent tests require unique exact official S2 source mapping.
- [ ] Locate **original 2023-study expression files for HS41/HS42/HS43**, with original accession(s), file hashes, ROI labels and full source provenance, or formally record an unrecoverable HOLD.
- [ ] Resolve or explicitly stratify/exclude the four HS12 stage-discordant records in any future model.
- [ ] Review source-paper PDF methods against LOQ, normalization, batch adjustment, statistical model and FDR methodology.
- [ ] Preserve member content and provenance in a durable **appropriately licensed** data archive before artifact expiration. Do not commit copyrighted paper PDFs or unlicensed raw matrices to the public repository.
- [ ] Re-run the 156-ROI author-equivalence gate **only if** the 24 quantitative ROI are validly recovered and the study design is reproducible.
- [ ] Require external/donor-independent validation and an appropriately preregistered analysis before any biomarker claim.

## Interpretation rules

- **PASS_METADATA** does not imply PASS_QUANTIFICATION.
- **PASS_EXPLORATORY_MODEL_QC** does not imply PASS_INFERENTIAL_STATISTICS.
- **HTTP_ERROR** is not evidence that a source does not exist.
- **HOLD_SOURCE** cannot be overridden with approximate or fabricated counts.
- Keep #25 in **Draft**, do not merge to main or market results as confirmed until the appropriate independent gates are satisfied.
