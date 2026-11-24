# Historical 24 HS GeoMx ROI — source request and controlled intake

**Date:** 2026-11-24  
**Current decision:** `HOLD_NO_VERIFIED_PAPER_LINKED_2023_ROI_QUANT`  
**Scope:** Historic HS41 (6 ROI), HS42 (6 ROI), HS43 (12 ROI) first described in *Suh-Yun Joh et al., British Journal of Dermatology* (2023), DOI [10.1093/bjd/ljad274](https://doi.org/10.1093/bjd/ljad274), PMID [37590958](https://pubmed.ncbi.nlm.nih.gov/37590958/).

## Source findings — do not conflate research cohorts

1. The 2023 study describes three Hurley I HS nodules from **three HS patients**, compared with three epidermal cysts. It is an earlier publication, not the 2026 Yonsei GEO study.
2. In the 2026 paper's official S2 ([GitHub run 45375886636, artifact 14910769107](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886636/artifacts/14910769107)), 24 historical HS ROIs are explicitly listed, with donor IDs HS41/42/43, author ROI labels, institutions and compartments. Their **GSM is blank**; expression values are absent.
3. The independent [2023 public-source scan 45375886768](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886768) checked Europe PMC, Crossref, PubMed-to-GDS links, NCBI GDS searches, the publisher and GEO references. No paper-specific raw DCC or 24-column count matrix was verified. Some endpoints returned **403 / 429 / 503**, making the search **incomplete**; this is not proof that the data do not exist.
4. The 2026 [GSE294009](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294009) DCC resources cover the separately identified public 2026 GEO cohort. A similar accession or different spatial-HS study is **not** an approved substitute for the historical 24 ROI.
5. Author/publisher source materials may have copyright, access conditions or human-study restrictions. Do not upload original matrices or full supplementary PDFs to this public GitHub repository.

## Requested items, in priority order

| Priority | Requested from corresponding author / authoritative source | Gate |
| --- | --- | --- |
| P0 | Original GeoMx **DCC** files or **gene/probe × 24 ROI raw counts** for HS41, HS42, HS43; direct archive/accession locator and study DOI linkage | File origin |
| P0 | Exact original ROI identifiers crosswalk to HS41:6, HS42:6, HS43:12 and compartment / panel identifiers | 24-to-24 identity |
| P0 | Genome/probe annotations, GeoMx assay panel and raw-vs-normalized count units | Platform consistency |
| P1 | ROI and negative-control QC, LOQ threshold, sequencing saturation, background/low-count filtering | Quality comparability |
| P1 | Original normalization, RUV/batch adjustment, statistical modeling and any transformed matrices if used | Methodological comparability |
| P1 | License, permitted reuse, ethics/consent and access controls; checksum of each file or downloadable archive | Governance and provenance |
| P2 | Any published GEO/ArrayExpress/BioStudies/Zenodo/other deposition accession that **explicitly** connects to DOI 10.1093/bjd/ljad274 | Independent discovery |

## Draft request — not sent

**Subject:** Research reproducibility inquiry: 2023 early-stage HS GeoMx ROI-level expression data (DOI 10.1093/bjd/ljad274)

Dear Corresponding Author,

I am reviewing the publicly available GeoMx spatial transcriptomics evidence associated with your 2023 *British Journal of Dermatology* study, “Spatial transcriptomic profiling reveals the pathophysiology of early-stage hidradenitis suppurativa.”

The supplementary ROI inventory in a subsequent HS study identifies 24 historical regions from three HS donors (HS41: six ROIs; HS42: six; HS43: twelve). We can verify those ROI labels and donor metadata, but have not located corresponding raw expression/DCC files in the public repositories checked.

Would you be willing to direct us to an official deposition or advise on an appropriate research-access procedure for these ROI-level raw counts, their exact ROI crosswalk, GeoMx gene/probe panel annotation, and associated QC/normalization methods? We would follow the applicable data-access, licensing, confidentiality and ethics requirements. The purpose is to assess reproducibility; no patient-level data or clinical predictions are requested.

Thank you for your consideration.

[Researcher / institution details]

**No outreach has been sent or authorized by this record.** Verify the actual corresponding author and preferred address from the original paper before contacting anyone.

## Intake protocol when a candidate file arrives

1. Obtain an authoritative, dated source URL or documented author transfer and capture **download/source checksum**, file name, scope and license **outside** the public repository.
2. Extract the source-verified `HISTORIC_24_EXPECTED_ROI_MANIFEST.csv` from [artifact 14910769107](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375886636/artifacts/14910769107), preserving SHA-256. If the file is regenerated, independently re-run the original S2 workbook SHA and GSM/donor verification first.
3. If and only if the candidate is a **raw nonnegative integer gene/probe × ROI CSV/TSV**, map its columns to exactly the 24 official `ROI_composite_key` values (e.g. `HS42:2`) using an author-provided explicit ROI mapping. **Never infer identity from column position.** Save the mapping separately with source proof.
4. Run the offline format-only gate below in a **private working area**, not GitHub CI against unlicensed data.

```bash
python projects/hs_spatial_gse294009/historic24_quant_candidate_gate.py \
  --manifest /private/HISTORIC_24_EXPECTED_ROI_MANIFEST.csv \
  --matrix /private/historical_24_raw_counts.csv.gz \
  --out /private/historical_24_gate
```

5. The gate outputs only bounded metadata/checksums in `HISTORIC_24_QUANT_GATE.json`. Even when the format passes, its only success-like status is `CANDIDATE_FORMAT_PASS__HOLD_INDEPENDENT_SOURCE_AND_METHOD_REVIEW`, with `source_authenticity_verified=false` and `combined_156_roi_analysis_allowed=false`.
6. If only normalized decimal values, unmatched DCC files or missing source information are available, retain `HOLD`; request a documented alternative validation design instead of coercing them into raw integer counts.
7. Independently verify provenance and **assay comparability** before any combined 156 ROI inference. Do not use a positive format test to infer biological validation.

## Artifact preservation deadlines

| Evidence | GitHub run / artifact | Expires (UTC) |
| --- | --- | --- |
| 2026 author S1-S3 XLSX and verified historic ROI manifest | 45375886636 / 14910769107 | 2026-12-24 10:34 |
| 2023 source endpoint audit (metadata only) | 45375886768 / 14910169962 | 2026-12-24 10:33 |
| Original paper methods PDF provenance index | 44908096673 / 14708172624 | 2026-12-21 13:18 |
| 2026 132-ROI exploratory model source JSON | 45375886224 / 14911168655 | 2027-01-08 10:44 |
| Producer-side independent QC result | 45375886224 / 14910914934 | 2027-01-08 10:44 |

Preserve small legal metadata/manifests in a durable access-controlled location, and preserve copyright-controlled or restricted-source files only in accordance with their terms. Artifact expiration is an evidence-retention risk, not proof of research source absence.

## Exit criteria

- **Metadata PASS:** exact 24 ROI/three donor identities from official S2 (already satisfied).
- **Candidate FORMAT PASS:** only if strictly checked against those 24 ROI and valid raw count units (new programmatic gate).
- **Source VERIFICATION PASS:** requires independently checked DOI/author/repository provenance, downloaded bytes, permissible use, ROI-level identity and assay protocol (NOT YET SATISFIED).
- **Analysis READY:** separate confirmatory design/method review; no pooled 156 ROI DE or biomarker claim before this stage.

This note changes neither original published results nor any scientific HOLD.
