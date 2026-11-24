# HS GSE294009: same-site patient-aware exploratory GeoMx sensitivity

Audit date: 2026-11-21. **Status: DESCRIPTIVE GO; confirmatory DE and therapy response HOLD.**

## Sources and pinned executable evidence

- Original manuscript: [Spatial transcriptomics maps early B-cell and T-cell activation in hidradenitis suppurativa](https://pmc.ncbi.nlm.nih.gov/articles/PMC13552446/), especially Fig. 4 and the Results subsection "HS inflammatory microenvironments are established in the deep dermis from the early stage of the disease".
- Public dataset: GEO GSE294009: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294009
- Explicit patient/ROI mapping: original NIHMS2203138 S2 workbook, verified 146/146 GSM same-row patient links.
- Source quantification: 204 GEO DCC -> WTA PKC RTS gene mapping, 3,876 ROI×gene targeted Q3 ratios (prior successful artifact Actions #44899730834).
- Executable CI: https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44906767437
- Machine-readable artifact: `hs-2026-132-roi-same-site-donor-sensitivity` (CSV, JSON, report).

## Cohort definition

- GSE294009 S2: 132 HS ROI from 15 unique HS donors, and 14 acne conglobata ROI from **1 donor**.
- Primary: Yonsei University FFPE only, exclude discordant HS12 (S1 Late vs S2 Early). Maximum **7 early and 3 late HS donors**; support decreases further by compartment.
- For each donor/gene/compartment: median of that donor's available ROI expression ratios, then compare early and late donor medians. Effect is `log2((Late_median+0.1)/(Early_median+0.1))`. This is **not the authors' official LOQ-normalized contrast**.
- Each leave-one-donor-out (LOO) check removes one independent patient, then recomputes the two group medians. Sign concordance does not imply significance.
- Sensitivity: include HS12 as S2 Early vs S1 Late, and allow all sites (mixes FFPE/frozen). The all-sites comparison is **not** batch corrected.

## Results (source-linked exploratory, no p values)

- **48** primary gene×compartment comparisons have both early and late donors.
- **13** comparisons were LOO direction-stable, had minimum support (early ≥3, late ≥2), and passed exploratory donor-level detectability thresholds; this is a screening heuristic, **not 13 validated hits**.
- **28/48** primary comparisons have direction concordant across alternate stage/site configurations.

### Deep dermis immune aggregates — primary vs all-site

| Gene | Primary Yonsei FFPE, HS12 excluded | All-site, HS12 excluded | Interpretation |
|---|---:|---:|---|
| CXCL13 | -0.9355 | -1.1334 | Early higher; direction stable |
| CXCL8 | -1.1926 | not evaluated in this table | Early higher in same-site LOO |
| JCHAIN | +1.8437 | +0.5463 | Late higher; substantial magnitude attenuation |
| MS4A1 | +0.6630 | +0.7004 | Late higher; direction stable |
| IL17F | +0.1952 | +0.1638 | Late slightly higher; direction stable, small effect |
| TNF | +0.3081 | +0.4406 | Late slightly higher; direction stable |
| CD79A | +0.2035 | -0.4884 | **Flips direction; unreliable as stage contrast** |
| CCL20 | -0.2224 | -0.2112 | Early higher; magnitude small |

These are only exploratory Q3 ratio differences. There is no formal hypothesis test, valid DEG FDR, or independent validation here.

### Figure 4 pre-specified superficial dermis direction check

The manuscript reports early-higher **CXCL13, CD79A** and late-higher **CCL20, IL17F, TNF** in superficial dermal inflammatory aggregates. In the primary within-site FFPE analysis:

| Gene | Paper direction | Primary log2 Late/Early | Concordance |
|---|---|---:|---|
| CXCL13 | Early higher | +0.3881 | **Discordant** |
| CD79A | Early higher | -1.1058 | Concordant |
| CCL20 | Late higher | +2.5850 | Concordant |
| IL17F | Late higher | +3.2224 | Concordant |
| TNF | Late higher | +3.0875 | Concordant |

**4 of 5** selected directional comparisons agree. Note superficial dermis has only **2 early vs 3 late** same-site donors. The CXCL13 sign becomes **-0.6674** in the all-site analysis, demonstrating sensitivity to cohort/site composition. A directional match cannot substitute for official author normalization, QC exclusions, multiple-testing adjustment, or precise patient-level replication.

## Final decision

- **GO**: public-data hypothesis generation from patient-mapped 132 HS ROI; B-cell/Type17/neutrophil candidate nomination with explicit stratum and missingness.
- **HOLD**: full 156-ROI quantitative replication (2023 historical 24 ROI DCC/matrix still unavailable); confirmatory differential expression or p/FDR; causal claims, therapy response prediction, and proof of new biomarkers.
- **Next scientific priority**: inspect raw probe detection, Q3/LOQ robustness and spatial compartment-specific donor coverage; independent patient or protein-level validation before an external collaborative claim.
- Controls: pre-existing main branch unchanged; no external outreach, SFDC/CRM write, or auto-merge.

The raw unadjusted Q3 ratios, `0.1` pseudocount, between-site handling, and donor-level sample sizes must always travel with the numeric results.
