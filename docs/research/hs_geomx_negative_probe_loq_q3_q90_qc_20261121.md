# HS GeoMx GSE294009: DCC negative-control probe LOQ/Q3/Q90 QC

Date: 2026-11-21 · **CI completed successfully** · Scientific status: **descriptive only; no confirmed differential expression**.

## Source and reproducible execution

- [Primary study, JID 2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13552446/), DOI `10.1016/j.jid.2026.01.041`.
- [GEO GSE294009](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294009).
- [NanoString/GeoMxTools WTA reference workflow](https://bioconductor.uib.no/packages/3.18/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html).
- Primary raw reconstructed provenance: Actions [#44907357144](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44907357144) SUCCESS, 204 DCC and PKC.
- Final pinned source QC: Actions [#44907538365](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44907538365), SUCCESS, artifact `hs-geomx-2026-loq-q3-q90-fast-qc`.
- Source code: `projects/hs_spatial_gse294009/reconstruct_roi_counts.py` and `projects/hs_spatial_gse294009/loq_q3_robustness.py`.
- Target: S2-mapped 132 HS GEO ROIs (15 HS donors). 14 acne conglobata ROIs stem from one donor and were not used as independent disease controls.
- Primary contrast: same-site Yonsei FFPE, HS12 excluded: up to 7 Early donors and 3 Late donors; within each compartment, donor medians, not ROI-level replicates. The other 24 historical S2 ROIs remain quantitatively unavailable.

## Method

The probe-class mapping comes from the published PKC; raw RTS counts are read from 204 GEO DCC files. Probe controls were identified from PKC `CodeClass`.

1. Q3: per-ROI 75th percentile across endogenous genes, denominator for the exploratory gene ratio.
2. Q90: 90th percentile across the same genes, alternative denominator.
3. Background: from every ROI's negative probes, clamp zeros to count=1 before log, compute **geometric mean** (GM) and **sample geometric standard deviation** (GSD). Exploratory threshold `max(2, GM × GSD²)`. **Not a confirmed exact reproduction of the authors' LOQ**, which may have different probe/outlier filtering.
4. Q3_LOQ×1 and ×2: set gene ratios below `1×` or `2×` exploratory LOQ to zero in a sensitivity analysis. A numerical contrast of zero after censoring **does not mean biological equivalence**.
5. Thresholds are compared only for genes mapping to exactly one RTS probe, because summing multiple gene probes would make a per-probe negative-control LOQ invalid.
6. Signal separation is measured as `Q3 / negative_probe_geomean`. Exploratory thresholds 3× and 5× generate separate sensitivity subsets; these are **not manufacturer-defined pass/fail limits**.
7. Donor/compartment median ratio contrast is `log₂((Late + 0.1)/(Early + 0.1))`. No p-values, FDR or official normalization claims.
8. Candidate remains an **exploratory validation candidate** only if effects are direction-consistent across Q3, Q90, LOQ×1 and LOQ×2 and if median patient-level LOQ detection fraction ≥0.5 in each stage, with at least 3 Early and 2 Late donors. This is an explicit heuristic, **not a prevalidated biological significance cutoff**.

## Original negative-control profile

All 132 S2-linked HS ROIs support the exploratory LOQ computation. Among the same-site primary subset:

| Quantity | Early, 57 ROIs, 7 donors | Late, 25 ROIs, 3 donors |
|---|---:|---:|
| Median negative-probe count per ROI | 119 | 121 |
| Median endogenous Q3 | 13.0 | 10.0 |
| Median negative-probe geometric mean | 2.3848 | 2.4495 |
| Median exploratory LOQ | 9.2437 | 9.3635 |
| Median Q3 / neg GM | 4.4955 | 3.2891 |
| ROI with Q3 / neg GM below 3 | 14 / 57 | 10 / 25 |
| ROI with Q3 / neg GM below 5 | 37 / 57 | **25 / 25** |

**Crucial:** at the stricter exploratory Q3/negative-GM ≥5 threshold, **zero Late ROIs remain**. The 5× subset cannot support any Early–Late comparison. Late-site background separation is weaker and may substantially influence the apparent contrasts.

## Deep dermal inflammatory aggregate: key genes

Effects are log₂(Late/Early); the published study's actual LOQ-normalized expression was not reproduced.

| Gene | Q3 | Q90 | Below LOQ×1 censored | Below LOQ×2 censored | Median donor LOQ detection (Early/Late) | Decision |
|---|---:|---:|---:|---:|---|---|
| **JCHAIN** | +1.8437 | +1.8662 | +2.9506 | +4.6049 | 0.50 / 1.00 | **RETAIN_EXPLORATORY_FOR_ORTHOGONAL_VALIDATION** |
| CXCL13 | −0.9355 | −0.5768 | −1.4483 | 0.0000 | 0.6667 / 0.50 | HOLD_QC_SENSITIVE_DIRECTION |
| MS4A1 | +0.6630 | +0.8301 | 0.0000 | 0.0000 | 0.00 / 0.00 | HOLD_LOW_DETECTABILITY |
| CD79A | +0.2035 | +0.1320 | −0.2370 | 0.0000 | 0.6667 / 0.50 | HOLD_QC_SENSITIVE_DIRECTION |

Across **24 predefined gene × compartment contrasts**, **1 retained** after the above exploratory screen (JCHAIN in deep dermis). **23 are HOLD** on low detectability, sensitivity to the threshold, or insufficient donor support.

JCHAIN is **not a validated biomarker**. Its effect inflates when low counts are zero-censored. It may be related to late-stage plasma cell programs already discussed by the original paper, rather than a novel finding. It should be checked by patient-level in situ RNA or protein (e.g. matched immunohistology) and independent cohorts.

## Scientific decision and next work

- **GO**: nominate JCHAIN as a targeted *technical verification* candidate, while carrying both LOQ detection fractions and donor sample sizes.
- **HOLD**: claims that CXCL13, CD79A or MS4A1 are robust stage biomarkers in this dataset. The gene counts at/near the exploratory LOQ cannot be interpreted as strong negative or positive biological evidence.
- **HOLD**: applying a stringent Q3/negative-GM ≥5 filter then comparing Late vs Early (no Late ROI survives).
- **HOLD**: confirmatory differential expression, valid adjusted p/FDR, prognostic or drug-response prediction, and full 156-ROI replication.
- **Priority next**: cross-check the authors' Supplementary Methods exactly for their probe QC, LOQ and normalization choices; consider independent molecular validation and non-ROI biological replicates.

The CI SUCCESS is the **success of a reproducible audit**, not successful validation of biomarkers or official GeoMx normalization.
