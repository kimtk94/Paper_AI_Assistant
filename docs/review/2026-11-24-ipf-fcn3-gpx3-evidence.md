# IPF FCN3 / GPX3 evidence review — 2026-11-24

**Decision: ARCHIVED_EVIDENCE_RECHECK__HOLD_CLINICAL_PROGNOSIS.**

Scope: public Finnish IPF/control plasma PXD010965 and independent lung RNA cohorts GSE135893 / GSE136831. They are not paired cohorts; neither disease-status association nor lung RNA establishes circulating-protein origin or 12-month FVC prediction.

## Verified input provenance

All five retrieved artifact ZIPs matched the GitHub API SHA-256. Four pinned inputs are used by the new offline replay:

| Dataset | Run | Artifact | Expiration (UTC) |
|---|---:|---:|---|
| PXD010965 | [45375552475](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375552475) | 14909794711 | 2027-01-23 10:32:25 |
| GSE135893 donor/source | [45375552102](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375552102) | 14910154655 | 2026-12-24 10:31:20 |
| GSE135893 donor cell types | [45375552581](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375552581) | 14911122272 | 2026-12-24 10:34:40 |
| GSE136831 donor cell types | [45375552605](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375552605) | 14912096747 | 2026-12-24 10:37:52 |
| GSE136831 metadata only (checked, not replay input) | [45375552346](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375552346) | 14909849236 | 2026-12-24 10:31:23 |

Replay validates the exact ZIP/member identities and the two matrix-entry-count gates, independently recomputes plasma log2 effects and leave-one-assay-out direction, and recomputes donor-level tissue medians. FDR values below are from the archived prior Welch/BH calculation, not newly recalculated. The source-file producer and runtime are unchanged. The original large MTX files have **not** been streamed again in this review.

## Plasma normalization sensitivity

PXD010965: 167 unique proteins; explicit CONTROL 20 and CASE 17 abundance columns; 6,179 positive numeric values. Source SHA-256: 1df37c4f3a2e130bd02ea87dc5da47ba1a6701efe24c657b358cc721bbd3f334.

| Gene | Raw log2FC | Raw archived FDR | Median-centered log2FC | Centered archived FDR | Centered LOO range | Direction |
|---|---:|---:|---:|---:|---|---|
| FCN3 | -0.845820 | 3.249e-7 | -0.337395 | 0.00222343 | -0.358703 to -0.305892 | 37/37 |
| GPX3 | -2.101865 | 2.496e-5 | -1.593439 | 0.00127845 | -1.703662 to -1.478060 | 37/37 |

Both also maintain negative raw log2FC in all 37 leave-one-assay-out folds. Per-sample median centering is a sensitivity analysis, **not** the original Progenesis normalization. Whole-proteome archived FDR<0.05 changes from 128 to 57, overlap 42, with group median intensity difference p=6.56e-7. No longitudinal endpoint, sample-level age/sex, smoking or treatment covariates in the processed table; no causal or clinical prediction claim. GPX3 is a reasonable candidate to prioritize for independent assay development because of its larger effect under both normalization choices; FCN3 remains a mechanistic candidate.

## Donor-level lung RNA triangulation

GSE135893 has 12 IPF / 10 control donors overall, and site confounding (NTI supplies IPF, DNA supplies controls; only Vanderbilt has both, 7/7 overall). GSE136831 metadata has 32 IPF / 28 control donors (plus excluded COPD). Their archived matrix gates passed. This review applies analyst-defined **20 cells per donor/celltype** and **5 donors per group** as conservative descriptive reporting thresholds, not a power calculation. All numbers below are log2[(median IPF CPM+1)/(median control CPM+1)], not inferential differential expression.

| Dataset and cell label | Gene | IPF/control donors | Descriptive log2 ratio |
|---|---|---:|---:|
| GSE135893 Vanderbilt endothelial | FCN3 | 6/6 | -1.63628 |
| GSE135893 Vanderbilt endothelial | GPX3 | 6/6 | -1.44626 |
| GSE135893 Vanderbilt fibroblasts | GPX3 | 1/2 | HOLD sparse |
| GSE135893 Vanderbilt myofibroblasts | GPX3 | 3/4 | HOLD sparse |
| GSE136831 Endothelial category | FCN3 | 28/19 | -1.49913 |
| GSE136831 Endothelial category | GPX3 | 28/19 | -0.26292 |
| GSE136831 Stromal category | GPX3 | 27/13 | -2.44475 |
| GSE136831 Fibroblast subclass | GPX3 | 21/10 | -1.90296 |
| GSE136831 Epithelial category | FCN3 | 32/20 | +0.54741 |
| GSE136831 Epithelial category | GPX3 | 32/20 | +0.01468 |

FCN3's endothelial direction is consistent across cohorts, but epithelial direction is not uniformly negative. GPX3's endothelial effect is weaker in the second cohort; its stromal/fibroblast signal is interesting but source-matched GSE135893 fibroblast donors are sparse. **No cross-cohort pooled estimate or formal replication significance.** The labels are not exactly harmonized. Lung transcript data are not paired plasma protein measurements.

## Reproduction and release gates

Place the four original downloaded ZIPs, with filenames github-actions-artifact-<id>.zip, in a local authorized directory. Then run:

```sh
python projects/ipf_evidence_review/replay.py --artifact-dir /path/to/archives --out /tmp/ipf-review.json
python -m unittest discover -s tests -p 'test_ipf_evidence_review.py' -v
```

Actual offline replay completed using retrieved artifact bytes; seven synthetic tests passed locally. The repository receives replay code and aggregate facts only, not assay rows or donor-level files. Preserve original archives only in authorized storage before expiration; checksum records alone do not preserve source bytes.

Five external-data workflows are now on-demand to avoid downloading large matrices on routine PR changes. A lightweight offline PR test workflow replaces this routine check; existing scientific source gates remain intact. On-demand workflows may need to exist on the default branch before manual dispatch is available.

Next gate: independent targeted protein assay validation with prespecified baseline collection, FVC outcome, follow-up window, covariates, and calibration plan; confirm actual eligible participant counts, preanalytics and protocol (#24), then evaluate independent predictive performance. No new FVC endpoint, clinical cohort, patient data, or outreach was introduced.
