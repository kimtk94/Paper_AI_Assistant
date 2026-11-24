# Migraine source and feasibility review — 2026-11-24

**Current decision: HOLD sample-linked quantitative values and participant/phenotype mapping.**

## Two different studies

- [Cho et al., J Headache Pain (2026)](https://doi.org/10.1186/s10194-026-02268-4): 476 people (190 episodic migraine, 178 chronic migraine, 108 healthy controls), serum targeted acetate/propionate/butyrate LC-MS/MS, cross-sectional associations with preventive treatment; participant data available on reasonable request from the corresponding author. No longitudinal predictive result from this study. Article license CC BY-NC-ND: check permissions for intended use.
- [Luo et al., iScience (2026)](https://doi.org/10.1016/j.isci.2026.115587): 30 PFO-associated migraine cases (arterial and venous plasma) and 17 healthy controls (venous plasma); reported 211 differential metabolites and an independent 220-person indoleacrylic acid cohort. These are **published results, not newly reproduced here**. Raw data statement points to MTBLS13630; original code not published. Not the same disease subgroup, sample matrix, participants or outcome as Cho. No direct validation or pooling.

## Verified earlier artifact

[Run 44876949912](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/44876949912), artifact 14696416712, ZIP SHA-256 84e683fe8d83696049f5944e6607dc44e75e107f8db3616f08eb929436b4f22f, expiration 2027-01-20 06:44:37 UTC. The retrieved bytes matched the GitHub API digest.

- 91 ISA rows (47 venous, 30 arterial, 14 QC); all 91 sample and source names distinct. This is not evidence of 91 distinct people. Variant field is blank in all 91 records. No verified clinical phenotype or arterial/venous subject crosswalk.
- GC-MS 91, LC-MS positive 90, LC-MS negative 90 records, 271/271 matched to sample sheet. Derived spectral file absent from all 271 records. Raw spectral file paths are present.
- Three MAF API projections: 37, 92, 82 annotation rows (211 total). 21 fields each, sample_abundance flag false; aggregate abundance fields blank; zero numeric abundance cells checked. A sample_abundance=true query did not add columns. This establishes only the behavior of the API responses; it does not prove original TSVs or repository subdirectories lack quantitative values.
- Root listing contained the standard ISA metadata, three MAF files, FILES/, HASHES/ and METADATA_REVISIONS/. Prior audit did not inspect FILES/ or direct original MAF content. The [official MAF guide](https://ebi-metabolights.github.io/guides/MAF/) describes quantitative sample columns to the right of annotation fields.

## New bounded direct-source audit

The new code reads original TSVs, requires exact sample/assay identities, inspects numeric cells only in exact sample-name columns (never infers columns by position), distinguishes aggregate annotation values, rejects ambiguous headers and malformed rows, inventories immediate FILES/ entries, and inspects the official article supplementary ZIP member inventory without extracting files. SHA-256 and access status are recorded. No raw instrument data, patient-level matrix, PDF or spreadsheet is downloaded from the listed file paths. A ZIP containing supplementary files is retrieved only for bounded inventory. Output is aggregate audit only, not a clinical data export. An inaccessible endpoint or schema mismatch means HOLD, not definitive source absence. Even a candidate numeric table does not authorize case-control analysis.

## Next approval conditions

1. Acquire original feature-by-sample quantitative table (separate GC-MS, LC-MS positive and negative mode) with feature annotations, QC/blanks, units and normalization definitions, or plan a separately authorized raw-instrument reprocessing project.
2. Obtain an authoritative deidentified participant/venous/arterial/QC crosswalk and explicit PFO/migraine/control labels, not guessed from filenames. Use venous-only case-control analysis and paired within-case arterial/venous contrast separately; do not pool as independent people.
3. Obtain the independent IA validation cohort data only if needed and permitted; do not claim the reported association or AUC is replicated by metadata.
4. For the actual Yonsei question, request the original SCFA values and covariates or plan a new prospective paired baseline/follow-up cohort with headache diary, medication exposure, preanalytical SOP and independent validation. Distinguish serum from plasma and targeted quantitative assays from untargeted discovery.
5. Preserve authorized source evidence before expiration. No author contact was sent, no private clinical data were retrieved, no clinical prediction or novel validated biomarker is asserted.

Prior numerical claims above come from the archived audit or publications, not from newly generated patient values. Live results from the new audit should be documented separately after the workflow completes.
