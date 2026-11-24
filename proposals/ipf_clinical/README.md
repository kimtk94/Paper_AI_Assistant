# IPF clinical biomarker collaboration feasibility audit

Read `RESEARCH_PROTOCOL.md` before use. This is a scientific proposal checkpoint, **not a prognostic model**.

`python -m unittest discover -s proposals/ipf_clinical -p 'test_*.py'`

Clinical data are private. Never commit or upload the manifest to GitHub, Actions artifacts, or public Drive folders. For local authorized review only:

`python proposals/ipf_clinical/feasibility.py --manifest /secure/deidentified_manifest.csv --out /secure/aggregate_audit.json`

Inputs must include all fields in `feasibility.COLUMNS`. Use only IRB-approved deidentified data. Output contains aggregate QC counts, without patient identifiers or clinical measurements. GO states mean only that exploratory recruitment/checklist criteria have been met; no assumption of assay detection or statistical power.

The synthetic CI uses no actual patient information, no external connection, and cannot trigger email, SFDC write, or quotation.
## Updated clinical design and v2 screening

Read [the 2026-11-24 design and readiness gates](CLINICAL_DESIGN_20261124.md). The legacy screen is preliminary only and excludes every row for a duplicated key. Use the v2 screen for explicit diagnosis, same-person linkage, specimen and endpoint adjudication. Both scripts are local-only. To generate planning assumptions without any patient data: `python proposals/ipf_clinical/study_readiness.py --out /secure/planning.json`. A v2 manifest must contain `study_readiness.V2_COLUMNS`; outcome codes and acceptance criteria are in the design and code. No output authorizes a study or verifies clinical validity.
