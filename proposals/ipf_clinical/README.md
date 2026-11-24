# IPF clinical biomarker collaboration feasibility audit

Read `RESEARCH_PROTOCOL.md` before use. This is a scientific proposal checkpoint, **not a prognostic model**.

`python -m unittest discover -s proposals/ipf_clinical -p 'test_*.py'`

Clinical data are private. Never commit or upload the manifest to GitHub, Actions artifacts, or public Drive folders. For local authorized review only:

`python proposals/ipf_clinical/feasibility.py --manifest /secure/deidentified_manifest.csv --out /secure/aggregate_audit.json`

Inputs must include all fields in `feasibility.COLUMNS`. Use only IRB-approved deidentified data. Output contains aggregate QC counts, without patient identifiers or clinical measurements. GO states mean only that exploratory recruitment/checklist criteria have been met; no assumption of assay detection or statistical power.

The synthetic CI uses no actual patient information, no external connection, and cannot trigger email, SFDC write, or quotation.