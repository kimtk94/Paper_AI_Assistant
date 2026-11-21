# Migraine public external-data feasibility (2026-11-21)

Primary target: Yonsei serum SCFA work, doi:10.1186/s10194-026-02268-4, 476 participants; metabolite values from original participants not available in this audit.

Public external study: MTBLS13630, doi:10.1016/j.isci.2026.115587, 30 PFO-migraine + 17 healthy controls; **not the same disease subtype** or treatment response.

Run: `python projects/migraine_mtbls13630/audit.py`. The workflow reads public metadata only, caps responses at 4MB and does not store raw data or sample-level clinical information. `HOLD` is the appropriate status until numerical LC-MS matrix and unambiguous clinical group mapping are validated.

2026 UK Biobank proteomics / NMR migraine multiomics (doi:10.1186/s10194-026-02493-x) makes generic plasma migraine association discovery low-novelty. Prioritize longitudinal preventive-treatment response, paired serum and headache diary, and validated SCFA with broad metabolite discovery, not cross-sectional clinical classifiers.