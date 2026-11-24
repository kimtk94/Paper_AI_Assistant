# Public data analysis: GSE283755 / endometriosis + Dienogest

## Scientific question and input
- Seok Kyo Seo et al., Yonsei Med J 2025 DOI: [10.3349/ymj.2024.0518](https://doi.org/10.3349/ymj.2024.0518), PMID 41145313.
- GEO [GSE283755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283755) / PRJNA1195352.
- 4 preoperative dienogest-exposed ovarian endometrioma tissues vs 4 untreated tissues **from different patients**.
- Public processed Excel `GSE283755_Expression_Profile.GRCh38.gene.xlsx` (11.9MB); SRA raw reads separately available.
- Original study reports 406 differentially expressed genes by DESeq2 and pathway enrichment.

## Critical boundaries
- This is an observational exposed/unexposed tissue comparison, **not a randomized intervention** or paired pretreatment/response study.
- If Excel contains only TPM/FPKM, DESeq2 must not be run using those values.
- The Python Welch + BH analysis, if sample values are available, is explicitly exploratory, **not a reproduction of DESeq2 or publication-grade inferential analysis**.
- Public atlas [GSE266265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265) contains 60 **eutopic endometrial** biopsies (35 endometriosis / 25 controls) collected without hormonal treatment, NOT the same tissue/site/treatment context. Contextualize cell type only, not an external Dienogest validation.

## Execution
`pip install -r projects/endometriosis_gse283755/requirements.txt`
`python projects/endometriosis_gse283755/run.py --output output/endometriosis_gse283755`

Workflow `.github/workflows/endometriosis-geo-audit.yml` has `pull_request` and `workflow_dispatch` triggers. All outputs retained as CI artifact. Raw GEO file is never committed.

## Gates
1. GEO file successfully downloaded, SHA256 recorded.
2. Original samples/columns mapped to the **8 published GSM IDs**; if ambiguous, `BLOCKED_MATRIX_MAPPING` and only column inventory generated.
3. Matrix QC, exploratory effect directions, PCA, targeted LRRC15 direction only if valid numeric sample columns.
4. DESeq2 reproduction only after raw counts confirmed; compare DEG overlap and direction with original Supplementary Table 2.
5. Cell-of-origin mapping with properly harmonized GSE266265 after raw-data evidence and tissue-context caveat.
6. Future clinical proteomic validation requires fresh pretreatment plasma and prospective clinical response, cannot be inferred from the GEO study.

No mail, SFDC, or external outreach actions are performed by this workflow.