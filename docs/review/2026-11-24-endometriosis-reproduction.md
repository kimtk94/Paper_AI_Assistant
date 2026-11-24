# GSE283755 reproduction review — 2026-11-24

**Decision: partial same-cohort reproduction; HOLD exact 406 DEG set and clinical prediction.**

Source: Chung et al., Yonsei Med J 2025; DOI https://doi.org/10.3349/ymj.2024.0518 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12568438/. The publication describes independent groups of four each, DESeq2 1.36.0, exclusion of genes with zero counts in more than three samples of either group, BH FDR <0.05 and twofold-change criterion. No batch correction was used (same experimental batch).

Audit inputs: run [45375547987](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375547987), artifact 14911890483, ZIP SHA256 2aecc6282820ad31ca49d104d5951864b789a93e8bb910e476374232ac420a3b; run [45375548423](https://github.com/kimtk94/Paper_AI_Assistant/actions/runs/45375548423), artifact 14912260210, ZIP SHA256 73ff1f2aa38fa9854820b1b12cc2d3df0abaaf55acd321419f46ccbf33aaaa71. Both ZIP digests were independently recalculated from retrieved bytes. Workbook input SHA256 reported by the producer: 0cff5c06c75a897b8fe16475226adf9dbe0c538a6ebc60b5eae431b9ccb5de2f. Workbook not included in this commit. Audited artifacts expire 2026-12-24.

| Metric | Independent prior run |
|---|---|
| Gene-symbol aggregated features | 46,425 |
| Conservative count filter | 20,602 features; 418 DEGs |
| Publication-like zero filter | 27,518 features; 421 DEGs (73 up, 348 down) |
| FDR-only count, publication-like filter | 503 |
| DESeq2 runtime | 1.46.0, not 1.36.0 |
| LOO DEG range (8 refits) | 58–1,424 |
| Published top five effects | 5/5 agree to two printed decimals; all 5 pass current threshold in all 8 LOO refits |

Table 3 anchors: HLA-DOB, ELOVL2, NSFP1, GSTM1, LRRC15. The five matching effects are informative because the same samples were used, but are not external replication. No exact 406-list intersection, pathways, or patient-level pretreatment response validation is asserted.

Official Supplementary Table 2 is listed by the article. Earlier download audit artifact 14910638208 (SHA256 41bfd98ec03b7d4c97f5b094532494d95c69be54847f2f140d915f793254bcc3) reports non-PDF/404 across six endpoints; this is an access failure, not proof the data are absent. Add the official Europe PMC supplement archive route with unique exact PDF member validation. Preserve only provenance metadata in CI artifacts, not PDF or extracted text. An available PDF still requires complete schema and 406-gene identity checks. The article is CC BY-NC; verify permitted use before redistribution or commercial use.

Next gates: (1) verify complete official Table 2 and compare unique gene identifiers, signed effects, thresholds and set overlap; (2) perform pinned DESeq2 1.36.0 run with the original annotation and feature definition (GEO includes Gene_ID and Gene_Symbol; aggregation can change the testing universe); (3) compare the gene-ID and gene-symbol matrices as explicitly labeled sensitivity runs; (4) plan independent donor validation before any predictive or therapeutic claim. Never choose settings simply to obtain 406.

Code changes in this review: exact public Read_Count header classification, conflicting sample-ID rejection, strict finite nonnegative integer validation before feature filtering, explicit gene-symbol aggregation row accounting, equally preferred matrix ambiguity block, official archive fallback, top-five comparison and synthetic-only regression tests. Prior audited numerical results precede these changes; no new full reanalysis is implied.
