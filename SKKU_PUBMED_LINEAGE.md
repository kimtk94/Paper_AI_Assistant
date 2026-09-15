# SKKU PubMed Research Lineage

## Goal

Build a reproducible PubMed-based map of how SKKU-associated research continues over time.

The workflow separates two questions:

1. **SKKU publication graph**  
   Find papers whose PubMed author affiliation contains Sungkyunkwan / SKKU, then connect papers by:
   - direct citation/reference relationships;
   - repeated SKKU-affiliated authors;
   - MeSH / keyword similarity.

2. **Researcher continuation graph**  
   Use ORCID found in the SKKU seed papers to follow the same researcher in PubMed even after their affiliation changes.

This distinction matters because a paper-level SKKU affiliation search alone will stop following a researcher after they move to another institution.

---

## Data source

Primary source: **PubMed / NCBI E-utilities**

Important PubMed fields and links used by this project:

- Affiliation: \`[ad]\`
- Author Identifier / ORCID: \`[auid]\`
- Publication Date: \`[dp]\`
- References: \`pubmed_pubmed_refs\`
- Cited-by: \`pubmed_pubmed_citedin\`

The scripts include NCBI \`tool\` and \`email\` parameters and throttle requests. If \`NCBI_API_KEY\` is supplied, a faster safe request interval is used.

---

## Files

### 1. SKKU paper graph

\`src/skku_pubmed_lineage.py\`

Example:

    export NCBI_EMAIL="your_email@example.com"

    python src/skku_pubmed_lineage.py \
      --start-year 2018 \
      --end-year 2026 \
      --max-results 500 \
      --output-dir outputs/skku_pubmed_lineage

Optional topic filter:

    python src/skku_pubmed_lineage.py \
      --start-year 2018 \
      --end-year 2026 \
      --topic "stroke OR cerebrovascular" \
      --max-results 500

Outputs:

- \`papers.json\`: full parsed PubMed records
- \`papers.csv\`: flat table for Excel / pandas
- \`edges.json\`: graph edges
- \`edges.csv\`: flat edge table
- \`lineage.md\`: readable chain summary
- \`lineage.html\`: interactive network visualization
- \`summary.json\`: run summary

### Edge meaning

| Relation | Weight | Meaning |
|---|---:|---|
| citation | 1.00 | Direct PubMed reference / cited-by relationship |
| shared_skku_author | 0.95 ORCID / 0.80 name | Same SKKU-affiliated researcher appears in both papers |
| topic_similarity | Jaccard score | MeSH + keyword overlap |

The graph is directed from an older paper to a newer paper where possible.

---

## 2. Follow the researcher after SKKU

\`src/skku_pubmed_author_followup.py\`

This script consumes the first script's \`papers.json\`.

High-confidence mode:

    python src/skku_pubmed_author_followup.py \
      --seed-json outputs/skku_pubmed_lineage/papers.json \
      --start-year 2002 \
      --end-year 2026 \
      --max-authors 100 \
      --max-per-author 200 \
      --with-citations \
      --output-dir outputs/skku_pubmed_followup

Outputs:

- \`author_registry.csv\`: researchers detected from SKKU seed papers
- \`lineage_papers.csv\`: all PubMed papers followed for those researchers
- \`lineage_papers.json\`
- \`continuation_edges.csv\`
- \`continuation_edges.json\`
- \`continuation.md\`
- \`summary.json\`

### Identity confidence

**High confidence**
- PubMed ORCID / Author Identifier exact match.

**Low confidence**
- Exact full-author-name match.
- Disabled by default.
- Enable only when needed:

    --allow-name-fallback

Name-only results require manual review because author-name collisions are possible.

---

## Recommended interpretation

A useful research lineage is not simply "papers with similar titles."

Priority should be:

1. **Direct citation**
2. **Same ORCID researcher**
3. **Same researcher name, manually reviewed**
4. **MeSH / keyword similarity**

This lets the output distinguish a defensible research continuation from a loose topical resemblance.

---

## Important limitation: "SKKU affiliation" vs "SKKU alumni"

PubMed records publication-time author affiliations. They do **not** establish whether someone graduated from SKKU.

Therefore:

- the first stage is an **SKKU-affiliation publication map**;
- the ORCID stage is a **researcher continuation map**;
- a true **SKKU alumni map** needs a separate verified alumni/researcher registry.

A later extension can add:

    config/skku_researchers.csv

with fields such as:

    name,orcid,skku_status,school_or_department,degree_year,verification_source

Then the same PubMed tracking code can follow only verified alumni.

---

## Recommended next extension

For practical literature intelligence, add three labels to each chain:

- **Research theme**: derived from MeSH / keywords
- **Method evolution**: e.g. GWAS -> MR -> multi-omics -> single-cell
- **Disease / phenotype**: normalized topic labels

Then a chain can be displayed as:

    Paper A (2019)
      -> Paper B (2021)
      -> Paper C (2023)
      -> Paper D (2026)

with the reason for each transition explicitly shown as citation, same researcher, or topic continuation.
