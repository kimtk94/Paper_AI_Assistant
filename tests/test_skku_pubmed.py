import json
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import skku_pubmed_lineage as lineage
import skku_pubmed_author_followup as followup


PUBMED_XML = """<PubmedArticle>
<MedlineCitation>
  <PMID>12345678</PMID>
  <Article>
    <ArticleTitle>Example multi-omics study.</ArticleTitle>
    <Abstract>
      <AbstractText Label="BACKGROUND">Background text.</AbstractText>
      <AbstractText Label="METHODS">Methods text.</AbstractText>
    </Abstract>
    <Journal>
      <Title>Example Journal</Title>
      <JournalIssue>
        <PubDate><Year>2026</Year><Month>Sep</Month><Day>15</Day></PubDate>
      </JournalIssue>
    </Journal>
    <AuthorList>
      <Author>
        <LastName>Kim</LastName>
        <ForeName>Tae Hoon</ForeName>
        <Initials>TH</Initials>
        <Identifier Source="ORCID">0000-0001-2345-6789</Identifier>
        <AffiliationInfo>
          <Affiliation>Department of Digital Health, Sungkyunkwan University, Seoul, Republic of Korea.</Affiliation>
        </AffiliationInfo>
      </Author>
      <Author>
        <LastName>Lee</LastName>
        <ForeName>Example</ForeName>
        <Initials>E</Initials>
        <AffiliationInfo>
          <Affiliation>Other University, Seoul, Republic of Korea.</Affiliation>
        </AffiliationInfo>
      </Author>
    </AuthorList>
    <PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
  </Article>
  <MeshHeadingList>
    <MeshHeading><DescriptorName>Stroke</DescriptorName></MeshHeading>
    <MeshHeading><DescriptorName>Multiomics</DescriptorName></MeshHeading>
  </MeshHeadingList>
  <KeywordList>
    <Keyword>single cell</Keyword>
  </KeywordList>
</MedlineCitation>
<PubmedData>
  <ArticleIdList>
    <ArticleId IdType="pubmed">12345678</ArticleId>
    <ArticleId IdType="doi">10.1000/example</ArticleId>
    <ArticleId IdType="pmc">PMC123456</ArticleId>
  </ArticleIdList>
</PubmedData>
</PubmedArticle>"""


def parsed_paper():
    return lineage.parse_paper(ET.fromstring(PUBMED_XML))


class TestPubMedParsing(unittest.TestCase):
    def test_parse_paper_extracts_core_fields(self):
        paper = parsed_paper()
        self.assertEqual(paper.pmid, "12345678")
        self.assertEqual(paper.year, 2026)
        self.assertEqual(paper.doi, "10.1000/example")
        self.assertEqual(paper.pmcid, "PMC123456")
        self.assertEqual(paper.journal, "Example Journal")
        self.assertIn("Stroke", paper.mesh_terms)
        self.assertIn("single cell", paper.keywords)
        self.assertIn("BACKGROUND: Background text.", paper.abstract)

    def test_skku_author_and_orcid_are_detected(self):
        paper = parsed_paper()
        self.assertEqual(paper.skku_authors, ["Tae Hoon Kim"])
        self.assertEqual(paper.skku_orcids, ["0000-0001-2345-6789"])
        self.assertEqual(len(paper.skku_affiliation_evidence), 1)
        self.assertIn("Sungkyunkwan University", paper.skku_affiliation_evidence[0])

    def test_query_contains_affiliation_date_and_topic(self):
        query = lineage.build_query(
            start_year=2020,
            end_year=2026,
            topic="stroke OR cerebrovascular",
            extra="humans[mh]",
        )
        self.assertIn("Sungkyunkwan[ad]", query)
        self.assertIn("SKKU[ad]", query)
        self.assertIn("2020:2026[dp]", query)
        self.assertIn("stroke OR cerebrovascular", query)
        self.assertIn("humans[mh]", query)


class NoNetworkClient:
    def links(self, pmid, linkname):
        raise AssertionError("network link lookup should not run in this test")


class TestGraphLogic(unittest.TestCase):
    def test_shared_author_and_topic_edges(self):
        p1 = parsed_paper()
        p1.pmid = "100"
        p1.year = 2024
        p1.title = "Earlier paper"

        p2 = parsed_paper()
        p2.pmid = "200"
        p2.year = 2026
        p2.title = "Later paper"
        p2.mesh_terms = ["Stroke", "Multiomics", "Humans"]

        edges = lineage.build_edges(
            papers=[p1, p2],
            client=NoNetworkClient(),
            include_citations=False,
            topic_threshold=0.20,
            max_pairs=100,
        )
        relations = {e.relation for e in edges}
        self.assertIn("shared_skku_author", relations)
        self.assertIn("topic_similarity", relations)

        shared = next(e for e in edges if e.relation == "shared_skku_author")
        self.assertEqual((shared.source, shared.target), ("100", "200"))
        self.assertEqual(shared.weight, 0.95)
        self.assertIn("ORCID=", shared.evidence)


class TestAuthorContinuation(unittest.TestCase):
    def test_registry_prefers_orcid(self):
        paper = parsed_paper()
        seed = [lineage.asdict(paper)]
        registry = followup.build_registry(seed, allow_name_fallback=False)
        self.assertEqual(len(registry), 1)
        self.assertEqual(registry[0].orcid, "0000-0001-2345-6789")
        self.assertEqual(registry[0].confidence, "high")
        self.assertEqual(registry[0].seed_pmids, ["12345678"])

    def test_orcid_author_query_uses_auid(self):
        tracked = followup.TrackedAuthor(
            key="orcid:0000-0001-2345-6789",
            name="Tae Hoon Kim",
            orcid="0000-0001-2345-6789",
            confidence="high",
            seed_pmids=["12345678"],
        )
        query = followup.author_query(tracked, 2020, 2026)
        self.assertEqual(
            query,
            '"orcid 0000-0001-2345-6789"[auid] AND 2020:2026[dp]',
        )

    def test_matching_author_by_orcid(self):
        paper = parsed_paper()
        tracked = followup.TrackedAuthor(
            key="orcid:0000-0001-2345-6789",
            name="Tae Hoon Kim",
            orcid="0000-0001-2345-6789",
            confidence="high",
            seed_pmids=[],
        )
        author = followup.matching_author(paper, tracked)
        self.assertIsNotNone(author)
        self.assertEqual(author.name, "Tae Hoon Kim")



class TestResearchLineageClassification(unittest.TestCase):
    def setUp(self):
        self.key = "orcid:0000-0001-2345-6789"
        self.registry = {
            self.key: followup.TrackedAuthor(
                key=self.key,
                name="Tae Hoon Kim",
                orcid="0000-0001-2345-6789",
                confidence="high",
                seed_pmids=["100"],
            )
        }
        self.paper_authors = {
            "100": {self.key},
            "200": {self.key},
            "300": {self.key},
        }

    def test_direct_citation_between_consecutive_papers_is_strongest(self):
        edges = [
            followup.ContinuationEdge(
                source="100",
                target="200",
                relation="author_continuation",
                tracked_author="Tae Hoon Kim",
                author_key=self.key,
                confidence="high",
                evidence="ORCID 0000-0001-2345-6789",
            ),
            followup.ContinuationEdge(
                source="100",
                target="200",
                relation="citation",
                tracked_author="",
                author_key="citation",
                confidence="high",
                evidence="200 cites 100",
            ),
        ]
        lineage_edges = followup.build_research_lineage_edges(
            edges, self.paper_authors, self.registry
        )
        self.assertEqual(len(lineage_edges), 1)
        self.assertEqual(lineage_edges[0].relation, "direct_citation_continuation")
        self.assertEqual(lineage_edges[0].score, 1.0)

    def test_orcid_only_continuation_is_distinguished_from_citation(self):
        edges = [
            followup.ContinuationEdge(
                source="100",
                target="200",
                relation="author_continuation",
                tracked_author="Tae Hoon Kim",
                author_key=self.key,
                confidence="high",
                evidence="ORCID 0000-0001-2345-6789",
            )
        ]
        lineage_edges = followup.build_research_lineage_edges(
            edges, self.paper_authors, self.registry
        )
        self.assertEqual(lineage_edges[0].relation, "author_continuation_only")
        self.assertEqual(lineage_edges[0].score, 0.75)

    def test_nonconsecutive_self_citation_is_preserved(self):
        edges = [
            followup.ContinuationEdge(
                source="100",
                target="200",
                relation="author_continuation",
                tracked_author="Tae Hoon Kim",
                author_key=self.key,
                confidence="high",
                evidence="ORCID 0000-0001-2345-6789",
            ),
            followup.ContinuationEdge(
                source="200",
                target="300",
                relation="author_continuation",
                tracked_author="Tae Hoon Kim",
                author_key=self.key,
                confidence="high",
                evidence="ORCID 0000-0001-2345-6789",
            ),
            followup.ContinuationEdge(
                source="100",
                target="300",
                relation="citation",
                tracked_author="",
                author_key="citation",
                confidence="high",
                evidence="300 cites 100",
            ),
        ]
        lineage_edges = followup.build_research_lineage_edges(
            edges, self.paper_authors, self.registry
        )
        relations = {e.relation for e in lineage_edges}
        self.assertIn("direct_citation_same_author", relations)
        direct = next(
            e for e in lineage_edges
            if e.relation == "direct_citation_same_author"
        )
        self.assertEqual(direct.score, 0.95)
        self.assertEqual(direct.tracked_authors, ["Tae Hoon Kim"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
