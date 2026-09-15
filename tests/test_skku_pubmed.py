import json
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import skku_pubmed_lineage as lineage
import skku_pubmed_author_followup as followup
import skku_pubmed_researcher_network as researcher_network
import skku_pubmed_research_communities as research_communities


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



class TestResearchProfileAnnotation(unittest.TestCase):
    def test_gwas_mr_stroke_profile(self):
        paper = parsed_paper()
        paper.title = "Genome-wide association study identifies stroke risk variants"
        paper.abstract = (
            "We performed GWAS and Mendelian randomization to evaluate "
            "genetic determinants and causal effects for ischemic stroke."
        )
        paper.mesh_terms = ["Stroke", "Genetic Association Studies"]
        paper.keywords = ["GWAS", "Mendelian randomization"]

        ann = followup.annotate_paper(paper)

        self.assertIn("stroke", ann.disease_terms)
        self.assertIn("GWAS", ann.methods)
        self.assertIn("Mendelian randomization", ann.methods)
        self.assertIn("genomics", ann.data_types)
        self.assertEqual(ann.research_stage, "causal inference")
        self.assertIn("causal relationship", ann.research_question)


    def test_clinical_registry_biomarker_profile(self):
        paper = parsed_paper()
        paper.title = (
            "Evaluation of the Correlation of Calprotectin and SES-CD Score "
            "in Pediatric Crohn's Disease"
        )
        paper.abstract = (
            "This multicenter registry-based inception cohort evaluated serum "
            "calprotectin and endoscopic SES-CD correlation during treatment."
        )
        paper.mesh_terms = ["Crohn Disease", "Cohort Studies"]
        paper.keywords = ["calprotectin", "endoscopy", "registry"]

        ann = followup.annotate_paper(paper)

        self.assertIn("correlation analysis", ann.methods)
        self.assertIn("registry/inception cohort", ann.methods)
        self.assertIn("endoscopic assessment", ann.methods)
        self.assertIn("laboratory biomarker", ann.data_types)
        self.assertIn("endoscopy", ann.data_types)
        self.assertIn("clinical registry/cohort", ann.data_types)

    def test_progression_detects_method_and_data_shift(self):
        source = followup.PaperAnnotation(
            pmid="100",
            disease_terms=["stroke"],
            methods=["GWAS"],
            data_types=["genomics"],
            research_stage="discovery/association",
            research_question="Q1",
            evidence_terms=["gwas"],
        )
        target = followup.PaperAnnotation(
            pmid="200",
            disease_terms=["stroke"],
            methods=["single-cell RNA-seq"],
            data_types=["transcriptomics", "single-cell"],
            research_stage="mechanism",
            research_question="Q2",
            evidence_terms=["single-cell rna"],
        )

        summary = followup.progression_summary(source, target)

        self.assertIn("stage: discovery/association → mechanism", summary)
        self.assertIn("new method: single-cell RNA-seq", summary)
        self.assertIn("new data:", summary)

    def test_lineage_edge_contains_profile_progression(self):
        key = "orcid:0000-0001-2345-6789"
        registry = {
            key: followup.TrackedAuthor(
                key=key,
                name="Tae Hoon Kim",
                orcid="0000-0001-2345-6789",
                confidence="high",
                seed_pmids=["100"],
            )
        }
        paper_authors = {"100": {key}, "200": {key}}
        edges = [
            followup.ContinuationEdge(
                source="100",
                target="200",
                relation="author_continuation",
                tracked_author="Tae Hoon Kim",
                author_key=key,
                confidence="high",
                evidence="ORCID 0000-0001-2345-6789",
            )
        ]
        annotations = {
            "100": followup.PaperAnnotation(
                pmid="100",
                disease_terms=["stroke"],
                methods=["GWAS"],
                data_types=["genomics"],
                research_stage="discovery/association",
                research_question="Q1",
                evidence_terms=[],
            ),
            "200": followup.PaperAnnotation(
                pmid="200",
                disease_terms=["stroke"],
                methods=["single-cell RNA-seq"],
                data_types=["transcriptomics", "single-cell"],
                research_stage="mechanism",
                research_question="Q2",
                evidence_terms=[],
            ),
        }

        lineage_edges = followup.build_research_lineage_edges(
            edges, paper_authors, registry, annotations
        )

        self.assertEqual(lineage_edges[0].source_stage, "discovery/association")
        self.assertEqual(lineage_edges[0].target_stage, "mechanism")
        self.assertIn("single-cell RNA-seq", lineage_edges[0].progression)



class TestResearchAnnotation(unittest.TestCase):
    def test_gwas_stroke_profile(self):
        paper = parsed_paper()
        paper.pmid = "9001"
        paper.title = "Genome-wide association study of ischemic stroke"
        paper.abstract = "We identify genetic variants associated with ischemic stroke risk."
        paper.mesh_terms = ["Stroke", "Genome-Wide Association Study"]
        paper.keywords = ["GWAS", "genetic association"]

        ann = followup.annotate_paper(paper)
        self.assertIn("stroke", ann.disease_terms)
        self.assertIn("GWAS", ann.methods)
        self.assertIn("genomics", ann.data_types)
        self.assertEqual(ann.research_stage, "discovery/association")
        self.assertTrue(ann.research_question)

    def test_single_cell_rna_is_not_mislabeled_as_bulk(self):
        paper = parsed_paper()
        paper.pmid = "9002"
        paper.title = "Single-cell RNA-seq reveals immune cell states in cancer"
        paper.abstract = "Single-cell RNA-seq was used to characterize tumor immune populations."
        paper.mesh_terms = ["Neoplasms"]
        paper.keywords = ["single-cell RNA-seq"]

        ann = followup.annotate_paper(paper)
        self.assertIn("single-cell RNA-seq", ann.methods)
        self.assertNotIn("bulk RNA-seq", ann.methods)
        self.assertIn("single-cell", ann.data_types)
        self.assertIn("transcriptomics", ann.data_types)

    def test_progression_summary_detects_method_and_stage_change(self):
        p1 = parsed_paper()
        p1.pmid = "9003"
        p1.title = "Genome-wide association study of stroke"
        p1.abstract = "We identify variants associated with stroke."
        p1.mesh_terms = ["Stroke"]
        p1.keywords = ["GWAS"]

        p2 = parsed_paper()
        p2.pmid = "9004"
        p2.title = "Mendelian randomization of genetic risk factors for stroke"
        p2.abstract = "Mendelian randomization was used for causal inference."
        p2.mesh_terms = ["Stroke"]
        p2.keywords = ["Mendelian randomization"]

        a1 = followup.annotate_paper(p1)
        a2 = followup.annotate_paper(p2)
        progression = followup.progression_summary(a1, a2)

        self.assertIn("stage:", progression)
        self.assertIn("Mendelian randomization", progression)
        self.assertEqual(a2.research_stage, "causal inference")



class TestResearcherTrajectory(unittest.TestCase):
    def _followed(self, pmid, year, stage, methods, data_types):
        return followup.FollowedPaper(
            pmid=pmid,
            year=year,
            title=f"Paper {pmid}",
            journal="Journal",
            doi="",
            pubmed_url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            tracked_authors=["Tae Hoon Kim"],
            tracked_author_keys=["orcid:0000-0001-2345-6789"],
            current_affiliations=["Sungkyunkwan University"],
            skku_current=True,
            disease_terms=["stroke"],
            methods=methods,
            data_types=data_types,
            research_stage=stage,
            research_question="Example question?",
        )

    def test_researcher_trajectory_summary(self):
        tracked = followup.TrackedAuthor(
            key="orcid:0000-0001-2345-6789",
            name="Tae Hoon Kim",
            orcid="0000-0001-2345-6789",
            confidence="high",
            seed_pmids=["100"],
        )
        papers = [
            self._followed("100", 2020, "discovery/association", ["GWAS"], ["genomics"]),
            self._followed("200", 2024, "causal inference", ["Mendelian randomization"], ["genomics"]),
        ]
        lineage_edges = [
            followup.ResearchLineageEdge(
                source="100",
                target="200",
                relation="direct_citation_continuation",
                score=1.0,
                tracked_authors=["Tae Hoon Kim"],
                confidence="high",
                evidence="200 cites 100",
                source_stage="discovery/association",
                target_stage="causal inference",
                progression="stage: discovery/association → causal inference",
            )
        ]

        trajectories = followup.build_researcher_trajectories(
            [tracked], papers, lineage_edges
        )
        self.assertEqual(len(trajectories), 1)
        t = trajectories[0]
        self.assertEqual(t.paper_count, 2)
        self.assertEqual((t.first_year, t.last_year), (2020, 2024))
        self.assertEqual(
            t.stage_path, ["discovery/association", "causal inference"]
        )
        self.assertEqual(t.strong_lineage_count, 1)
        self.assertIn("GWAS", t.top_methods)
        self.assertIn("Mendelian randomization", t.top_methods)

    def test_trajectory_html_contains_researcher_and_graph_payload(self):
        tracked = followup.TrackedAuthor(
            key="orcid:0000-0001-2345-6789",
            name="Tae Hoon Kim",
            orcid="0000-0001-2345-6789",
            confidence="high",
            seed_pmids=["100"],
        )
        papers = [
            self._followed("100", 2020, "discovery/association", ["GWAS"], ["genomics"]),
            self._followed("200", 2024, "causal inference", ["Mendelian randomization"], ["genomics"]),
        ]
        lineage_edges = [
            followup.ResearchLineageEdge(
                source="100",
                target="200",
                relation="direct_citation_continuation",
                score=1.0,
                tracked_authors=["Tae Hoon Kim"],
                confidence="high",
                evidence="200 cites 100",
                source_stage="discovery/association",
                target_stage="causal inference",
                progression="stage: discovery/association → causal inference",
            )
        ]
        trajectories = followup.build_researcher_trajectories(
            [tracked], papers, lineage_edges
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "trajectory.html"
            followup.render_research_trajectory_html(
                path, papers, lineage_edges, trajectories
            )
            html = path.read_text(encoding="utf-8")

        self.assertIn("vis-network", html)
        self.assertIn("Tae Hoon Kim", html)
        self.assertIn("direct_citation_continuation", html)
        self.assertIn("strong citation lineage only", html)



class TestResearcherNetwork(unittest.TestCase):
    def test_collaboration_and_topic_overlap_are_combined(self):
        researchers = [
            {
                "key": "orcid:a",
                "name": "Researcher A",
                "orcid": "a",
                "paper_count": 4,
                "first_year": 2020,
                "last_year": 2025,
                "top_diseases": ["stroke"],
                "top_methods": ["GWAS", "Mendelian randomization"],
                "top_data_types": ["genomics"],
                "stage_path": ["discovery/association", "causal inference"],
                "strong_lineage_count": 1,
            },
            {
                "key": "orcid:b",
                "name": "Researcher B",
                "orcid": "b",
                "paper_count": 3,
                "first_year": 2021,
                "last_year": 2026,
                "top_diseases": ["stroke"],
                "top_methods": ["GWAS"],
                "top_data_types": ["genomics"],
                "stage_path": ["discovery/association"],
                "strong_lineage_count": 0,
            },
        ]
        papers = [
            {
                "pmid": "100",
                "tracked_author_keys": ["orcid:a", "orcid:b"],
                "tracked_authors": ["Researcher A", "Researcher B"],
            }
        ]
        lineage_edges = []

        nodes, edges = researcher_network.build_researcher_network(
            researchers, papers, lineage_edges, topic_threshold=0.20
        )

        self.assertEqual(len(nodes), 2)
        self.assertEqual(len(edges), 1)
        edge = edges[0]
        self.assertEqual(edge.relation, "collaboration")
        self.assertEqual(edge.shared_papers, 1)
        self.assertGreater(edge.topic_similarity, 0)
        self.assertIn("stroke", edge.shared_diseases)

    def test_cross_researcher_citation_is_high_confidence_network_edge(self):
        researchers = [
            {
                "key": "orcid:a",
                "name": "Researcher A",
                "orcid": "a",
                "paper_count": 2,
                "first_year": 2020,
                "last_year": 2024,
                "top_diseases": ["stroke"],
                "top_methods": ["GWAS"],
                "top_data_types": ["genomics"],
                "stage_path": ["discovery/association"],
                "strong_lineage_count": 1,
            },
            {
                "key": "orcid:b",
                "name": "Researcher B",
                "orcid": "b",
                "paper_count": 2,
                "first_year": 2021,
                "last_year": 2025,
                "top_diseases": ["cancer"],
                "top_methods": ["single-cell RNA-seq"],
                "top_data_types": ["transcriptomics"],
                "stage_path": ["mechanism"],
                "strong_lineage_count": 1,
            },
        ]
        papers = []
        lineage_edges = [
            {
                "source": "100",
                "target": "200",
                "relation": "cross_researcher_citation",
                "score": 0.85,
                "tracked_authors": ["Researcher A", "Researcher B"],
            }
        ]

        _, edges = researcher_network.build_researcher_network(
            researchers, papers, lineage_edges, topic_threshold=0.30
        )

        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0].relation, "citation")
        self.assertEqual(edges[0].direct_citations, 1)
        self.assertGreaterEqual(edges[0].score, 0.85)



class TestResearchCommunities(unittest.TestCase):
    def test_two_dense_groups_form_two_communities(self):
        nodes = [
            {
                "key": "a", "name": "A", "orcid": "a", "paper_count": 5,
                "first_year": 2020, "last_year": 2026,
                "top_diseases": ["stroke"], "top_methods": ["GWAS"],
                "top_data_types": ["genomics"],
                "stage_path": ["discovery/association"], "strong_lineage_count": 1,
            },
            {
                "key": "b", "name": "B", "orcid": "b", "paper_count": 4,
                "first_year": 2021, "last_year": 2026,
                "top_diseases": ["stroke"], "top_methods": ["Mendelian randomization"],
                "top_data_types": ["genomics"],
                "stage_path": ["causal inference"], "strong_lineage_count": 2,
            },
            {
                "key": "c", "name": "C", "orcid": "c", "paper_count": 6,
                "first_year": 2019, "last_year": 2025,
                "top_diseases": ["cancer"], "top_methods": ["single-cell RNA-seq"],
                "top_data_types": ["transcriptomics", "single-cell"],
                "stage_path": ["mechanism"], "strong_lineage_count": 1,
            },
            {
                "key": "d", "name": "D", "orcid": "d", "paper_count": 3,
                "first_year": 2022, "last_year": 2026,
                "top_diseases": ["cancer"], "top_methods": ["spatial transcriptomics"],
                "top_data_types": ["transcriptomics", "spatial"],
                "stage_path": ["mechanism"], "strong_lineage_count": 0,
            },
        ]
        edges = [
            {
                "source": "a", "target": "b", "relation": "collaboration",
                "score": 0.93, "shared_papers": 3, "direct_citations": 0,
                "topic_similarity": 0.6, "shared_diseases": ["stroke"],
                "shared_methods": [], "shared_data_types": ["genomics"], "evidence": "",
            },
            {
                "source": "c", "target": "d", "relation": "collaboration",
                "score": 0.91, "shared_papers": 2, "direct_citations": 0,
                "topic_similarity": 0.5, "shared_diseases": ["cancer"],
                "shared_methods": [], "shared_data_types": ["transcriptomics"], "evidence": "",
            },
            {
                "source": "b", "target": "c", "relation": "thematic_overlap",
                "score": 0.50, "shared_papers": 0, "direct_citations": 0,
                "topic_similarity": 0.11, "shared_diseases": [],
                "shared_methods": [], "shared_data_types": [], "evidence": "",
            },
        ]

        communities, members, community_edges = research_communities.build_research_communities(
            nodes, edges, min_edge_score=0.65
        )

        self.assertEqual(len(communities), 2)
        self.assertEqual(sorted(x.researcher_count for x in communities), [2, 2])
        self.assertEqual(len(members), 4)
        self.assertEqual(sum(m.role == "hub" for m in members), 2)
        labels = " | ".join(x.label for x in communities)
        self.assertIn("stroke", labels)
        self.assertIn("cancer", labels)
        self.assertEqual(len(community_edges), 1)

    def test_live_style_two_researcher_cluster_has_hub_and_theme(self):
        nodes = [
            {
                "key": "orcid:1", "name": "Researcher 1", "orcid": "1",
                "paper_count": 10, "first_year": 2020, "last_year": 2026,
                "top_diseases": ["Crohn Disease"],
                "top_methods": ["pharmacokinetics/TDM"],
                "top_data_types": ["clinical registry/cohort"],
                "stage_path": ["pharmacokinetics"], "strong_lineage_count": 1,
            },
            {
                "key": "orcid:2", "name": "Researcher 2", "orcid": "2",
                "paper_count": 8, "first_year": 2021, "last_year": 2026,
                "top_diseases": ["Crohn Disease"],
                "top_methods": ["biomarker analysis"],
                "top_data_types": ["clinical registry/cohort"],
                "stage_path": ["discovery/association"], "strong_lineage_count": 0,
            },
        ]
        edges = [
            {
                "source": "orcid:1", "target": "orcid:2", "relation": "collaboration",
                "score": 0.9357, "shared_papers": 5, "direct_citations": 0,
                "topic_similarity": 0.5714, "shared_diseases": ["Crohn Disease"],
                "shared_methods": [], "shared_data_types": ["clinical registry/cohort"],
                "evidence": "shared papers=5",
            }
        ]

        communities, members, _ = research_communities.build_research_communities(
            nodes, edges, min_edge_score=0.65
        )
        self.assertEqual(len(communities), 1)
        self.assertEqual(communities[0].researcher_count, 2)
        self.assertIn("Crohn Disease", communities[0].label)
        self.assertTrue(communities[0].hub_researcher)
        self.assertEqual(sum(m.role == "hub" for m in members), 1)


class FakePagedPubMedClient:
    def __init__(self, resolver):
        self.resolver = resolver
        self.calls = []

    def search_page(self, query, retstart=0, retmax=500):
        self.calls.append((query, retstart, retmax))
        ids = list(self.resolver(query))
        count = len(ids)
        if retmax == 0:
            return count, []
        return count, ids[retstart:retstart + retmax]


class TestPartitionedPubMedCrawl(unittest.TestCase):
    def test_year_partition_pagination_preserves_total_hit_count_with_cap(self):
        def resolver(query):
            if "2024:2025[dp]" in query:
                return [f"ALL{i}" for i in range(12)]
            if "2025:2025[dp]" in query:
                return [f"25{i}" for i in range(7)]
            if "2024:2024[dp]" in query:
                return [f"24{i}" for i in range(5)]
            return []

        client = FakePagedPubMedClient(resolver)
        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "crawl_checkpoint.json"
            total, pmids, windows = lineage.search_partitioned_pubmed(
                client,
                start_year=2024,
                end_year=2025,
                page_size=3,
                max_results=9,
                checkpoint_path=checkpoint,
            )

            self.assertEqual(total, 12)
            self.assertEqual(len(pmids), 9)
            self.assertEqual(pmids[:3], ["250", "251", "252"])
            self.assertEqual(len(windows), 2)
            self.assertTrue(checkpoint.exists())
            state = json.loads(checkpoint.read_text(encoding="utf-8"))
            self.assertEqual(len(state["retrieved_pmids"]), 9)
            self.assertEqual(
                [x["label"] for x in state["completed_windows"]],
                ["2025", "2024"],
            )

    def test_resume_skips_completed_year_windows(self):
        def resolver(query):
            if "2024:2025[dp]" in query:
                return [f"ALL{i}" for i in range(6)]
            if "2025:2025[dp]" in query:
                return ["a", "b", "c"]
            if "2024:2024[dp]" in query:
                return ["d", "e", "f"]
            return []

        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "crawl_checkpoint.json"
            first = FakePagedPubMedClient(resolver)
            total1, pmids1, _ = lineage.search_partitioned_pubmed(
                first,
                start_year=2024,
                end_year=2025,
                page_size=2,
                max_results=0,
                checkpoint_path=checkpoint,
            )
            second = FakePagedPubMedClient(resolver)
            total2, pmids2, windows2 = lineage.search_partitioned_pubmed(
                second,
                start_year=2024,
                end_year=2025,
                page_size=2,
                max_results=0,
                checkpoint_path=checkpoint,
                resume=True,
            )

            self.assertEqual((total1, total2), (6, 6))
            self.assertEqual(pmids1, pmids2)
            self.assertEqual(len(windows2), 2)
            # Resume still performs one overall count request but no year pages.
            self.assertEqual(len(second.calls), 1)
            self.assertIn("2024:2025[dp]", second.calls[0][0])

    def test_year_above_pubmed_10k_limit_is_split_by_month(self):
        def resolver(query):
            if "2026:2026[dp]" in query:
                return [f"Y{i}" for i in range(10001)]
            if "2026/12/01:2026/12/31[dp]" in query:
                return [f"D{i}" for i in range(6000)]
            return []

        client = FakePagedPubMedClient(resolver)
        total, pmids, windows = lineage.search_partitioned_pubmed(
            client,
            start_year=2026,
            end_year=2026,
            page_size=4,
            max_results=10,
        )

        self.assertEqual(total, 10001)
        self.assertEqual(len(pmids), 10)
        self.assertEqual(windows[0]["label"], "2026-12")
        self.assertTrue(
            any("2026/12/01:2026/12/31[dp]" in query for query, _, _ in client.calls)
        )

if __name__ == "__main__":
    unittest.main(verbosity=2)
