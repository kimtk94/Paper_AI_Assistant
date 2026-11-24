"""Network-free regressions for incomplete PubMed responses and resume state."""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from skku_pubmed_lineage import PubMedClient, search_partitioned_pubmed


def article_xml(*pmids):
    articles = "".join(
        f"<PubmedArticle><MedlineCitation><PMID>{pmid}</PMID>"
        "<Article><ArticleTitle>Synthetic record</ArticleTitle></Article>"
        "</MedlineCitation></PubmedArticle>"
        for pmid in pmids
    )
    return f"<PubmedArticleSet>{articles}</PubmedArticleSet>".encode()


class TestPubMedResponseValidation(unittest.TestCase):
    def setUp(self):
        self.client = PubMedClient("synthetic@example.invalid")

    def test_search_rejects_error_or_malformed_response(self):
        for response in (
            {"error": "upstream failure"},
            {},
            [],
            {"esearchresult": {"error": "query failed"}},
            {"esearchresult": {"idlist": []}},
            {"esearchresult": {"count": "bad", "idlist": []}},
            {"esearchresult": {"count": "0", "idlist": ""}},
        ):
            with self.subTest(response=response):
                self.client._request = Mock(return_value=json.dumps(response).encode())
                with self.assertRaises(RuntimeError):
                    self.client.search("synthetic", 2)
        self.client._request = Mock(return_value=b"<html>unavailable</html>")
        with self.assertRaises(RuntimeError):
            self.client.search("synthetic", 2)

    def test_search_rejects_incomplete_and_duplicate_pages(self):
        for ids in (["1"], ["1", "1"]):
            with self.subTest(ids=ids):
                self.client._request = Mock(return_value=json.dumps(
                    {"esearchresult": {"count": "3", "idlist": ids}}
                ).encode())
                with self.assertRaisesRegex(RuntimeError, "incomplete or duplicate"):
                    self.client.search("synthetic", 2)

    def test_search_accepts_true_empty_results_and_count_only_queries(self):
        for count in (0, 19):
            with self.subTest(count=count):
                self.client._request = Mock(return_value=json.dumps(
                    {"esearchresult": {"count": str(count), "idlist": []}}
                ).encode())
                self.assertEqual(self.client.search_page("synthetic", retmax=0), (count, []))

    def test_fetch_rejects_invalid_or_incomplete_xml_before_caching(self):
        for payload in (
            b"<broken",
            b"<html>unavailable</html>",
            b"<eFetchResult><ERROR>upstream failure</ERROR></eFetchResult>",
            b"<PubmedArticleSet><ERROR>upstream failure</ERROR></PubmedArticleSet>",
            article_xml(),
            article_xml("2"),
            article_xml("1", "1"),
        ):
            with self.subTest(payload=payload), tempfile.TemporaryDirectory() as tmp:
                cache = Path(tmp)
                self.client._request = Mock(return_value=payload)
                with self.assertRaises(RuntimeError):
                    self.client.fetch(["1"], cache_dir=cache)
                self.assertEqual(list(cache.iterdir()), [])

    def test_fetch_uses_valid_cache_without_network(self):
        self.client._request = Mock(return_value=article_xml("1"))
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp)
            self.client.fetch(["1"], cache_dir=cache)
            self.client._request.reset_mock()
            papers = self.client.fetch(["1"], cache_dir=cache, resume=True)
            self.assertEqual([p.pmid for p in papers], ["1"])
            self.client._request.assert_not_called()

    def test_parseable_invalid_cache_is_replaced_by_verified_download(self):
        self.client._request = Mock(return_value=article_xml("1"))
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp)
            self.client.fetch(["1"], cache_dir=cache)
            saved = next(cache.glob("*.xml"))
            for bad_payload in (
                b"<eFetchResult><ERROR>upstream failure</ERROR></eFetchResult>",
                article_xml("2"),
            ):
                with self.subTest(payload=bad_payload):
                    saved.write_bytes(bad_payload)
                    self.client._request.reset_mock()
                    papers = self.client.fetch(["1"], cache_dir=cache, resume=True)
                    self.assertEqual([p.pmid for p in papers], ["1"])
                    self.client._request.assert_called_once()
                    self.assertEqual(saved.read_bytes(), article_xml("1"))

    def test_fetch_validates_batch_size(self):
        for size in (0, -1):
            with self.subTest(size=size), self.assertRaises(ValueError):
                self.client.fetch(["1"], batch_size=size)


class TestIncompleteCrawlCheckpoints(unittest.TestCase):
    def test_legacy_checkpoint_cannot_hide_a_previously_truncated_window(self):
        client = Mock()
        client.search_page.side_effect = lambda query, retstart=0, retmax=500: (
            3, ["1", "2", "3"][retstart:retstart + retmax]
        )
        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "checkpoint.json"
            checkpoint.write_text(json.dumps({
                "version": 1,
                "spec": {"start_year": 2025, "end_year": 2025,
                         "topic": "", "extra": "", "page_size": 2, "max_results": 0},
                "retrieved_pmids": ["1"],
                "completed_windows": [{"label": "2025", "pubmed_hits": 3,
                                       "retrieved_target": 3, "retrieved_total_after_window": 1}],
            }))
            _, ids, _ = search_partitioned_pubmed(
                client, start_year=2025, end_year=2025, page_size=2,
                checkpoint_path=checkpoint, resume=True,
            )
            self.assertEqual(ids, ["1", "2", "3"])
            self.assertEqual(json.loads(checkpoint.read_text())["version"], 2)

    def test_truncated_window_preserves_only_previously_completed_progress(self):
        class PartialClient:
            def search_page(self, query, retstart=0, retmax=500):
                if "2024:2025" in query:
                    return 4, []
                if "2025:2025" in query:
                    return 1, ["done"][:retmax]
                if retmax == 0:
                    return 3, []
                return 3, ["partial"] if retstart == 0 else []

        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "checkpoint.json"
            with self.assertRaisesRegex(RuntimeError, "incomplete window"):
                search_partitioned_pubmed(
                    PartialClient(), start_year=2024, end_year=2025,
                    page_size=2, checkpoint_path=checkpoint,
                )
            state = json.loads(checkpoint.read_text())
            self.assertEqual(state["retrieved_pmids"], ["done"])
            self.assertEqual([w["label"] for w in state["completed_windows"]], ["2025"])
            self.assertFalse(checkpoint.with_suffix(".json.tmp").exists())

    def test_duplicate_page_never_completes_window(self):
        class DuplicateClient:
            def search_page(self, query, retstart=0, retmax=500):
                return (3, []) if retmax == 0 else (3, ["1"])

        with tempfile.TemporaryDirectory() as tmp:
            checkpoint = Path(tmp) / "checkpoint.json"
            with self.assertRaisesRegex(RuntimeError, "duplicate"):
                search_partitioned_pubmed(
                    DuplicateClient(), start_year=2025, end_year=2025,
                    page_size=1, checkpoint_path=checkpoint,
                )
            self.assertFalse(checkpoint.exists())

    def test_overlapping_years_do_not_undershoot_unique_result_cap(self):
        class OverlappingClient:
            def search_page(self, query, retstart=0, retmax=500):
                ids = ["1"] if "2025:2025" in query else ["1", "2", "3"]
                return len(ids), ids[retstart:retstart + retmax]

        _, pmids, windows = search_partitioned_pubmed(
            OverlappingClient(), start_year=2024, end_year=2025,
            page_size=2, max_results=3,
        )
        self.assertEqual(pmids, ["1", "2", "3"])
        self.assertEqual(windows[-1]["retrieved_target"], 3)


if __name__ == "__main__":
    unittest.main()
