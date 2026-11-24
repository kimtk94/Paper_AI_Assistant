import sys
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
from urllib.error import URLError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import paper_assistants as papers


class PdfDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.outdir = Path(self.temp.name)
        self.downloader = papers.PdfDownloader(outdir=self.outdir, unpaywall_email="test@example.com")
        self.destination = self.outdir / "10.1234_example.pdf"
        self.original = b"%PDF-1.4\noriginal content\n%%EOF"
        self.destination.write_bytes(self.original)

    def test_html_oa_landing_page_is_reported_missing_and_preserves_pdf(self):
        with (
            patch.object(self.downloader, "_id_convert", return_value={"doi": "10.1234/example"}),
            patch.object(
                self.downloader, "_resolve_unpaywall_pdf_url", return_value="https://example.invalid/landing"
            ),
            patch.object(papers, "urlopen", return_value=BytesIO(b"<html>Landing page</html>")),
        ):
            result = self.downloader.download("10.1234/example")
        self.assertFalse(result.downloaded)
        self.assertIsNone(result.filename)
        self.assertEqual(self.destination.read_bytes(), self.original)

    def test_non_pdf_responses_do_not_create_files(self):
        target = self.outdir / "missing.pdf"
        for response in (b"", b"{\"error\": \"blocked\"}", b"<html>Paywall</html>"):
            with self.subTest(response=response):
                with patch.object(papers, "urlopen", return_value=BytesIO(response)):
                    self.assertFalse(self.downloader._download_file("https://example.invalid/pdf", target))
                self.assertFalse(target.exists())

    def test_network_failure_preserves_existing_pdf(self):
        with patch.object(papers, "urlopen", side_effect=URLError("unavailable")):
            self.assertFalse(
                self.downloader._download_file("https://example.invalid/pdf", self.destination)
            )
        self.assertEqual(self.destination.read_bytes(), self.original)

    def test_valid_pdf_replaces_target_without_leaving_temporary_file(self):
        updated = b"%PDF-1.7\nupdated content\n%%EOF"
        with patch.object(papers, "urlopen", return_value=BytesIO(updated)):
            self.assertTrue(
                self.downloader._download_file("https://example.invalid/pdf", self.destination)
            )
        self.assertEqual(self.destination.read_bytes(), updated)
        self.assertEqual(list(self.outdir.iterdir()), [self.destination])

    def test_replacement_failure_preserves_existing_pdf_and_cleans_up(self):
        updated = b"%PDF-1.7\nupdated content\n%%EOF"
        with (
            patch.object(papers, "urlopen", return_value=BytesIO(updated)),
            patch.object(Path, "replace", side_effect=OSError("replacement failed")),
        ):
            with self.assertRaises(OSError):
                self.downloader._download_file("https://example.invalid/pdf", self.destination)
        self.assertEqual(self.destination.read_bytes(), self.original)
        self.assertEqual(list(self.outdir.iterdir()), [self.destination])


class MissingPdfCandidateTests(unittest.TestCase):
    def test_doi_url_and_case_variants_match_without_pmid(self):
        for paper_doi, result_doi in (
            ("https://doi.org/10.1234/AbC", "10.1234/abc"),
            ("10.1234/abc", "https://doi.org/10.1234/ABC"),
            ("https://dx.doi.org/10.1234/ABC", "10.1234/AbC"),
        ):
            with self.subTest(paper_doi=paper_doi, result_doi=result_doi):
                paper = papers.RetrievalAssistant(show_progress=False)._parse_work(
                    {"id": "https://openalex.org/W1", "title": "Example", "doi": paper_doi}
                )
                candidates = papers._collect_missing_pdf_candidates(
                    [paper],
                    pdf_payload={
                        "results": [{"downloaded": False, "doi": result_doi, "pmid": None, "reason": "No OA"}]
                    },
                )
                self.assertEqual(len(candidates), 1)
                self.assertEqual(candidates[0].paper_id, paper.paper_id)
                self.assertEqual(candidates[0].identifier, "10.1234/abc")
                self.assertEqual(candidates[0].doi, paper_doi)
                self.assertIsNone(candidates[0].pmid)

    def test_unrelated_doi_is_not_assigned_to_a_paper(self):
        paper = papers.RetrievalAssistant(show_progress=False)._parse_work(
            {"id": "https://openalex.org/W1", "doi": "https://doi.org/10.1234/abc"}
        )
        candidates = papers._collect_missing_pdf_candidates(
            [paper], pdf_payload={"results": [{"downloaded": False, "doi": "10.1234/def"}]}
        )
        self.assertEqual(candidates, [])


if __name__ == "__main__":
    unittest.main()
