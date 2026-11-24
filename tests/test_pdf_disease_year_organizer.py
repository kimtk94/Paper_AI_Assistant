import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pdf_disease_year_organizer as organizer


class FilenameYearTests(unittest.TestCase):
    def test_pmc_identifier_does_not_override_publication_year(self):
        self.assertEqual(organizer.infer_year(Path("PMC2012345.pdf"), "Published in 2024"), 2024)

    def test_identifiers_without_publication_evidence_have_unknown_year(self):
        for filename in ("PMC2012345.pdf", "PMC2012.pdf", "PMID20241234.pdf", "20241234.pdf", "id2024x.pdf"):
            with self.subTest(filename=filename):
                self.assertIsNone(organizer.infer_year(Path(filename), "No publication date available"))

    def test_standalone_filename_year_still_has_priority(self):
        for filename in ("study_2022_final.pdf", "study-2022.pdf", "2022-study.pdf", "study (2022).pdf"):
            with self.subTest(filename=filename):
                self.assertEqual(organizer.infer_year(Path(filename), "Reference published in 2024"), 2022)


if __name__ == "__main__":
    unittest.main()
