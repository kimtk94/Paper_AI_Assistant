"""Offline fail-closed checks for official S2 patient/ROI provenance."""
import unittest

from projects.hs_spatial_gse294009.hs_s1s2_site_stage_qc import verify_explicit_s2_mapping


class SourceS2Tests(unittest.TestCase):
    def setUp(self):
        self.found = {
            "GSM1111111": [{"sheet": "xl/worksheets/sheet2.xml", "row": 3, "same_row_patient": True}],
            "GSM2222222": [{"sheet": "xl/worksheets/sheet2.xml", "row": 4, "same_row_patient": True}],
        }
        self.gsms = {"GSM1111111", "GSM2222222"}

    def test_unique_explicit_source_rows_pass(self):
        self.assertEqual(verify_explicit_s2_mapping(self.found, self.gsms), 2)

    def test_missing_gsm_fails(self):
        self.found.pop("GSM2222222")
        with self.assertRaisesRegex(ValueError, "not independently found"):
            verify_explicit_s2_mapping(self.found, self.gsms)

    def test_duplicate_gsm_rows_fail(self):
        self.found["GSM1111111"].append({"sheet": "xl/worksheets/sheet2.xml", "row": 6, "same_row_patient": True})
        with self.assertRaisesRegex(ValueError, "duplicated"):
            verify_explicit_s2_mapping(self.found, self.gsms)

    def test_wrong_patient_same_row_fails(self):
        self.found["GSM1111111"][0]["same_row_patient"] = False
        with self.assertRaisesRegex(ValueError, "not explicit"):
            verify_explicit_s2_mapping(self.found, self.gsms)

    def test_wrong_source_sheet_fails(self):
        self.found["GSM1111111"][0]["sheet"] = "xl/worksheets/sheet3.xml"
        with self.assertRaisesRegex(ValueError, "outside"):
            verify_explicit_s2_mapping(self.found, self.gsms)


if __name__ == "__main__":
    unittest.main()
