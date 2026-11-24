"""Synthetic fixtures for source-only feasibility checks."""
import io
import unittest
import zipfile
from projects.migraine_mtbls13630 import direct_source_audit as q


class DirectSourceTests(unittest.TestCase):
    def test_source_names_are_not_donors(self):
        data=b"Sample Name\tSource Name\tCharacteristics[Sample type]\tCharacteristics[Variant]\nS_A\tA\tVenous Blood\t\nS_B\tB\tarterial blood\t\n"
        ids,summary=q.sample_scope(data)
        self.assertEqual(ids,{"S_A","S_B"})
        self.assertFalse(summary["independent_participant_mapping_verified"])
        self.assertEqual(summary["blank_variant_rows"],2)

    def test_duplicate_sample_rejected(self):
        with self.assertRaises(ValueError):
            q.sample_scope(b"Sample Name\tSource Name\tCharacteristics[Sample type]\tCharacteristics[Variant]\nS_A\tA\tvenous\t\nS_A\tB\tarterial\t\n")

    def test_assay_link_must_be_exact(self):
        with self.assertRaises(ValueError):
            q.assay_scope(b"Sample Name\tMetabolite Assignment File\nS_WRONG\tm.tsv\n",{"S_A"})

    def test_aggregate_annotation_not_quant(self):
        r=q.maf_scope(b"database_identifier\tsmallmolecule_abundance_sub\nHMDB1\t123.5\n",{"S_A"},{"S_A"})
        self.assertEqual(r["finite_numeric_sample_cells"],0)
        self.assertEqual(r["status"],"NO_EXACT_LINKED_NUMERIC_VALUES")

    def test_exact_columns_only(self):
        r=q.maf_scope(b"database_identifier\tS_A\tS_B\nHMDB1\t0\t3.4\nHMDB2\tNA\t4\n",{"S_A","S_B"},{"S_A","S_B"})
        self.assertEqual(r["finite_numeric_sample_cells"],3)
        self.assertEqual(r["missing_sample_cells"],1)
        self.assertIn("NEED_UNITS",r["status"])

    def test_duplicate_and_wrong_assay_columns_rejected(self):
        for h in ("S_A\tS_A","S_A\tS_B"):
            with self.assertRaises(ValueError):
                q.maf_scope(("database_identifier\t"+h+"\nM\t1\t2\n").encode(),{"S_A","S_B"},{"S_A"})

    def test_nonfinite_or_nonnumeric_held(self):
        r=q.maf_scope(b"database_identifier\tS_A\nM\tinf\nN\ttext\n",{"S_A"},{"S_A"})
        self.assertEqual(r["invalid_sample_cells"],2)
        self.assertIn("HOLD",r["status"])

    def test_ragged_or_html_not_table(self):
        for b in (b"a\tb\nx\n",b"<html>error</html>"):
            with self.assertRaises(ValueError):q.tsv(b)

    def test_directory_only_bounded_inventory(self):
        r=q.directory_scope(b'<a href="../">up</a><a href="FILE.raw.zip">raw</a><a href="quant.csv">q</a><a href="sub/">sub</a><a href="https://x/a">x</a>')
        self.assertEqual(r["immediate_files"],2)
        self.assertEqual(r["potential_processed_file_count"],1)
        self.assertEqual(r["immediate_subdirectories"],1)

    def test_supplement_inventory_not_contents(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:z.writestr("supplement.xlsx",b"synthetic")
        r=q.supplement_scope(b.getvalue())
        self.assertEqual(r["tabular_member_count"],1)
        self.assertFalse(r["contents_verified"])


if __name__=="__main__":unittest.main()
