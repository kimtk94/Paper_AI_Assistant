"""Synthetic regression tests; not research results."""
import csv, io, tempfile, unittest, zipfile
from pathlib import Path
from projects.endometriosis_gse283755.count_evidence import resolve_sample, classify_metric, validate_raw_counts
from projects.endometriosis_gse283755_validation.check_supplement import pdf_from_archive, MEMBER
from projects.endometriosis_gse283755_validation.compare_top5 import compare_results, TOP5

GROUPS={"GSM111":("untreated",1),"GSM222":("dienogest",1)}
LOOKUP={"GSM111":"12345","GSM222":"23456"}

class EvidenceTests(unittest.TestCase):
    def test_unique_sample_and_read_count(self):
        self.assertEqual(resolve_sample("12345_Read_Count",GROUPS,LOOKUP)[0],"GSM111")
        self.assertEqual(classify_metric("12345_Read_Count.1",GROUPS,LOOKUP),"count")
    def test_conflicting_ids(self):
        with self.assertRaises(ValueError):
            resolve_sample("GSM111_23456_Read_Count",GROUPS,LOOKUP)
    def test_normalized_and_expected_counts_not_raw(self):
        for name in ("12345_normalized_count","12345_RSEM_expected_count"):
            self.assertNotEqual(classify_metric(name,GROUPS,LOOKUP),"count")
    def test_bad_raw_counts(self):
        for value in (float("nan"),float("inf"),-1,100000.1):
            with self.assertRaises(ValueError):
                validate_raw_counts([[1,value],[1,2]])
        with self.assertRaises(ValueError):
            validate_raw_counts([[0,2],[0,2]])
        validate_raw_counts([[1,2],[0,2]])
    def test_archive_requires_unique_exact_member(self):
        payload=b"%PDF-1.4\n"+b"x"*30000
        for names in (["wrong.pdf"],[MEMBER,"other/"+MEMBER]):
            out=io.BytesIO()
            with zipfile.ZipFile(out,"w") as z:
                for n in names:z.writestr(n,payload)
            with self.assertRaises(ValueError):pdf_from_archive(out.getvalue())
        out=io.BytesIO()
        with zipfile.ZipFile(out,"w") as z:z.writestr(MEMBER,payload)
        self.assertEqual(pdf_from_archive(out.getvalue()),payload)
    def test_top5_never_authorizes_full_reproduction(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"de.csv"
            with f.open("w",newline="") as o:
                w=csv.writer(o);w.writerow(["gene","log2FoldChange","padj"])
                for gene,lfc in TOP5.items():w.writerow([gene,lfc,"0.001"])
            r=compare_results(f)
            self.assertEqual(r["matched"],5)
            self.assertFalse(r["complete_406_gene_reproduction"])
            f.write_text("gene,log2FoldChange,padj\nHLA-DOB,1.81,0.001\n")
            self.assertEqual(compare_results(f)["matched"],0)
