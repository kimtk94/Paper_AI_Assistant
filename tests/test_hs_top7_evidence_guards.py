"""Synthetic-only checks; no actual patient expression values."""
import gzip
import io
import tarfile
import unittest
from projects.top7_remaining4.evidence_guards import verify_archive_mapping, paired_lesions, nmosd_track_mapping


def fixture(v2=7, member_name="GSM2_19-508-008.txt.gz", header="X19.508.008"):
    out=io.BytesIO()
    files={"GSM1_GTG-360-B1-47.txt.gz":b"GTG.360.B1.47\n100 1\n200 2\n",
           member_name:(f"{header}\n100 3\n200 {v2}\n").encode()}
    with tarfile.open(fileobj=out,mode="w") as t:
        for name,data in files.items():
            payload=gzip.compress(data,mtime=0)
            info=tarfile.TarInfo(name);info.size=len(payload)
            t.addfile(info,io.BytesIO(payload))
    return out.getvalue()


def verify(blob):
    vals={"GTG.360.B1.47":[1,2],"X19.508.008":[3,7]}
    return verify_archive_mapping(blob,list(vals),["100","200"],lambda c:vals[c],["GSM1","GSM2"])


def pairs():
    return [{"gsm":f"GSM{i+1}","column":f"c{i}","patient":p,"tissue":"lesional","visit":v,"responder":r}
            for i,(p,r,v) in enumerate([("a","YES","V1"),("a","YES","V12"),("b","NO","V1"),("b","NO","V12")])]


class Top7Guards(unittest.TestCase):
    def test_full_count_vectors(self):
        self.assertEqual(verify(fixture()),{"GTG.360.B1.47":"GSM1","X19.508.008":"GSM2"})
    def test_count_difference_rejected(self):
        with self.assertRaisesRegex(ValueError,"full source count vector"):verify(fixture(v2=8))
    def test_header_conflict_rejected(self):
        with self.assertRaisesRegex(ValueError,"header conflicts"):verify(fixture(header="other"))
    def test_unknown_filename_rejected(self):
        with self.assertRaisesRegex(ValueError,"ambiguous or unknown"):verify(fixture(member_name="GSM2_19-508-009.txt.gz"))
    def test_duplicate_gsm_rejected(self):
        with self.assertRaises(ValueError):verify(fixture(member_name="GSM1_19-508-008.txt.gz"))
    def test_pairs_and_response(self):
        all_pairs,lesion=paired_lesions(pairs(),expected_patients=2,expected_each=1)
        self.assertEqual(len(all_pairs),2);self.assertEqual(len(lesion),2)
    def test_duplicate_visit_rejected(self):
        r=pairs();r[1]["visit"]="V1"
        with self.assertRaisesRegex(ValueError,"duplicate patient"):paired_lesions(r,2,1)
    def test_response_mismatch_rejected(self):
        r=pairs();r[1]["responder"]="NO"
        with self.assertRaisesRegex(ValueError,"inconsistent"):paired_lesions(r,2,1)
    def test_unknown_response_rejected(self):
        r=pairs();r[1]["responder"]="UNKNOWN"
        with self.assertRaises(ValueError):paired_lesions(r,2,1)
    def test_missing_track_is_explicit(self):
        s=[{"accession":"GSM1","title":"T-AQP4_01 [L1_Track-1]"},
           {"accession":"GSM2","title":"T-AQP4_01 [L2_Track-2]"}]
        mapped,missing=nmosd_track_mapping(s,["","Gene_id","gene_biotype","description","L1_Track-1"])
        self.assertEqual(len(mapped),1);self.assertEqual(missing,["L2_Track-2"])
    def test_duplicate_or_unknown_tracks_rejected(self):
        s=[{"accession":"GSM1","title":"T-AQP4_01 [L1_Track-1]"}]
        h=["","Gene_id","gene_biotype","description"]
        for cols in (["L1_Track-1","L1_Track-1"],["L9_Track-9"]):
            with self.assertRaises(ValueError):nmosd_track_mapping(s,h+cols)
        with self.assertRaises(ValueError):
            nmosd_track_mapping(s+[{"accession":"GSM2","title":"T-HD_01 [L1_Track-1]"}],h+["L1_Track-1"])


if __name__=="__main__":unittest.main()
