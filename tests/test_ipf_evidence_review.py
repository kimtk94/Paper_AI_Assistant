"""Synthetic tests for bounded archived-evidence replay; no patient data."""
import tempfile
import unittest
from pathlib import Path
from projects.ipf_evidence_review import replay


def row(donor, group, cpm=1, cells=20):
    return {'gene':'FCN3','celltype':'Endothelial Cells','donor':donor,
            'diagnosis':group,'n_cells':str(cells),'pseudobulk_cpm':str(cpm)}


class IPFReplayTests(unittest.TestCase):
    def test_loo_direction(self):
        r=replay.loo_effect([4,5,6,1,2,3],3)
        self.assertEqual(r['same_nonzero_direction'],6)
        self.assertEqual(r['log2fc'],-3)
    def test_loo_invalid_group(self):
        with self.assertRaises(ValueError):replay.loo_effect([1,2,3],2)
    def test_sparse_donors_suppress_effect(self):
        r=replay.donor_summary([row('a','IPF'),row('b','Control')],'FCN3','Endothelial Cells')
        self.assertEqual(r['status'],'HOLD_SPARSE_DONORS')
        self.assertIsNone(r['log2_median_ratio_plus1'])
    def test_source_match_and_cell_threshold(self):
        rows=[row(f'{g}{i}',g,1 if g=='IPF' else 3) for g in ('IPF','Control') for i in range(5)]
        rows.extend([row('other','IPF',99),row('few','Control',99,19)])
        sources={r['donor']:(r['diagnosis'],'Vanderbilt') for r in rows}
        sources['other']=('IPF','NTI')
        r=replay.donor_summary(rows,'FCN3','Endothelial Cells',sources=sources,site='Vanderbilt')
        self.assertEqual((r['ipf_donors'],r['control_donors']),(5,5))
        self.assertEqual(r['log2_median_ratio_plus1'],-1)
    def test_duplicate_donor_rejected(self):
        rows=[row('a','IPF'),row('a','Control')]
        with self.assertRaises(ValueError):replay.donor_summary(rows,'FCN3','Endothelial Cells')
    def test_source_diagnosis_conflict_rejected(self):
        with self.assertRaises(ValueError):replay.donor_summary([row('a','IPF')],'FCN3','Endothelial Cells',sources={'a':('Control','Vanderbilt')},site='Vanderbilt')
    def test_wrong_artifact_digest_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            aid,_=replay.ARTIFACTS['plasma']
            (Path(d)/f'github-actions-artifact-{aid}.zip').write_bytes(b'not the pinned artifact')
            with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                replay.member(d,'plasma','test')


if __name__=='__main__':unittest.main()
