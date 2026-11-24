"""Synthetic-only fixtures: no asserted or fabricated study expression values."""
import csv
import tempfile
import unittest
from pathlib import Path

from projects.hs_spatial_gse294009.historic24_quant_candidate_gate import (
    evaluate, read_expected_manifest, scan_candidate,
)


class Historic24CandidateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = self.root/'historic.csv'
        self.keys = [f'{donor}:{i}' for donor, n in [('HS41',6),('HS42',6),('HS43',12)] for i in range(1,n+1)]
        with self.manifest.open('w',newline='') as f:
            w = csv.DictWriter(f,fieldnames=['ROI_composite_key','patient_id','author_ROI_label',
                                             'source_paper_doi','expression_value_available'])
            w.writeheader()
            for key in self.keys:
                donor, roi = key.split(':')
                w.writerow({'ROI_composite_key':key,'patient_id':donor,'author_ROI_label':roi,
                            'source_paper_doi':'10.1093/bjd/ljad274',
                            'expression_value_available':'NO'})
        self.matrix = self.root/'candidate.csv'
        self.write_counts()

    def write_counts(self, override=None, heading=None, nrows=2):
        with self.matrix.open('w',newline='') as f:
            w=csv.writer(f)
            w.writerow(['gene_id',*(heading or self.keys)])
            for i in range(nrows):
                values=['1']*len(self.keys)
                if override is not None: values[0]=override
                w.writerow([f'gene_{i}',*values])

    def test_manifest_only_remains_hold(self):
        v=evaluate(self.manifest)
        self.assertIn('HOLD',v['decision'])
        self.assertFalse(v['combined_156_roi_analysis_allowed'])

    def test_valid_synthetic_counts_do_not_prove_source(self):
        self.write_counts(nrows=1000)
        v=scan_candidate(self.matrix,read_expected_manifest(self.manifest))
        self.assertEqual(v['n_roi_columns'],24)
        self.assertEqual(len(v['candidate_input_sha256']),64)
        evaluated=evaluate(self.manifest,self.matrix)
        self.assertEqual(evaluated['decision'],'CANDIDATE_FORMAT_PASS__HOLD_INDEPENDENT_SOURCE_AND_METHOD_REVIEW')
        self.assertFalse(evaluated['source_authenticity_verified'])
        self.assertFalse(evaluated['combined_156_roi_analysis_allowed'])

    def test_wrong_roi_rejected(self):
        names=list(self.keys)
        names[0]='HS44:1'
        self.write_counts(heading=names)
        with self.assertRaisesRegex(ValueError,'24 official HS ROI'):
            scan_candidate(self.matrix,read_expected_manifest(self.manifest),minimum_genes=2)

    def test_duplicate_column_rejected(self):
        names=list(self.keys)
        names[0]=names[1]
        self.write_counts(heading=names)
        with self.assertRaisesRegex(ValueError,'24 official HS ROI'):
            scan_candidate(self.matrix,read_expected_manifest(self.manifest),minimum_genes=2)

    def test_normalized_fraction_rejected(self):
        self.write_counts(override='0.75')
        with self.assertRaisesRegex(ValueError,'noninteger'):
            scan_candidate(self.matrix,read_expected_manifest(self.manifest),minimum_genes=2)

    def test_negative_rejected(self):
        self.write_counts(override='-5')
        with self.assertRaisesRegex(ValueError,'noninteger'):
            scan_candidate(self.matrix,read_expected_manifest(self.manifest),minimum_genes=2)

    def test_too_few_genes_rejected(self):
        with self.assertRaisesRegex(ValueError,'too small'):
            scan_candidate(self.matrix,read_expected_manifest(self.manifest))

    def test_zero_library_rejected(self):
        with self.matrix.open('w',newline='') as f:
            w=csv.writer(f)
            w.writerow(['gene_id',*self.keys])
            for i in range(2):w.writerow([f'g{i}',0,*([1]*23)])
        with self.assertRaisesRegex(ValueError,'zero counts'):
            scan_candidate(self.matrix,read_expected_manifest(self.manifest),minimum_genes=2)

    def test_wrong_manifest_counts_rejected(self):
        bad=self.root/'bad.csv'
        bad.write_text(self.manifest.read_text().splitlines()[0]+'\n'+'\n'.join(self.manifest.read_text().splitlines()[1:24])+'\n')
        with self.assertRaisesRegex(ValueError,'manifest ROI/donor count'):
            read_expected_manifest(bad)


if __name__=='__main__':unittest.main()
