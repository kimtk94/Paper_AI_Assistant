import unittest
import datetime as dt
from feasibility import audit

BASE=dt.date(2023,1,1)
def row(i,override=None):
    r={'subject_id':f'SYNTH{i:04d}','irb_permits_reuse':'yes','diagnosis':'IPF',
    'baseline_plasma_available':'yes','plasma_date':BASE.isoformat(),
    'baseline_fvc_date':BASE.isoformat(), 'baseline_fvc_percent':'78',
    'followup_fvc_date':(BASE+dt.timedelta(days=365)).isoformat(),
    'followup_fvc_percent':'73','age':'68','sex':'M','smoking':'ever',
    'antifibrotic':'nintedanib','baseline_dlco_percent':'51','batch':'SYN1','hemolysis_flag':'no'}
    if override:r.update(override)
    return r
class Tests(unittest.TestCase):
    def test_80_go_exploratory_not_prognosis(self):
        x=audit([row(i) for i in range(80)])
        self.assertEqual(x['status'],'GO_EXPLORATORY_FEASIBILITY_ONLY')
        self.assertEqual(x['evaluable_n'],80)
    def test_60_pilot(self):
        self.assertEqual(audit([row(i) for i in range(60)])['status'],'GO_SMALL_EXPLORATORY_ONLY')
    def test_59_hold(self):
        self.assertTrue(audit([row(i) for i in range(59)])['status'].startswith('HOLD'))
    def test_duplicates_not_counted(self):
        x=audit([row(i) for i in range(80)]+[row(4)])
        self.assertEqual(x['evaluable_n'],80)
        self.assertEqual(x['exclusion_reasons']['missing_or_duplicate_subject_key'],1)
    def test_invalid_clinical_timing_and_irb(self):
        a=[row(i) for i in range(75)]+[row(90,{'followup_fvc_date':(BASE+dt.timedelta(days=600)).isoformat()}),row(91,{'irb_permits_reuse':'no'})]
        x=audit(a)
        self.assertEqual(x['evaluable_n'],75)
        self.assertEqual(x['exclusion_reasons']['followup_not_9_to_15_months'],1)
        self.assertEqual(x['exclusion_reasons']['IRB_or_consent_not_confirmed'],1)
    def test_no_patient_ids_in_output(self):
        x=audit([row(1)])
        self.assertNotIn('SYNTH0001',str(x))
if __name__=='__main__':unittest.main()