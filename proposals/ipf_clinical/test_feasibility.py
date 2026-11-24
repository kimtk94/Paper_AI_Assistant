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
        self.assertEqual(x['evaluable_n'],79)
        self.assertEqual(x['duplicate_subject_rows_excluded'],2)
        self.assertEqual(x['exclusion_reasons']['missing_or_duplicate_subject_key'],2)
    def test_invalid_clinical_timing_and_irb(self):
        a=[row(i) for i in range(75)]+[row(90,{'followup_fvc_date':(BASE+dt.timedelta(days=600)).isoformat()}),row(91,{'irb_permits_reuse':'no'})]
        x=audit(a)
        self.assertEqual(x['evaluable_n'],75)
        self.assertEqual(x['exclusion_reasons']['followup_not_9_to_15_months'],1)
        self.assertEqual(x['exclusion_reasons']['IRB_or_consent_not_confirmed'],1)
    def test_unknown_hemolysis_cannot_qualify_cohort(self):
        for value in ('', None, 'unknown', 'NA', 'unreviewed', 'false-ish', '2'):
            with self.subTest(hemolysis=value):
                x=audit([row(i,{'hemolysis_flag':value}) for i in range(80)])
                self.assertTrue(x['status'].startswith('HOLD'))
                self.assertEqual(x['evaluable_n'],0)
                self.assertEqual(x['exclusion_reasons']['hemolysis_not_confirmed'],80)
    def test_explicit_negative_hemolysis_is_accepted(self):
        for value in ('no', 'n', 'false', '0', ' NO ', 'False', False, 0):
            with self.subTest(hemolysis=value):
                x=audit([row(i,{'hemolysis_flag':value}) for i in range(80)])
                self.assertEqual(x['status'],'GO_EXPLORATORY_FEASIBILITY_ONLY')
                self.assertEqual(x['evaluable_n'],80)
    def test_positive_hemolysis_requires_review(self):
        for value in ('yes', 'y', 'true', '1', ' YES ', True, 1):
            with self.subTest(hemolysis=value):
                x=audit([row(1,{'hemolysis_flag':value})])
                self.assertEqual(x['evaluable_n'],0)
                self.assertEqual(x['exclusion_reasons']['hemolysis_review_required'],1)
    def test_age_must_be_finite_numeric(self):
        for value in ('sixty-eight', 'inf', '-inf', '1e309', float('inf')):
            with self.subTest(age=value):
                x=audit([row(i,{'age':value}) for i in range(80)])
                self.assertTrue(x['status'].startswith('HOLD'))
                self.assertEqual(x['evaluable_n'],0)
                self.assertEqual(x['exclusion_reasons']['age_not_finite_numeric'],80)
    def test_missing_covariate_tokens_cannot_qualify_cohort(self):
        for field in ('age','sex','smoking','antifibrotic','batch'):
            for value in ('unknown',' N/A ','nan',float('nan'),None,' ','not recorded'):
                with self.subTest(field=field,value=value):
                    x=audit([row(i,{field:value}) for i in range(80)])
                    self.assertTrue(x['status'].startswith('HOLD'))
                    self.assertEqual(x['evaluable_n'],0)
                    self.assertEqual(x['exclusion_reasons']['missing_critical_covariate'],80)
    def test_explicit_absence_of_smoking_or_treatment_is_recorded(self):
        for smoking in ('never','none','no'):
            for treatment in ('none','no','nintedanib'):
                with self.subTest(smoking=smoking,treatment=treatment):
                    x=audit([row(i,{'age':'68.5','smoking':smoking,'antifibrotic':treatment}) for i in range(80)])
                    self.assertEqual(x['status'],'GO_EXPLORATORY_FEASIBILITY_ONLY')
                    self.assertEqual(x['evaluable_n'],80)
    def test_unknown_qc_cannot_cross_minimum_pilot_size(self):
        x=audit([row(i) for i in range(59)]+[row(59,{'hemolysis_flag':'unknown'})])
        self.assertEqual(x['status'],'HOLD_FEWER_THAN_60_EVALUABLE')
        self.assertEqual(x['evaluable_n'],59)
    def test_no_patient_ids_in_output(self):
        x=audit([row(1)])
        self.assertNotIn('SYNTH0001',str(x))
if __name__=='__main__':unittest.main()
