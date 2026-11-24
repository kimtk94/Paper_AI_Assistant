"""Synthetic-only v2 clinical screening tests; no actual patient information."""
import tempfile
import unittest
from pathlib import Path
from test_feasibility import row
from feasibility import audit,load_manifest
from study_readiness import screen,planning,sample_size_for_correlation


def v2(i, overrides=None):
    r=row(i,{'diagnosis_adjudicated':'yes','same_person_linkage_verified':'yes',
             'specimen_matrix':'plasma','anticoagulant':'EDTA','preanalytics_reviewed':'yes',
             'baseline_treatment_timing_verified':'yes','endpoint_status':'OBSERVED_NATIVE_LUNG_FVC',
             'endpoint_adjudicated':'yes','native_lung_fvc_confirmed':'yes',
             'fvc_quality_reviewed':'yes','terminal_event_date':''})
    if overrides:r.update(overrides)
    return r


class ReadinessTests(unittest.TestCase):
    def test_80_complete_not_authorized(self):
        a=screen([v2(i) for i in range(80)])
        self.assertEqual(a['evaluable_observed_native_lung_fvc_n'],80)
        self.assertTrue(a['preferred_80_complete_case_target_met'])
        for k in ('clinical_analysis_authorized','assay_technical_validation_verified','prognosis_validated'):
            self.assertFalse(a[k])

    def test_all_conflicting_duplicates_held(self):
        a=screen([v2(i) for i in range(60)]+[v2(0,{'irb_permits_reuse':'no'})])
        self.assertEqual(a['evaluable_observed_native_lung_fvc_n'],59)
        self.assertEqual(a['baseline_exclusion_reasons']['missing_or_duplicate_subject_key'],2)
        self.assertTrue(a['status'].startswith('HOLD'))

    def test_adult_plausibility(self):
        for value in ('-1','0','17','121','1000'):
            self.assertEqual(audit([row(1,{'age':value})])['evaluable_n'],0)

    def test_source_and_diagnosis_checks(self):
        for k in ('diagnosis_adjudicated','same_person_linkage_verified','preanalytics_reviewed','baseline_treatment_timing_verified'):
            self.assertEqual(screen([v2(1,{k:'unknown'})])['baseline_eligible_n'],0)

    def test_matrix_not_silently_pooled(self):
        for over in ({'specimen_matrix':'serum'},{'anticoagulant':'unknown'}):
            self.assertEqual(screen([v2(1,over)])['baseline_eligible_n'],0)

    def test_terminal_outcomes_not_stable(self):
        for event in ('DEATH_BEFORE_FVC','TRANSPLANT_BEFORE_FVC'):
            a=screen([v2(1,{'endpoint_status':event,'terminal_event_date':'2023-06-01',
                            'followup_fvc_date':'','followup_fvc_percent':''})])
            self.assertEqual(a['baseline_eligible_outcome_states'][event],1)
            self.assertEqual(a['evaluable_observed_native_lung_fvc_n'],0)

    def test_terminal_date_required(self):
        a=screen([v2(1,{'endpoint_status':'DEATH_BEFORE_FVC','followup_fvc_date':'','followup_fvc_percent':''})])
        self.assertEqual(a['baseline_eligible_outcome_states']['HOLD_TERMINAL_EVENT_DATE_OUTSIDE_REVIEW_WINDOW'],1)

    def test_postterminal_fvc_held(self):
        a=screen([v2(1,{'terminal_event_date':'2023-06-01'})])
        self.assertEqual(a['baseline_eligible_outcome_states']['HOLD_FVC_NOT_PRETERMINAL_FOLLOWUP'],1)

    def test_missing_not_zero_change(self):
        a=screen([v2(1,{'endpoint_status':'MISSING_FVC','followup_fvc_date':'','followup_fvc_percent':''})])
        self.assertEqual(a['baseline_eligible_outcome_states']['MISSING_FVC'],1)
        self.assertEqual(a['evaluable_observed_native_lung_fvc_n'],0)

    def test_contradictory_missing_held(self):
        a=screen([v2(1,{'endpoint_status':'MISSING_FVC'})])
        self.assertEqual(a['baseline_eligible_outcome_states']['HOLD_CONFLICTING_MISSING_AND_FVC_RECORD'],1)

    def test_fvc_quality_and_adjudication(self):
        for key in ('endpoint_adjudicated','native_lung_fvc_confirmed','fvc_quality_reviewed'):
            self.assertEqual(screen([v2(1,{key:'unknown'})])['evaluable_observed_native_lung_fvc_n'],0)

    def test_no_identifiers_or_clinical_values(self):
        a=screen([v2(987654321,{'baseline_fvc_percent':'78.123456'})])
        self.assertNotIn('SYNTH987654321',str(a))
        self.assertNotIn('78.123456',str(a))

    def test_invalid_csv_schema(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'x.csv'
            for s in ('subject_id,subject_id\nx,y\n','subject_id,other\nx\n','subject_id\nx,y\n'):
                path.write_text(s)
                with self.assertRaises(ValueError):load_manifest(path,['subject_id'])

    def test_sample_size_benchmark(self):
        self.assertEqual(sample_size_for_correlation(.30,tests=1),85)
        self.assertEqual(sample_size_for_correlation(.30),103)
        scenarios=planning()['scenarios']
        self.assertEqual(scenarios[1]['screen_n_at_80pct_baseline_and_80pct_followup'],161)
        self.assertGreater(scenarios[0]['complete_n'],scenarios[1]['complete_n'])
        self.assertLess(scenarios[2]['complete_n'],scenarios[1]['complete_n'])


if __name__=='__main__':unittest.main()
