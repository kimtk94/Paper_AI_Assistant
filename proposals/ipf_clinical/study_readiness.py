#!/usr/bin/env python3
"""Local-only aggregate clinical readiness and illustrative planning calculations.

A manifest is never an approval, and no software result permits recruitment,
sample transfer, clinical testing or prognostic claims. No measurements/IDs emitted.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter
from pathlib import Path
from statistics import NormalDist
from feasibility import (COLUMNS,baseline_reason,followup_reason,load_manifest,
                         normalize_token,parse_date,subject_key,yes)

V2_COLUMNS=COLUMNS+[
    'diagnosis_adjudicated','same_person_linkage_verified','specimen_matrix',
    'anticoagulant','preanalytics_reviewed','baseline_treatment_timing_verified',
    'endpoint_status','endpoint_adjudicated','native_lung_fvc_confirmed',
    'fvc_quality_reviewed','terminal_event_date',
]
OUTCOME_STATES={'OBSERVED_NATIVE_LUNG_FVC','DEATH_BEFORE_FVC','TRANSPLANT_BEFORE_FVC',
                'MISSING_FVC','UNKNOWN'}
TERMINAL={'DEATH_BEFORE_FVC','TRANSPLANT_BEFORE_FVC'}


def sample_size_for_correlation(r, alpha=0.05, power=0.80, tests=2):
    """Approximate Fisher-z two-sided correlation benchmark, Bonferroni family alpha.

    Not multivariable regression power; r is assumed, not estimated from disease
    case-control log fold-change or the lung RNA cohorts.
    """
    if not 0<abs(r)<1 or not 0<alpha<1 or not 0.5<power<1 or not isinstance(tests,int) or tests<1:
        raise ValueError('invalid planning assumptions')
    z=NormalDist().inv_cdf(1-alpha/(2*tests))+NormalDist().inv_cdf(power)
    return math.ceil(3+(z/math.atanh(abs(r)))**2)


def planning():
    rows=[]
    for r in (0.25,0.30,0.40):
        n=sample_size_for_correlation(r)
        rows.append({'assumed_abs_correlation':r,'complete_n':n,
                     'screen_n_at_80pct_baseline_and_80pct_followup':math.ceil(n/(0.8*0.8))})
    return {'scope':'ILLUSTRATIVE_UNADJUSTED_CORRELATION_NOT_FORMAL_POWER',
            'family_alpha':0.05,'two_sided_per_test_alpha':0.025,'power':0.80,
            'primary_markers':['FCN3','GPX3'],'scenarios':rows,
            'single_marker_abs_r_0_3_n':sample_size_for_correlation(0.30,tests=1),
            'n_screened_for_80_complete_at_64pct_coverage':125,
            'expected_complete_from_107_if_64pct_coverage':107*0.8*0.8,
            'limitations':['Independent complete cases assumed; actual attrition may be informative.',
                           'Adjusted partial effect, residual variance, covariate degrees of freedom, batch and death/transplant require separate planning.',
                           'No power estimate from cross-sectional abundance effects.']}


def baseline_v2_reason(row):
    reason=baseline_reason(row)
    if reason:return reason
    for field in ('diagnosis_adjudicated','same_person_linkage_verified',
                  'preanalytics_reviewed','baseline_treatment_timing_verified'):
        if not yes(row.get(field)):return field+'_not_confirmed'
    # Proposed EDTA-plasma primary cohort; alternative matrices require a revised
    # method and comparability review, not silent pooling with EDTA plasma.
    if normalize_token(row.get('specimen_matrix'))!='plasma' or normalize_token(row.get('anticoagulant'))!='edta':
        return 'specimen_matrix_or_anticoagulant_requires_review'
    return None


def outcome_category(row):
    state=str(row.get('endpoint_status') or '').strip().upper()
    if state not in OUTCOME_STATES:
        return 'HOLD_UNKNOWN_ENDPOINT_STATE'
    if state=='UNKNOWN':
        return 'HOLD_UNKNOWN_ENDPOINT_STATE'
    if not yes(row.get('endpoint_adjudicated')):
        return 'HOLD_ENDPOINT_NOT_ADJUDICATED'
    b=parse_date(row.get('baseline_fvc_date'));p=parse_date(row.get('plasma_date'))
    event_raw=str(row.get('terminal_event_date') or '').strip()
    event=parse_date(event_raw) if event_raw else None
    if event_raw and event is None:
        return 'HOLD_INVALID_TERMINAL_EVENT_DATE'
    if event and event<=max(b,p):
        return 'HOLD_TERMINAL_EVENT_BEFORE_OR_AT_BASELINE'
    if state in TERMINAL:
        if not event or not 0<(event-b).days<=450:
            return 'HOLD_TERMINAL_EVENT_DATE_OUTSIDE_REVIEW_WINDOW'
        f=parse_date(row.get('followup_fvc_date'))
        # Even if an earlier FVC exists, do not let contradictory qualifying
        # native-lung followup/event records be accepted without adjudication.
        if f and 270<=(f-b).days<=450 and not followup_reason(row):
            return 'HOLD_CONFLICTING_TERMINAL_AND_FVC_RECORD'
        return state
    if state=='MISSING_FVC':
        if event and 0<(event-b).days<=450:
            return 'HOLD_TERMINAL_EVENT_REQUIRES_ADJUDICATION'
        if not followup_reason(row):
            return 'HOLD_CONFLICTING_MISSING_AND_FVC_RECORD'
        return 'MISSING_FVC'
    reason=followup_reason(row)
    if reason:
        return 'HOLD_'+reason.upper()
    f=parse_date(row.get('followup_fvc_date'))
    if f<=max(b,p) or (event and event<=f):
        return 'HOLD_FVC_NOT_PRETERMINAL_FOLLOWUP'
    if not yes(row.get('native_lung_fvc_confirmed')) or not yes(row.get('fvc_quality_reviewed')):
        return 'HOLD_FVC_QUALITY_OR_NATIVE_LUNG_NOT_CONFIRMED'
    return 'OBSERVED_NATIVE_LUNG_FVC'


def screen(rows):
    rows=list(rows);keys=Counter(subject_key(r) for r in rows if subject_key(r))
    excluded=Counter();outcomes=Counter()
    for row in rows:
        sid=subject_key(row)
        if not sid or keys[sid]!=1:
            excluded['missing_or_duplicate_subject_key']+=1;continue
        reason=baseline_v2_reason(row)
        if reason:
            excluded[reason]+=1;continue
        outcomes[outcome_category(row)]+=1
    baseline_n=sum(outcomes.values())
    n=outcomes['OBSERVED_NATIVE_LUNG_FVC']
    return {'status':('ELIGIBLE_FOR_EXPLORATORY_DESIGN_REVIEW_ONLY' if n>=60 else 'HOLD_FEWER_THAN_60_OBSERVED_NATIVE_LUNG_FVC'),
            'screened_rows':len(rows),'unique_nonblank_subject_keys':len(keys),
            'baseline_eligible_n':baseline_n,'evaluable_observed_native_lung_fvc_n':n,
            'preferred_80_complete_case_target_met':n>=80,
            'baseline_exclusion_reasons':dict(sorted(excluded.items())),
            'baseline_eligible_outcome_states':dict(sorted(outcomes.items())),
            'cohort_counts_are_manifest_reported_not_independently_verified':True,
            'clinical_analysis_authorized':False,'assay_technical_validation_verified':False,
            'prognosis_validated':False,
            'next_gate':'Independent data-access/IRB, specimen SOP, assay qualification, outcome/attrition adjudication, statistical and external validation review.',
            'limitations':['Only observed native-lung FVC records enter the preliminary complete-case count.',
                           'Death/transplant and missing/unknown outcomes are not stable FVC, zero change or verified non-events.',
                           'Aggregate counts are still sensitive research information; keep output in an authorized local area.']}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',help='Optional authorized local v2 CSV; never commit or upload')
    p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    result={'planning':planning(),'clinical_screen':screen(load_manifest(args.manifest,V2_COLUMNS)) if args.manifest else
            {'status':'HOLD_NO_CLINICAL_MANIFEST','clinical_analysis_authorized':False,'prognosis_validated':False}}
    # Write only to the explicit path chosen by the authorized operator.
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['clinical_screen']['status'],'clinical_analysis_authorized':False}))


if __name__=='__main__':main()
