#!/usr/bin/env python3
"""Strict, aggregate-only, optional clinical specimen feasibility audit.

DO NOT COMMIT PATIENT DATA. Use only secure local inputs; CI runs synthetic tests.
No external API connections, no subject-level output. This is not a diagnosis/prognosis model.
"""
import argparse
import csv
import datetime as dt
import json
import math
from collections import Counter
from pathlib import Path

COLUMNS = [
    'subject_id', 'irb_permits_reuse', 'diagnosis', 'baseline_plasma_available',
    'plasma_date', 'baseline_fvc_date', 'baseline_fvc_percent',
    'followup_fvc_date', 'followup_fvc_percent', 'age', 'sex',
    'smoking', 'antifibrotic', 'baseline_dlco_percent',
    'batch', 'hemolysis_flag',
]
YES = {'yes', 'y', 'true', '1'}
NO = {'no', 'n', 'false', '0'}
MISSING_TOKENS = {
    '', 'na', 'n/a', 'n.a.', 'nan', 'null', 'none', 'unknown', 'unk',
    'not known', 'not available', 'not recorded', 'not assessed',
    'not collected', 'missing', 'pending', 'tbd', '-', '--', '?',
}

def normalize_token(x):
    return '' if x is None else ' '.join(str(x).strip().lower().split())

def yes(x):
    return normalize_token(x) in YES

def missing_covariate(field, value):
    token = normalize_token(value)
    # An explicit absence of smoking/treatment is a recorded covariate, not NA.
    if field in {'smoking', 'antifibrotic'} and token == 'none':
        return False
    return token in MISSING_TOKENS

def parse_date(x):
    try: return dt.date.fromisoformat(str(x).strip())
    except (ValueError,TypeError): return None

def as_float(x):
    try: return float(x)
    except (ValueError,TypeError): return None

def subject_key(row):
    return (row.get('subject_id') or '').strip()


def baseline_reason(r):
    """Preliminary baseline screen; flags are not independent IRB or assay proof."""
    if not yes(r.get('irb_permits_reuse')):
        return 'IRB_or_consent_not_confirmed'
    if (r.get('diagnosis') or '').strip().upper() != 'IPF':
        return 'diagnosis_not_IPF'
    if not yes(r.get('baseline_plasma_available')):
        return 'baseline_plasma_unavailable'
    p, b = [parse_date(r.get(n)) for n in ('plasma_date','baseline_fvc_date')]
    if not p or not b:
        return 'missing_date'
    if abs((p-b).days)>30:
        return 'plasma_FVC_timing_over_30d'
    base=as_float(r.get('baseline_fvc_percent'))
    if base is None or not 0 < base <= 150:
        return 'FVC_missing_or_implausible'
    if any(missing_covariate(x,r.get(x)) for x in ('age','sex','smoking','antifibrotic','batch')):
        return 'missing_critical_covariate'
    age=as_float(r.get('age'))
    if age is None or not math.isfinite(age):
        return 'age_not_finite_numeric'
    # Proposed adult cohort / plausibility guard, not a clinical diagnostic rule.
    if not 18 <= age <= 120:
        return 'age_outside_proposed_adult_range'
    hemolysis=normalize_token(r.get('hemolysis_flag'))
    if hemolysis in YES:
        return 'hemolysis_review_required'
    if hemolysis not in NO:
        return 'hemolysis_not_confirmed'
    return None


def followup_reason(r):
    b,f=[parse_date(r.get(n)) for n in ('baseline_fvc_date','followup_fvc_date')]
    if not b or not f:
        return 'missing_date'
    if not 270 <= (f-b).days <= 450:
        return 'followup_not_9_to_15_months'
    final=as_float(r.get('followup_fvc_percent'))
    if final is None or not 0 < final <= 150:
        return 'FVC_missing_or_implausible'
    return None


def audit(rows, min_pilot=60):
    rows=list(rows)
    if not isinstance(min_pilot,int) or not 1 <= min_pilot <= 80:
        raise ValueError('invalid exploratory minimum')
    reasons=Counter()
    keys=Counter(subject_key(r) for r in rows if subject_key(r))
    n=0
    for r in rows:
        sid=subject_key(r)
        if not sid or keys[sid]!=1:
            # A conflicting duplicate invalidates ALL rows for that key until
            # independently resolved; never prefer the first clinical record.
            reasons['missing_or_duplicate_subject_key']+=1
            continue
        reason=baseline_reason(r) or followup_reason(r)
        if reason:
            reasons[reason]+=1
            continue
        n+=1
    total=len(rows)
    if total==0: status='HOLD_NO_DATA'
    elif n < min_pilot: status='HOLD_FEWER_THAN_60_EVALUABLE' if min_pilot==60 else 'HOLD_BELOW_CONFIGURED_MINIMUM'
    elif n < 80: status='GO_SMALL_EXPLORATORY_ONLY'
    else: status='GO_EXPLORATORY_FEASIBILITY_ONLY'
    return {'status':status,'screened_n':total,'unique_nonblank_subject_keys':len(keys),
            'duplicate_subject_rows_excluded':sum(v for v in keys.values() if v>1),
            'evaluable_n':n,'evaluable_percent':round(100*n/total,1) if total else 0,
            'target_exploratory_min_n':min_pilot,'preferred_exploratory_n':80,
            'exclusion_reasons':dict(sorted(reasons.items())),
            'screen_scope':'PRELIMINARY_MANIFEST_ONLY_USE_V2_FOR_OUTCOME_REVIEW',
            'clinical_analysis_authorized':False,'assay_technical_validation_verified':False,
            'scientific_limitations':['No patient-level data exported','Not a power calculation',
            'Candidate protein detectability unknown','No prognostic model validation',
            'Sample screening is not proof that plasma and FVC data may lawfully be reused'],
            'required_next_checks':['assay pilot FCN3/GPX3','source/covariate and treatment balance',
            'IRB and data handling','formal power/event-rate review','external cohort availability']}


def load_manifest(path, required_columns):
    with open(path,encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f,strict=True)
        header=r.fieldnames or []
        if len(header)!=len(set(header)) or any(not x.strip() for x in header):
            raise ValueError('Duplicate or blank CSV headers')
        if not set(required_columns).issubset(header):
            raise ValueError('Missing required clinical manifest columns')
        rows=list(r)
        if any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError('Ragged clinical manifest rows')
        return rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',required=True,help='Local, de-identified CSV only; never upload to GitHub')
    p.add_argument('--out',default='audit_aggregate.json')
    args=p.parse_args()
    result=audit(load_manifest(args.manifest,COLUMNS))
    Path(args.out).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','screened_n','evaluable_n','evaluable_percent')},ensure_ascii=False))

if __name__=='__main__':main()
