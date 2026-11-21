#!/usr/bin/env python3
"""Strict, aggregate-only, optional clinical specimen feasibility audit.

DO NOT COMMIT PATIENT DATA. Use only secure local inputs; CI runs synthetic tests.
No external API connections, no subject-level output. This is not a diagnosis/prognosis model.
"""
import argparse
import csv
import datetime as dt
import json
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

def yes(x):
    return str(x or '').strip().lower() in YES

def parse_date(x):
    try: return dt.date.fromisoformat(str(x).strip())
    except (ValueError,TypeError): return None

def as_float(x):
    try: return float(x)
    except (ValueError,TypeError): return None

def audit(rows, min_pilot=60):
    rows=list(rows); reasons=Counter(); valid=[]; seen=set()
    for r in rows:
        sid=(r.get('subject_id') or '').strip()
        if not sid or sid in seen:
            reasons['missing_or_duplicate_subject_key'] += 1
            continue
        seen.add(sid)
        if not yes(r.get('irb_permits_reuse')):
            reasons['IRB_or_consent_not_confirmed'] += 1; continue
        if (r.get('diagnosis') or '').strip().upper() != 'IPF':
            reasons['diagnosis_not_IPF']+=1; continue
        if not yes(r.get('baseline_plasma_available')):
            reasons['baseline_plasma_unavailable']+=1; continue
        p, b, f = [parse_date(r.get(n)) for n in ['plasma_date','baseline_fvc_date','followup_fvc_date']]
        if not p or not b or not f:
            reasons['missing_date']+=1; continue
        if abs((p-b).days)>30:
            reasons['plasma_FVC_timing_over_30d']+=1; continue
        follow=(f-b).days
        if not 270 <= follow <= 450:
            reasons['followup_not_9_to_15_months']+=1; continue
        base, final=as_float(r.get('baseline_fvc_percent')), as_float(r.get('followup_fvc_percent'))
        if base is None or final is None or not(0 < base <= 150 and 0 < final <= 150):
            reasons['FVC_missing_or_implausible']+=1; continue
        missing=[x for x in ['age','sex','smoking','antifibrotic','batch'] if not str(r.get(x) or '').strip()]
        if missing:
            reasons['missing_critical_covariate']+=1;continue
        if yes(r.get('hemolysis_flag')):
            reasons['hemolysis_review_required']+=1;continue
        valid.append((final-base)*365.25/follow)
    n=len(valid); total=len(rows)
    if total==0: status='HOLD_NO_DATA'
    elif n < min_pilot: status='HOLD_FEWER_THAN_60_EVALUABLE'
    elif n < 80: status='GO_SMALL_EXPLORATORY_ONLY'
    else: status='GO_EXPLORATORY_FEASIBILITY_ONLY'
    return {'status':status,'screened_n':total,'evaluable_n':n,'evaluable_percent':round(100*n/total,1) if total else 0,
            'target_exploratory_min_n':min_pilot,'preferred_exploratory_n':80,
            'exclusion_reasons':dict(sorted(reasons.items())),
            'scientific_limitations':['No patient-level data exported','Not a power calculation',
            'Candidate protein detectability unknown','No prognostic model validation',
            'Sample screening is not proof that plasma and FVC data may lawfully be reused'],
            'required_next_checks':['assay pilot FCN3/GPX3','source/covariate and treatment balance',
            'IRB and data handling','formal power/event-rate review','external cohort availability']}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',required=True,help='Local, de-identified CSV only; never upload to GitHub')
    p.add_argument('--out',default='audit_aggregate.json')
    args=p.parse_args()
    with open(args.manifest,encoding='utf-8-sig',newline='') as f:
        r=csv.DictReader(f)
        missing=sorted(set(COLUMNS)-set(r.fieldnames or []))
        if missing: raise ValueError('Missing required columns: '+','.join(missing))
        result=audit(list(r))
    Path(args.out).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','screened_n','evaluable_n','evaluable_percent')},ensure_ascii=False))

if __name__=='__main__':main()