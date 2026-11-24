#!/usr/bin/env python3
"""Audit publicly annotated GSE213761 clinical covariate imbalance, no individual rows published."""
import argparse
import itertools
import math
import re
from pathlib import Path
import pandas as pd
import numpy as np

FIELDS={
    'age':'age', 'bmi':'bmi', 'hurley stage':'hurley_stage',
    'abcess count':'abscess_count','abscess count':'abscess_count',
    'node count':'nodule_count','fistula count':'fistula_count',
    'smoking status':'smoking_status','gender':'gender','responder':'responder'}

def src_metadata_to_patient(df, matched_patient_ids):
    body=[]
    for _,row in df.iterrows():
        m=re.match(r'^(\d+)_Lesional_V1$',str(row['title']),flags=re.I)
        if not m:continue
        pat=m.group(1)
        if pat not in matched_patient_ids:continue
        record={'patient':pat}
        for col in df.columns:
            if 'characteristics' not in col:continue
            val=str(row[col]);parts=val.split(':',1)
            if len(parts)!=2:continue
            k=parts[0].strip().lower();value=parts[1].strip()
            if k in FIELDS:record[FIELDS[k]]=value
        body.append(record)
    out=pd.DataFrame(body)
    assert len(out)==14 and out.patient.nunique()==14
    assert out.responder.value_counts().to_dict()=={'YES':7,'NO':7}
    return out

def exact_p(vals,yes):
    n=len(vals);k=int(sum(yes));obs=float(np.mean(vals[yes])-np.mean(vals[~yes]));tail=0;tot=0
    for c in itertools.combinations(range(n),k):
        mask=np.zeros(n,dtype=bool);mask[list(c)]=True
        val=float(np.mean(vals[mask])-np.mean(vals[~mask]))
        if abs(val)>=abs(obs)-1e-10:tail+=1
        tot+=1
    return obs,tail/tot

def main():
    p=argparse.ArgumentParser();p.add_argument('--design',required=True);p.add_argument('--paired');p.add_argument('--out',required=True);a=p.parse_args()
    df=pd.read_csv(a.design,dtype=str)
    if a.paired:
        paired=pd.read_csv(a.paired,dtype={'patient':str})
        matched_ids=set(paired.patient)
    else:
        pre=set(df.title.str.extract(r'^(\d+)_Lesional_V1$',flags=re.I)[0].dropna())
        post=set(df.title.str.extract(r'^(\d+)_Lesional_V12$',flags=re.I)[0].dropna())
        matched_ids=pre&post
    clinical=src_metadata_to_patient(df,matched_ids)
    numbers=['age','bmi','abscess_count','nodule_count','fistula_count']
    out=[]
    yes=(clinical.responder=='YES').values
    for key in numbers:
        if key not in clinical.columns:continue
        ser=pd.to_numeric(clinical[key],errors='coerce')
        if ser.isnull().any():
            out.append({'variable':key,'status':'MISSING_NOT_TESTED','yes_n':int(ser[yes].notna().sum()),'no_n':int(ser[~yes].notna().sum())});continue
        vals=ser.to_numpy();delta,pval=exact_p(vals,yes)
        out.append({'variable':key,'status':'DESCRIPTIVE_ONLY','yes_n':7,'no_n':7,'yes_median':round(float(np.median(vals[yes])),3),
                    'no_median':round(float(np.median(vals[~yes])),3),'mean_delta_yes_minus_no':round(delta,3),
                    'unadjusted_exact_p':round(pval,5)})
    for key in ['gender','hurley_stage','smoking_status']:
        if key not in clinical.columns:continue
        for level in sorted(clinical[key].dropna().unique()):
            out.append({'variable':key,'level':level,'status':'CATEGORICAL_COUNTS_ONLY',
                        'yes_n':int(((clinical[key]==level)&yes).sum()),
                        'no_n':int(((clinical[key]==level)&(~yes)).sum())})
    target=Path(a.out);target.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(out).to_csv(target/'hs_response_public_clinical_balance.csv',index=False)
    md=['# GSE213761 baseline patient covariate audit','',
        '- All 14 people are represented by one baseline lesional biopsy.',
        '- 7 with GEO responder=YES, 7 with responder=NO.',
        '- Clinical variables transcribed from public GEO characteristics without exposing individual-level rows.',
        '- Permutation p values are unadjusted descriptive checks and do NOT establish randomized balance.',
        '- Small-sample baseline group differences may confound post-hoc response association.','',
        '|Variable/level|YES median (or n)|NO median (or n)|mean YES-NO / exact p|',
        '|---|---:|---:|---|']
    for a in out:
        var=a['variable']+(' : '+str(a['level']) if 'level' in a else '')
        y=a.get('yes_median',a.get('yes_n','–'));n=a.get('no_median',a.get('no_n','–'))
        v=f"{a.get('mean_delta_yes_minus_no','–')} / {a.get('unadjusted_exact_p','–')}" if 'unadjusted_exact_p' in a else 'descriptive'
        md.append(f'|{var}|{y}|{n}|{v}|')
    (target/'HS_CLINICAL_BALANCE_REPORT.md').write_text('\n'.join(md)+'\n',encoding='utf8')
    print(pd.DataFrame(out).to_string(index=False))

if __name__=='__main__':main()