#!/usr/bin/env python3
"""Audit GSE135893 donor/source/disease confounding before interpreting candidate gene shifts."""
from __future__ import annotations
import argparse,csv,json,datetime as dt
from collections import Counter,defaultdict
from pathlib import Path
from gse135893_metadata import fetch,read_gz_table

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='output/ipf_gse135893_source');o=Path(p.parse_args().out);o.mkdir(parents=True,exist_ok=True)
    raw,meta=fetch('GSE135893_IPF_metadata.csv.gz')
    hdr,rows=read_gz_table(raw,',',has_header=True); col={n:i for i,n in enumerate(hdr)}
    required=['Diagnosis','Sample_Name','Sample_Source','orig.ident']
    if any(x not in col for x in required):raise ValueError(f'Missing {set(required)-set(col)}')
    donors=defaultdict(lambda:{'cells':0,'libraries':set()})
    for row in rows:
        k=(row[col['Sample_Name']],row[col['Diagnosis']],row[col['Sample_Source']])
        donors[k]['cells']+=1
        donors[k]['libraries'].add(row[col['orig.ident']])
    by_donor=defaultdict(set)
    for donor,diag,source in donors:by_donor[donor].add((diag,source))
    collision={k:sorted(v) for k,v in by_donor.items() if len(v)>1}
    if collision:raise ValueError('Donor diagnosis/source collision '+str(collision)[:500])
    summarized=[]; source_count=Counter(); cell_count=Counter()
    for (donor,diag,src),v in sorted(donors.items()):
        summarized.append({'donor':donor,'diagnosis':diag,'source':src,'n_cells':v['cells'],'n_libraries':len(v['libraries'])})
        source_count[(diag,src)]+=1;cell_count[(diag,src)]+=v['cells']
    with (o/'donor_source_mapping.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(summarized[0]));wr.writeheader();wr.writerows(summarized)
    status='GO_SOURCE_STRATIFIABLE' if any(source_count[('IPF',s)]>=3 and source_count[('Control',s)]>=3 for _,s in source_count) else 'HOLD_SOURCE_STRATIFICATION'
    result={'accession':'GSE135893','source_sha256':meta['sha256'],'checked_at_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
        'diagnosis_source_donors':[{'diagnosis':d,'source':s,'n_donors':n,'n_cells':cell_count[(d,s)]} for (d,s),n in sorted(source_count.items())],
        'donor_total':len(summarized),'collision_n':len(collision),'status':status,
        'note':'This source audit is not proof that technical batch effects are fully corrected; no age/sex in public annotation columns.'}
    (o/'source_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    report=['# GSE135893 collection source confounding audit', '',f"- Status: **{status}**",f"- Unique donors: {len(summarized)}",'','| Diagnosis | Source | Donors | Cells |','|---|---|---:|---:|']
    for (diag,src),v in sorted(source_count.items()): report.append(f'| {diag} | {src} | {v} | {cell_count[(diag,src)]} |')
    report+=['','## Interpretation','- Prioritized IPF vs Control comparison should be stratified by Sample_Source if >=3 donors each within shared source.','- Unmeasured center / biological differences can remain even with source stratification.']
    (o/'SOURCE_REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__': main()