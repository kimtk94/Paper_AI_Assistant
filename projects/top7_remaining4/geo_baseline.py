#!/usr/bin/env python3
"""Reproducible *public-data* GEO feasibility/QC for the remaining four Novogene prospects.
No service outreach, no patients' private data, no inferential claims without donor mapping.
"""
import argparse, collections, csv, datetime, gzip, io, json, math, re, time, urllib.request, urllib.error, hashlib
from pathlib import Path

DATASETS = {
 'hs_therapy': {'accessions':['GSE213761'],'caveat':'Adalimumab/anti-TNF longitudinal skin RNA, not secukinumab/anti-IL17. Tissue and patient timepoint pairing must be checked.'},
 'musk_mg': {'accessions':['GSE227835','GSE297547'],'caveat':'Neither dataset can be called a MuSK+ treatment-response cohort without explicit MuSK status; GSE227835 includes AChR+ and seronegative MG.'},
 'nmosd': {'accessions':['GSE242232'],'caveat':'Stimulated antigen-enriched CD154+ memory T cells, not unsorted PBMC/plasma and not relapse prediction.'},
 'hs_spatial': {'accessions':['GSE294009','GSE249729'],'caveat':'Spatial technologies and sample units differ; 156 ROIs are not 156 patients; GSE249729 is only 2 HS lesions and 1 control.'},
}
SUPPS = {'GSE213761':['GSE213761_RawCountsGEO.csv.gz'], 'GSE242232':['GSE242232_RawExpressionCountTable.csv.gz']}
GENES={'hs_therapy':['IL17A','IL17F','IL1B','TNF','CXCL8','CXCL13','S100A8','S100A9','KRT17','CCL20'],
       'nmosd':['FOXP3','PDCD1','CTLA4','LAG3','TIGIT','TOX','CD40LG','IL21'],
       'musk_mg':['CD22','MS4A1','CD79A','CD74','CFD','IL21','CXCR5'],
       'hs_spatial':['CXCL13','CD79A','MS4A1','IL17A','IL17F','CCL20','COL1A1','KRT17']}
UA='PublicResearchDataAudit/1.0 (+https://www.ncbi.nlm.nih.gov/geo/)'

def urlopen(url, maximum=32_000_000, timeout=50, retries=2):
    error=None
    for i in range(retries+1):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':UA, 'Accept':'*/*'})
            with urllib.request.urlopen(req,timeout=timeout) as fh:
                data=fh.read(maximum+1)
                if len(data)>maximum: return None, {'status':'TOO_LARGE','gt_bytes':maximum,'url':url}
                return data, {'status':'OK','bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'url':url}
        except (urllib.error.HTTPError,urllib.error.URLError,OSError,TimeoutError) as e:
            error=str(e)
            time.sleep(0.6*(i+1))
    return None,{'status':'ERROR','error':error[:300],'url':url}

def series_path(acc):
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'

def read_series(acc):
    base=series_path(acc)
    url=f'{base}/matrix/{acc}_series_matrix.txt.gz'
    blob,q=urlopen(url,maximum=32_000_000)
    rows=[]
    if blob:
        try:
            txt=gzip.decompress(blob).decode('utf8','replace')
            for line in txt.splitlines():
                if not line.startswith('!Sample_'):continue
                row=next(csv.reader([line],delimiter='\t'))
                if len(row)>1:rows.append((row[0],row[1:]))
        except Exception as e:q.update({'parse_error':str(e)[:250]})
    values=collections.defaultdict(list)
    for key,items in rows:values[key].append(items)
    accessions=values.get('!Sample_geo_accession',[[]])[0]
    titles=values.get('!Sample_title',[[]])[0]
    n=max(len(accessions),len(titles))
    cols={'accession':accessions,'title':titles}
    for k,items in values.items():
        if k not in ('!Sample_geo_accession','!Sample_title') and any(w in k.lower() for w in ('characteristics','source_name','description','relation')):
            for i,item in enumerate(items):cols[f'{k}_{i+1}']=item
    samples=[]
    for idx in range(n):
        samples.append({k:(v[idx] if len(v)>idx else '') for k,v in cols.items()})
    return samples,{'fetch':q,'n_samples':len(samples),'metadata_keys':list(cols)[:100]}

def get_series_page(acc):
    url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}'
    data,q=urlopen(url,maximum=3_000_000,timeout=35)
    if not data:return {},q
    text=data.decode('utf8','replace')
    patterns=[r'Samples\s*\((\d+)\)',r'Overall design\s*([^<]{0,200})',r'Supplementary file\s*Size']
    return {'gsm_in_html':len(set(re.findall(r'GSM\d{5,}',text))),
            'sample_count_text':re.findall(patterns[0],text)[:1],
            'small_files':sorted(set(re.findall(r'\b'+acc+r'[^<"\s]{0,90}\.(?:csv\.gz|txt\.gz|tar|zip|xlsx|h5ad|rds)',text,re.I)))[:50],
            'mentions_MuSK':bool(re.search('MuSK',text,re.I))},q

def group_label(sample):
    t=(sample.get('title','')+' '+' '.join(str(v) for k,v in sample.items() if 'characteristics' in k)).lower()
    if 'seronegative' in t:return 'seronegative'
    if 'achr-positive' in t or 'achr positive' in t:return 'achr_positive'
    if 'healthy control' in t or 'normal' in t:return 'control'
    if 'nmosd' in t or 'aquaporin' in t:return 'nmosd'
    if 'musk' in t:return 'musk'
    if 'hs' in t or 'lesional' in t:return 'hs_or_lesional'
    return 'other_or_unspecified'

def analyze_small_matrix(acc,role,out):
    rows={}
    for name in SUPPS.get(acc,[]):
        url=f'{series_path(acc)}/suppl/{name}'
        blob,status=urlopen(url,maximum=17_000_000,timeout=60)
        item={'download':status,'filename':name}
        if not blob:
            rows[name]=item;continue
        try:
            txt=gzip.decompress(blob).decode('utf-8-sig','replace')
            header=next(csv.reader(io.StringIO(txt[:100_000]),delimiter=','))
            item['header_preview']=header[:30]
            item['header_count']=len(header)
            item['line_count']=txt.count('\n')
            item['uncompressed_bytes']=len(txt)
            if len(header)>1:
                vals=header[1:]
                item['gsm_col_matches']=sum(bool(re.fullmatch(r'GSM\d+',v.strip().strip('"'))) for v in vals)
            # Controlled subset of candidate genes, ONLY descriptive when sample or group mapping is unproven.
            cand=set(GENES[role])
            matches=[]
            r=csv.reader(io.StringIO(txt),delimiter=',')
            h=next(r,None)
            for row in r:
                if row and any(row[0].strip('"').upper()==g.upper() or row[0].strip('"').upper().endswith('|'+g.upper()) for g in cand):
                    matches.append({'gene':row[0][:100],'values_first_8':row[1:9]})
                    if len(matches)>=30:break
            item['target_gene_preview']=matches
            # Small numeric matrices retained in CI memory; do not upload the raw matrix artifact.
            if acc=='GSE213761':
                from io import StringIO
                try:
                    import pandas as pd
                    head=pd.read_csv(StringIO(txt),nrows=4)
                    item['column_names_first']=list(head.columns[:30])
                    item['numeric_columns_first']=list(head.select_dtypes(include='number').columns[:12])
                    item['row_shape']=list(head.shape)
                except Exception as e:item['small_matrix_parse_warning']=str(e)[:200]
        except Exception as e:item['parse_error']=str(e)[:350]
        rows[name]=item
    return rows

def write_csv(path, rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def run(role,out):
    conf=DATASETS[role]
    results={'role':role,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'scientific_caveat':conf['caveat'],'accessions':{},'outreach_authorized':False}
    for acc in conf['accessions']:
        samples,q=read_series(acc)
        page,pq=get_series_page(acc)
        group_counts=dict(collections.Counter(group_label(s) for s in samples))
        title_counts=collections.Counter(s.get('title','').split('_')[0] for s in samples)
        info={'series_matrix':q,'html_page':pq,'html_inventory':page,
              'group_titles':group_counts,'sample_title_head':[x.get('title','') for x in samples[:12]],
              'distinct_title_prefixes':len(title_counts),
              'small_public_expression':analyze_small_matrix(acc,role,out)}
        # HS patient/time design observed from titles only, no causal or efficacy claim.
        if acc=='GSE213761':
            parsed=[]
            for s in samples:
                m=re.match(r'^(\d+)_(Lesional|Nonlesional)_(V1|V12)$',s.get('title',''),re.I)
                if m:parsed.append((m.group(1),m.group(2).lower(),m.group(3).upper()))
            pairindex=collections.defaultdict(set)
            for patient,tissue,timepoint in parsed:pairindex[(patient,tissue)].add(timepoint)
            info['design_qc']={'valid_titles':len(parsed),'unique_patient_ids':len({x[0] for x in parsed}),
                'tissue_prepost_pairs':sum(v=={'V1','V12'} for v in pairindex.values()),
                'lesional_prepost_pairs':sum(t=='lesional' and v=={'V1','V12'} for (_,t),v in pairindex.items()),
                'nonlesional_prepost_pairs':sum(t=='nonlesional' and v=={'V1','V12'} for (_,t),v in pairindex.items())}
        if samples:write_csv(out/f'{acc}_sample_design.csv',samples)
        results['accessions'][acc]=info
    results['decision']='READY_FOR_FOCUSED_COMPARISON' if role=='hs_therapy' and results['accessions']['GSE213761'].get('design_qc',{}).get('lesional_prepost_pairs',0)>=5 else 'HOLD_FOR_PROPER_MATCHING_OR_PROCESSED_MATRIX'
    if role=='musk_mg':results['decision']='HOLD_NO_CONFIRMED_MUSK_SPECIFIC_COHORT'
    if role=='nmosd':results['decision']='MECHANISM_CONTEXT_ONLY_NOT_RELAPSE_PREDICTION'
    if role=='hs_spatial':results['decision']='HOLD_DIFFERENT_PLATFORM_PATIENT_ROI_MAPPING'
    (out/'audit.json').write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding='utf8')
    report=[f'# Remaining Top-7 public-data baseline: {role}', '',f"Decision: **{results['decision']}**",'',f"Design caution: {conf['caveat']}", '']
    for acc,v in results['accessions'].items():
        report.extend([f'## {acc}',f"- sample metadata: {v['series_matrix']['n_samples']} rows; {v['series_matrix']['fetch']['status']}",f"- group/title buckets: {v['group_titles']}",f"- direct processed candidates: {[(k,z['download']['status'],z.get('header_count'),z.get('line_count')) for k,z in v['small_public_expression'].items()]}",f"- study: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}",''])
        if 'design_qc' in v:report.append(f"- matched patient/tissue pairs: {v['design_qc']}")
    report.extend(['## Interpretation constraints','- Public cohorts are different from the named Korean clinician cohort.','- Cells and ROI are not independent patients; distinguish biological donors and assay replicates.','- Different drugs, mechanisms, tissue, and treatment timepoints require their own analysis design.','- No claim about predictive model, causal treatment effect, or clinical biomarker without matched outcomes and independent validation.','- Territory and SFDC ownership must be checked before outreach.',''])
    (out/'REPORT.md').write_text('\n'.join(report),encoding='utf8')
    print(json.dumps({'role':role,'decision':results['decision'],'counts':{k:v['series_matrix']['n_samples'] for k,v in results['accessions'].items()}},ensure_ascii=False))
    return results

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--role',choices=list(DATASETS),required=True);ap.add_argument('--out',required=True)
    args=ap.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True);run(args.role,out)