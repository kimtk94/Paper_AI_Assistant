#!/usr/bin/env python3
"""GSE294009 ROI/donor and file access audit. Metadata only; no ROI pseudo-replication."""
import csv, gzip, html, io, json, os, re, collections, urllib.request, urllib.error, time, argparse, hashlib
from pathlib import Path

ACC='GSE294009'
BASE=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{ACC[:-3]}nnn/{ACC}'
UA='HS-GeoMx-Scientific-Audit/2026; public GEO metadata'

def get(url,limit=35_000_000, timeout=75):
    error=None
    for attempt in range(2):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':UA,'Accept':'*/*'})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body=resp.read(limit+1)
                if len(body)>limit:return None,{'status':'LARGE','lower_bound_bytes':limit,'url':url,'http':resp.status}
                return body,{'status':'OK','bytes':len(body),'sha256':hashlib.sha256(body).hexdigest(),'url':url,'http':resp.status}
        except Exception as e:
            error=f'{type(e).__name__}: {e}'[:230];time.sleep(attempt*.8)
    return None,{'status':'ERROR','detail':error,'url':url}

def head(url):
    try:
        req=urllib.request.Request(url,method='HEAD',headers={'User-Agent':UA})
        with urllib.request.urlopen(req, timeout=30) as resp:
            return {'status':'OK','size':resp.headers.get('Content-Length'),'type':resp.headers.get('Content-Type'),'http':resp.status}
    except Exception as e: return {'status':'ERROR','detail':str(e)[:150]}

def parse_series_matrix(b):
    rows={}
    for line in gzip.decompress(b).decode('utf8','replace').splitlines():
        if not line.startswith('!Sample_'):continue
        arr=next(csv.reader([line],delimiter='\t'))
        if len(arr)>1:rows.setdefault(arr[0],[]).append(arr[1:])
    ids=rows.get('!Sample_geo_accession',[[]])[0]
    results=[{'gsm':v} for v in ids]
    for k,vv in rows.items():
        if k=='!Sample_geo_accession':continue
        if not any(z in k for z in ('title','characteristics','source_name','description','supplementary_file','relation','data_row_count')):continue
        for j,vals in enumerate(vv):
            key=k.replace('!Sample_','').replace('_ch1','')+('_'+str(j) if len(vv)>1 else '')
            for i,item in enumerate(vals[:len(results)]):results[i][key]=item
    return results

def parse_soft(b):
    samples={};gsm=None
    for line in gzip.decompress(b).decode('utf8','replace').splitlines():
        if line.startswith('^SAMPLE = '):
            gsm=line.split('=',1)[1].strip();samples[gsm]={}
        elif gsm and line.startswith('!Sample_'):
            key,value=line.split(' = ',1) if ' = ' in line else (line,'')
            if any(t in key for t in ('supplementary_file','relation','characteristics_ch1','description','title','source_name','data_row_count','table_begin','table_end')):
                samples[gsm].setdefault(key[8:],[]).append(value[:800])
    return samples

def write_csv(path, data):
    if not data:return
    keys=list(dict.fromkeys(k for row in data for k in row))
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(data)

def main():
    args=argparse.ArgumentParser();args.add_argument('--out',default='output/hs_spatial_gse294009');a=args.parse_args()
    out=Path(a.out);out.mkdir(exist_ok=True,parents=True)
    audit={'accession':ACC,'base':BASE,'results':{}}
    matrix,sm=get(BASE+'/matrix/'+ACC+'_series_matrix.txt.gz')
    audit['results']['series_matrix']=sm
    if matrix:
        rows=parse_series_matrix(matrix)
        write_csv(out/'roi_metadata.csv',rows)
        audit['results']['n_GSM']=len(rows)
        audit['results']['keys']=list(rows[0]) if rows else []
        audit['results']['titles_first']= [v.get('title','') for v in rows[:12]]
        attrs=collections.defaultdict(collections.Counter)
        for r in rows:
            for k,v in r.items():
                if k.startswith('characteristics'):attrs[k][v]+=1
        audit['results']['characteristics_counts']={k:{'nunique':len(v),'top':v.most_common(16)} for k,v in attrs.items()}
    else:rows=[]
    soft,ss=get(BASE+'/soft/'+ACC+'_family.soft.gz',limit=40_000_000)
    audit['results']['soft_download']=ss
    if soft:
        sf=parse_soft(soft)
        audit['results']['soft_n_samples']=len(sf)
        file_refs=[];table_counts=collections.Counter()
        sample_details=[]
        for gsm,info in sf.items():
            for k in info:
                if 'table_begin' in k:table_counts[gsm]+=1
            supp=info.get('supplementary_file',[])
            for u in supp:
                if u and u!='NONE':file_refs.append({'gsm':gsm,'url':u,'ext':u.split('/')[-1].split('.')[-1]})
            if len(sample_details)<8:sample_details.append({'gsm':gsm,'details':info})
        audit['results']['soft_sample_details_first']=sample_details
        audit['results']['soft_supplement_links_count']=len(file_refs)
        audit['results']['soft_supplement_ext_counts']=dict(collections.Counter(x['ext'] for x in file_refs))
        audit['results']['soft_table_begin_count']=len(table_counts)
        write_csv(out/'gsm_supplement_links.csv',file_refs)
    for sub in ['suppl/','miniml/','matrix/']:
        u=BASE+'/'+sub
        body,info=get(u,limit=1_500_000,timeout=35)
        if body:
            text=body.decode('utf8','replace')
            links=re.findall(r'href=["\']([^"\']+)["\']',text,re.I)
            info['links']=[html.unescape(x) for x in links if x not in ('../','/','?C=N;O=D','?C=M;O=A')][:100]
        audit['results'][sub.rstrip('/')+'_listing']=info
    # Standard GEO Series processed candidates and raw archive, HEAD only.
    for name in [ACC+'_RAW.tar',ACC+'_counts.txt.gz',ACC+'_matrix.csv.gz',ACC+'_processed.csv.gz']:
        audit['results']['head_'+name]=head(BASE+'/suppl/'+name)
    # Public small GSM page preview to locate supplementary resources, no data modification.
    if rows:
        gsm=rows[0]['gsm']
        u='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+gsm
        body,info=get(u,limit=2_000_000,timeout=35)
        if body:
            t=body.decode('utf8','replace')
            info['supp_candidates']=list(dict.fromkeys(re.findall(r'(?:https?://[^"<>\s]+|ftp://[^"<>\s]+|GSM\d+[^"<>\s]*)\.(?:dcc|pkc|txt\.gz|csv\.gz|xlsx|tar)',t,re.I)))[:25]
            info['excerpt_count_mentions']=len(re.findall('read counts',t,re.I))
        audit['results']['first_GSM_page']=info
    audit['decision']='HOLD_PATIENT_ID_AND_QUANTITATIVE_LINKAGE'
    audit['note']='GeoMx regions are nested in patients. Without verified patient IDs and probe-to-gene mapping, do not compute ROI-independent DEGs.'
    (out/'AUDIT.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf8')
    (out/'README.md').write_text('# GSE294009 first probe\n\nDo not call sample prefixes patients. GEO counts 204 registrations; original paper used 156 filtered ROI.\n\nSee AUDIT.json and linkage CSV; all conclusions gated on explicit donor and gene-count mapping.\n')
    print(json.dumps({'n_gsm':audit['results'].get('n_GSM'),'soft':audit['results'].get('soft_n_samples'),'n_supp':audit['results'].get('soft_supplement_links_count'),'raw_tar':audit['results']['head_'+ACC+'_RAW.tar'],'listing':audit['results']['suppl_listing'].get('links',[])[:30], 'decision':audit['decision']},ensure_ascii=False))

if __name__=='__main__':main()