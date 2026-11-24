#!/usr/bin/env python3
"""Retrieve one small public processed Proteins.csv and audit quantitative reuse.

No raw MS download. No assumptions about IPF labels, no clinical predictions.
"""
import argparse, csv, datetime, hashlib, io, json, pathlib, re, urllib.request
BASE='https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD010965'
CANDIDATE='Proteins.csv'
MAX=3_000_000

def fetch(url):
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'PXD010965-Open-Data-Audit/1.0'})
        with urllib.request.urlopen(req,timeout=25) as r:
            data=r.read(MAX+1)
            if len(data)>MAX:return {'status':'TOO_LARGE','size_at_least':len(data),'url':url},None
            return {'status':'OK','bytes':len(data),'content_type':r.headers.get('content-type',''),'url':url},data
    except Exception as e:return {'status':'ERROR','reason':str(e)[:230],'url':url},None

def obj_files(o):
    if isinstance(o,list):return [x for x in o if isinstance(x,dict)]
    if isinstance(o,dict):
        for k in ('files','projectFiles','_embedded','content'):
            if k in o:
                p=obj_files(o[k]);
                if p:return p
    return []

def links(obj):
    results=[]
    def walk(o,depth=0):
        if depth>5:return
        if isinstance(o,dict):
            for k,v in o.items():
                if isinstance(v,str) and v.startswith(('https://','http://','ftp://')):results.append({'field':k,'url':v})
                elif isinstance(v,(dict,list)):walk(v,depth+1)
        elif isinstance(o,list):
            for v in o:walk(v,depth+1)
    walk(obj)
    return results

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='output/ipf_pxd010965_quant');args=p.parse_args()
    out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
    res={'accession':'PXD010965','checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'target':'Proteins.csv','status':'HOLD_NO_QUANT','attempts':[]}
    st,body=fetch(BASE+'/files?pageSize=200');res['file_api']=st
    entries=obj_files(json.loads(body)) if body else []
    f=next((x for x in entries if x.get('fileName')==CANDIDATE),None)
    if f:
        res['candidate_info']={k:f.get(k) for k in ('fileName','fileSizeBytes','fileCategory','checksum','fileSize','fileSizeInBytes') if k in f}
        res['file_url_fields']=[{'field':x['field'],'host':x['url'].split('/')[2]} for x in links(f)]
    discovered=[x['url'] for x in links(f)] if f else []
    # Try exact advertised project location first, then the canonical public archive path.
    fallback='https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/05/PXD010965/Proteins.csv'
    candidates=[x for x in discovered if re.search(r'\.csv|Proteins',x,re.I)] + [fallback]
    seen=set()
    for url in candidates[:6]:
        if url in seen:continue
        seen.add(url)
        status,data=fetch(url);res['attempts'].append(status)
        if not data or b'<'==data[:1]:continue
        try:text=data.decode('utf-8-sig')
        except Exception:continue
        if not text.count('\n')>=3:continue
        try:dia=csv.Sniffer().sniff(text[:4000],delimiters=',;\t').delimiter
        except csv.Error:dia=','
        reader=csv.reader(io.StringIO(text),delimiter=dia)
        hdr=next(reader,[]);rows=list(reader)
        if len(hdr)<2 or len(rows)<2:continue
        nrows=len(rows);ncols=len(hdr)
        counts=[]
        for j in range(ncols):
            val=[row[j].strip() for row in rows if len(row)>j]
            cnt=sum(bool(re.fullmatch(r'-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',v)) for v in val)
            if cnt>=max(5,int(nrows*.5)):counts.append(j)
        sample_like=[(i,name) for i,name in enumerate(hdr) if re.search(r'IPF|control|healthy|^\d+[A-Za-z]?|^Abundance',name,re.I)]
        res['processed_table']={'rows':nrows,'columns':ncols,'numeric_majority_columns':len(counts),
                                'column_headers':hdr[:80],
                                'sample_like_column_count':len(sample_like),
                                'delimiter':repr(dia),'sha256':hashlib.sha256(data).hexdigest(),
                                'source_url':url}
        (out/'Proteins_public_source.csv').write_bytes(data)
        res['status']='GO_NUMERIC_MATRIX_HOLD_GROUP_MAPPING' if len(counts)>=4 else 'HOLD_PROTEINS_TABLE_HAS_NO_SAMPLE_MATRIX'
        break
    if 'processed_table' not in res:res['status']='HOLD_PROCESSED_CSV_NOT_ACCESSIBLE'
    res['warnings']=['A 37-person Finnish 2020 cross-sectional discovery study cannot validate longitudinal Korean FVC progression.',
      'Do not group sample column names as IPF vs control without provenance for raw files and subject/group coding.',
      'Even a populated per-sample protein matrix does not independently verify label-free normalization and batch controls.']
    (out/'audit.json').write_text(json.dumps(res,indent=2,ensure_ascii=False),encoding='utf-8')
    (out/'REPORT.md').write_text('\n'.join(['# IPF public Proteins.csv feasibility','',
       f"Decision: **{res['status']}**",f"- API file index: {st['status']}",
       f"- Processed source accessible: {'YES' if 'processed_table' in res else 'NO'}",
       f"- Matrix shape: {res.get('processed_table',{}).get('rows','?')} rows × {res.get('processed_table',{}).get('columns','?')} columns",
       f"- Numeric-majority columns: {res.get('processed_table',{}).get('numeric_majority_columns','?')}",
       '- No clinical progression is available in this public assay and no label-based statistical comparison has been performed.',
       '- Dataset https://doi.org/10.6019/PXD010965', '']),encoding='utf-8')
    print(json.dumps({'status':res['status'],'processed_table':res.get('processed_table'),'attempts':res['attempts']},ensure_ascii=False,indent=2))
if __name__=='__main__':main()