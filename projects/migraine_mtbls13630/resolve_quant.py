#!/usr/bin/env python3
"""MTBLS13630 public-data patient/assay/quantification feasibility audit.

No inference of treatment response or merging participants across studies. Produces
aggregate counts, hashed ISA row IDs, and evidence-backed GO/HOLD statuses.
"""
from __future__ import annotations
import argparse, collections, csv, datetime as dt, hashlib, html.parser, json, re, urllib.error, urllib.parse, urllib.request
from pathlib import Path

STUDY = 'MTBLS13630'
API = f'https://www.ebi.ac.uk/metabolights/ws/studies/{STUDY}'
FTP = f'https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{STUDY}/'
ASSAYS = [
 'a_MTBLS13630_GC-MS___metabolite_profiling.txt',
 'a_MTBLS13630_LC-MS_positive_reverse-phase_metabolite_profiling.txt',
 'a_MTBLS13630_LC-MS_negative_reverse-phase_metabolite_profiling.txt',
]
MAFS = [
 'm_MTBLS13630_GC-MS___metabolite_profiling_v2_maf.tsv',
 'm_MTBLS13630_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv',
 'm_MTBLS13630_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv',
]
MAX_REQUEST=8_000_000
UA='Research-reproducibility-audit/2.0 public metadata; no patient-identifiable output'

def fetch(url: str, max_bytes=MAX_REQUEST):
    try:
        request=urllib.request.Request(url, headers={'User-Agent':UA, 'Accept':'application/json,text/plain,text/html,*/*'})
        with urllib.request.urlopen(request,timeout=35) as resp:
            data=resp.read(max_bytes+1)
            if len(data)>max_bytes: return {'status':'TOO_LARGE','at_least':len(data),'url':url}, None
            return {'status':'OK','bytes':len(data),'mime':resp.headers.get('Content-Type',''),'url':url}, data
    except Exception as exc:
        return {'status':'ERROR','reason':str(exc)[:240],'url':url},None

def jparse(data):
    if data is None: return None
    try:return json.loads(data)
    except Exception:return None

def unwrap(value):
    if value is None: return ''
    if isinstance(value, (str,int,float,bool)):return str(value)
    if isinstance(value,dict):
        for k in ('value','text','label','name','term','data'):
            if k in value and isinstance(value[k],(str,int,float,bool)):
                return str(value[k])
        return json.dumps(value,ensure_ascii=False,sort_keys=True)[:250]
    return json.dumps(value,ensure_ascii=False)[:250]

def as_table(data):
    doc=jparse(data)
    if not isinstance(doc,dict):return [],[],{'error':'NOT_JSON_OBJECT'}
    hdr=doc.get('header',{})
    cols=list(hdr) if isinstance(hdr,dict) else (hdr if isinstance(hdr,list) else [])
    raw=doc.get('data',[])
    if isinstance(raw,dict):raw=raw.get('rows',[])
    if isinstance(raw,dict):raw=list(raw.values())
    if not isinstance(raw,list):return cols,[],{'error':'NO_ROWS','data_type':type(raw).__name__}
    rows=[]
    for item in raw:
        if not isinstance(item,dict):continue
        rows.append({k:unwrap(v) for k,v in item.items() if k!='index'})
    return cols,rows,{'row_count':len(rows),'has_sample_abundance':doc.get('sample_abundance'),'original_row_count':len(raw)}

def counts(rows,key, max_categories=18):
    counter=collections.Counter(r.get(key,'').strip() for r in rows)
    return {'distinct':len(counter),'missing':counter.get('',0),'top':dict(counter.most_common(max_categories))}

def signed_label(rows):
    candidates=[k for k in (rows[0] if rows else {}) if re.search('variant|group|disease|condition|sample type|organism part|factor',k,re.I)]
    return {k:counts(rows,k,20) for k in candidates}

def row_digest(s):return hashlib.sha256(('metabolights-public/'+s).encode('utf-8')).hexdigest()[:16]

class LinkParser(html.parser.HTMLParser):
    def __init__(self):super().__init__();self.hrefs=[]
    def handle_starttag(self,tag,attrs):
        if tag=='a':
            href=dict(attrs).get('href')
            if href:self.hrefs.append(href)

def list_http_dir(url,limit=1_500_000):
    stat, data=fetch(url,limit)
    if not data:return stat,[]
    p=LinkParser();p.feed(data.decode('utf-8','replace'))
    files=[urllib.parse.unquote(x) for x in p.hrefs if not x.startswith('?') and x not in ('../','..','/')]
    return stat,files[:4000]

def write_csv(path,rows,columns):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='output/migraine_mtbls13630_quant');args=parser.parse_args()
    out=Path(args.out);out.mkdir(exist_ok=True,parents=True)
    audit={'study':STUDY,'checked_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
           'questions':['Are true participant groups present in public ISA sheet?',
                        'Can individual analytical assays be linked without pseudoreplication?',
                        'Is a numeric metabolite-by-sample matrix available publicly?'],
           'file_access':{},'assay_tables':{},'maf_tables':{},'ftp':{},'warnings':[]}
    stat,data=fetch(API+'/s_MTBLS13630.txt');audit['file_access']['sample']=stat
    cols, samples, structure=as_table(data)
    audit['sample_metadata']={'structure':structure,'columns':cols,'group_candidates':signed_label(samples),
                              'source_name':counts(samples,'Source Name'),'sample_name':counts(samples,'Sample Name')}
    write_csv(out/'sample_manifest_hashed.csv',
              [{'sample_hash':row_digest(r.get('Sample Name','')),
                'source_hash':row_digest(r.get('Source Name','')),
                'variant':r.get('Characteristics[Variant]',''),
                'sample_type':r.get('Characteristics[Sample type]',''),
                'organism_part':r.get('Characteristics[Organism part]','')} for r in samples],
              ['sample_hash','source_hash','variant','sample_type','organism_part'])
    sample_keys=set(r.get('Sample Name','') for r in samples)
    maf_names=set()
    derived_names=set()
    all_assay=[]
    for name in ASSAYS:
        st,body=fetch(API+'/'+urllib.parse.quote(name));audit['file_access'][name]=st
        cols,rows,shape=as_table(body)
        sample_col='Sample Name'
        matched=sum(1 for r in rows if r.get(sample_col,'') in sample_keys)
        missing_names=sorted(set(r.get(sample_col,'') for r in rows)-sample_keys)
        for r in rows:
            maf_names.add(r.get('Metabolite Assignment File',''))
            for k in ('Derived Spectral Data File','Raw Spectral Data File'):
                if r.get(k,''):derived_names.add(r[k])
            all_assay.append({'assay':name,'sample_hash':row_digest(r.get(sample_col,'')),
                              'matched_to_sample_sheet':r.get(sample_col,'') in sample_keys,
                              'derived_file':r.get('Derived Spectral Data File',''),
                              'maf_file':r.get('Metabolite Assignment File','')})
        audit['assay_tables'][name]={'shape':shape,'matched_sample_rows':matched,'unmatched_sample_name_count':len(missing_names),
                                    'derived':counts(rows,'Derived Spectral Data File',12),
                                    'assignment':counts(rows,'Metabolite Assignment File',12),
                                    'raw_files':counts(rows,'Raw Spectral Data File',3)}
    write_csv(out/'assay_link_manifest.csv',all_assay,['assay','sample_hash','matched_to_sample_sheet','derived_file','maf_file'])
    numeric_by_maf={}
    for name in MAFS:
        st,body=fetch(API+'/'+urllib.parse.quote(name));audit['file_access'][name]=st
        cols,rows,shape=as_table(body)
        columns=sorted(set(k for r in rows for k in r.keys()))
        abund_columns=[x for x in columns if re.search(r'abundance|concentration|peak.area|intensity',x,re.I)]
        abund_numeric={}
        for k in abund_columns:
            vals=[r.get(k,'') for r in rows]
            numeric=sum(bool(re.fullmatch(r'-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?',v.strip())) for v in vals)
            abund_numeric[k]={'numbers':numeric,'nonempty':sum(bool(v.strip()) for v in vals),'rows':len(rows)}
        audit['maf_tables'][name]={'shape':shape,'header_count':len(cols),'keys_count':len(columns),
                                   'abundance_columns':abund_numeric,'unknown_columns':len(set(columns)-set(cols)),
                                   'abundance_flag':shape.get('has_sample_abundance')}
        numeric_by_maf[name]=sum(v['numbers'] for v in abund_numeric.values())
    # Does the explicit API query expose extra sample-wise abundance columns?
    maf_query=MAFS[0]
    st,data=fetch(API+'/'+urllib.parse.quote(maf_query)+'?sample_abundance=true')
    cols,rows,shape=as_table(data)
    audit['maf_sample_abundance_query']={'response':st,'structure':shape,'col_count':len(cols),
                                        'comparison_to_default_col_count':audit['maf_tables'][maf_query]['header_count']}
    # Public repository listing (no large instrument data download).
    for url in (FTP, 'https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/'+STUDY+'/'):
        if url in audit['ftp']:continue
        stat,names=list_http_dir(url)
        audit['ftp'][url]={'status':stat,'n_listing_links':len(names),
            'candidate_processed_files':[n for n in names if re.search(r'processed|result|matrix|peak|\.csv$|\.xlsx$|\.tsv$|\.txt$',n,re.I)][:100],
            'files_total':len(names)}
        if names:
            (out/'public_directory_file_names.txt').write_text('\n'.join(names)+'\n',encoding='utf-8')
    audit['linked_data']={'sample_metadata_rows':len(samples),'distinct_sample_names':len(sample_keys),
        'distinct_source_names':len(set(r.get('Source Name','') for r in samples)),
        'assay_rows':len(all_assay),'assay_matched_to_sample_rows':sum(x['matched_to_sample_sheet'] for x in all_assay),
        'raw_and_derived_n_distinct_filenames':len(derived_names),'maf_n_distinct_filenames':len(maf_names),
        'numeric_abundance_cells_via_metadata_api':numeric_by_maf,
        'possible_repeated_subjects':len(samples)>len(set(r.get('Source Name','') for r in samples))}
    # Only release automated case-control calculations when participant identity,
    # phenotype, and numerical specimen matrix are all explicitly linked.
    numerical=any(v>0 for v in numeric_by_maf.values())
    group_columns=[k for k,v in audit['sample_metadata']['group_candidates'].items() if v['distinct']>=2 and v['missing']<len(samples)]
    if numerical and group_columns and len(sample_keys)==len(samples) and all(x['matched_to_sample_sheet'] for x in all_assay):
        audit['decision']='REVIEW_NUMERIC_MATRIX_PRESENT_NEEDS_QC'
    elif not numerical:
        audit['decision']='HOLD_NO_LINKED_PER_SAMPLE_QUANT_MATRIX'
    else:
        audit['decision']='HOLD_SUBJECT_OR_GROUP_MAPPING'
    audit['warnings']=[
      'Sample records can be repeated venous/arterial specimens. Rows are not independent participants.',
      'PFO-positive migraine vs healthy controls is NOT preventative-treatment response.',
      'Do not infer clinic group from sample identifiers unless explicitly defined by ISA metadata.',
      'MAF rows represent annotations and cannot be substituted for patient-level quantitative measurements.',
      'Do not claim reanalysis of 4,858 features solely from 211 annotation rows.'
    ]
    (out/'audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    lines=['# MTBLS13630 — sample/assay/quantitative matrix reusability audit','',
       f"Decision: **{audit['decision']}**",'',
       f"- Sample rows: {len(samples)}, distinct sample names: {len(sample_keys)}; source names: {audit['linked_data']['distinct_source_names']}",
       f"- Assay rows: {len(all_assay)}, mapped to sample: {audit['linked_data']['assay_matched_to_sample_rows']}",
       f"- MAF numerical abundance cell counts: {json.dumps(numeric_by_maf)}",
       '- External data scope: PFO-associated migraine (doi:10.1016/j.isci.2026.115587), not Yonsei preventive therapy.',
       '- Next: verify large *processed* files via FTP and patient-level phenotype mapping; absent this, stop at feasibility and switch to IPF.',
       '', '## Metadata group and sample types']
    for k,v in audit['sample_metadata']['group_candidates'].items():
        lines.append(f"- {k}: {v['top']}")
    lines.extend(['','## Data dictionary / metadata file references','- https://ebi-metabolights.github.io/guides/Files/',
       '- https://ebi-metabolights.github.io/guides/MAF/',
       '- https://www.ebi.ac.uk/metabolights/MTBLS13630',''])
    (out/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({'decision':audit['decision'],'linked':audit['linked_data'],
      'group_candidates':audit['sample_metadata']['group_candidates'],
      'ftp':audit['ftp']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()