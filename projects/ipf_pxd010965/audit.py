#!/usr/bin/env python3
"""Audit whether public PRIDE PXD010965 IPF plasma files are reusable.
Metadata only: no raw patient-data download, no disease-effect inference.
"""
import argparse, csv, datetime, json, pathlib, re, urllib.request
ROOT='https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD010965'
BASE='https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD010965'

def get(url,maxbytes=2_000_000):
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'PublicPRIDE-Feasibility/1.0','Accept':'application/json,*/*'}),timeout=30) as r:
            raw=r.read(maxbytes+1)
            if len(raw)>maxbytes:return {'status':'TOO_LARGE','bytes':len(raw)},None
            try:obj=json.loads(raw)
            except ValueError:return {'status':'NON_JSON','bytes':len(raw)},None
            return {'status':'OK','bytes':len(raw),'keys':list(obj)[:30] if isinstance(obj,dict) else [],'url':url},obj
    except Exception as e:return {'status':'ERROR','reason':str(e)[:220],'url':url},None

def collect_files(blob):
    if isinstance(blob,list):return blob
    if isinstance(blob,dict):
        for key in ('files','projectFiles','fileList','content','_embedded','_items','data'):
            if key in blob:
                nested=collect_files(blob[key])
                if nested:return nested
    return []

def scalar(obj,*keys):
    for key in keys:
        if key in obj and obj[key] is not None:
            v=obj[key]
            if isinstance(v,dict):
                return str(v.get('value',v.get('name','')))
            return str(v)
    return ''

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='output/ipf_pxd010965');a=parser.parse_args()
    out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
    audit={'accession':'PXD010965','checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'design':'Public Finnish 2020 baseline IPF plasma protein discovery, NOT Severance longitudinal FVC outcomes',
           'endpoints':{},'file_manifest':[]}
    for ver,base in [('v2',ROOT),('v3',BASE)]:
        st,pr=get(base);audit['endpoints'][ver+'_project']=st
        if pr and isinstance(pr,dict):
            audit[ver+'_metadata']={k:pr.get(k) for k in ('accession','title','publicationDate','projectDescription','sampleProcessingProtocol','dataProcessingProtocol','totalFiles') if k in pr}
        st,raw=get(base+'/files?pageSize=200');audit['endpoints'][ver+'_files']=st
        files=collect_files(raw)
        for f in files:
            if not isinstance(f,dict):continue
            name=scalar(f,'fileName','name','file_name')
            category=scalar(f,'fileCategory','category','fileType')
            size=scalar(f,'fileSizeBytes','fileSize','sizeInBytes','fileSizeInBytes')
            checksum=scalar(f,'checksum','md5','fileChecksum')
            audit['file_manifest'].append({'api_version':ver,'filename':name,'category':category,
                                          'size_bytes':size,'checksum':checksum,
                                          'likely_quant_candidate':bool(re.search(r'protein|quant|intens|abundan|result|xlsx|\.tsv$|\.csv$',name,re.I) and not re.search(r'\.raw|\.mzml$',name,re.I))})
    seen={}
    for f in audit['file_manifest']:
        seen[(f['filename'],f['category'])]=f
    files=list(seen.values())
    with (out/'file_manifest.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['api_version','filename','category','size_bytes','checksum','likely_quant_candidate']);writer.writeheader();writer.writerows(files)
    audit['files_summary']={'num_visible_files':len(files),'quant_candidates':[x['filename'] for x in files if x['likely_quant_candidate']][:40],
                            'raw_files':sum(bool(re.search(r'\.(raw|mzml)(?:\.|$)',x['filename'],re.I)) for x in files),
                            'category_counts':{x:sum(z['category']==x for z in files) for x in set(z['category'] for z in files)}}
    audit['decision']='HOLD_DOWNLOAD_AND_SAMPLE_MAPPING' if files else 'HOLD_PRIDE_FILE_API_OR_NO_LIST'
    audit['warnings']=['A PRIDE file named RESULT is not necessarily a sample-by-protein matrix.',
       'No links between Finnish 2020 donor samples and Korean Severance progression metadata.',
       'Disease-vs-control proteins cannot be called 12-month FVC progression biomarkers without timed clinical follow-up.',
       'Do not claim direct plasma protein-to-lung-cell links without protein-gene ID harmonization and proper cell enrichment controls.']
    (out/'audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'REPORT.md').write_text('\n'.join(['# IPF PXD010965 public-data feasibility', '',
      'CI execution status is separate from the scientific readiness gate.',
      '**'+audit['decision']+'**',
      f"- PRIDE file count: {len(files)}",f"- Candidate result files: {len(audit['files_summary']['quant_candidates'])}",
      '- Next: only small processed results with explicit sample columns and metadata should be downloaded; never automate large RAW reprocessing as a prerequisite for contacting Korean PIs.',
      '- Context: PXD010965 Finnish IPF vs control plasma, doi:10.1038/s41598-020-64759-x; scRNA GSE135893 different unpaired lung donors.',
      '- Next Korean proposal: baseline plasma with 12-month FVC followup; need documented availability and owner authorization.',
      '','## Files','\n'.join('- '+x['filename']+' ['+x['category']+']' for x in files[:60])]),encoding='utf-8')
    print(json.dumps({'decision':audit['decision'],'summary':audit['files_summary'],'endpoints':audit['endpoints']},ensure_ascii=False,indent=2))
if __name__=='__main__':main()