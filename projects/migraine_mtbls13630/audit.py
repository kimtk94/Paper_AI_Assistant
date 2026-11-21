#!/usr/bin/env python3
"""Public-access metadata and small processed-table audit for MTBLS13630.

Limits: PFO-associated migraine (30 cases, 17 controls), not preventive-treatment
response, not a replication of Yonsei serum SCFA measurements.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import io
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ACCESSION = 'MTBLS13630'
BASE = 'https://www.ebi.ac.uk/metabolights/ws'
MAX_BYTES = 4_000_000
ENDPOINTS = {
  'study': f'{BASE}/studies/{ACCESSION}',
  'files': f'{BASE}/studies/{ACCESSION}/files',
  'sample': f'{BASE}/studies/{ACCESSION}/s_{ACCESSION}.txt',
}
USER_AGENT = 'Novogene-Research-Feasibility/1.0 (public ISA metadata only)'

def get(url: str, limit: int = MAX_BYTES):
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT, 'Accept': 'application/json,text/plain,*/*'})
    try:
        with urllib.request.urlopen(req, timeout=25) as res:
            body = res.read(limit + 1)
            if len(body) > limit:
                return {'status':'TOO_LARGE','http_status':res.status,'bytes_gt':limit,'url':url}, None
            return {'status':'OK','http_status':res.status,'bytes':len(body),'content_type':res.headers.get('Content-Type',''),'url':url}, body
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        return {'status':'ERROR','reason':str(exc)[:280], 'url':url}, None

def maybe_json(body):
    try: return json.loads(body)
    except (json.JSONDecodeError, UnicodeError, TypeError): return None

def walk(obj, prefix='', depth=0):
    if depth>8: return
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = prefix + '/' + str(k)
            if isinstance(v, (dict,list)):
                yield from walk(v,p,depth+1)
            elif isinstance(v,str) and re.search(r'\.(tsv|txt|csv|xlsx|mzml|mzxml|raw|zip)(?:\?|$)',v,re.I):
                yield p, v
    elif isinstance(obj, list):
        for i, v in enumerate(obj[:2000]):
            yield from walk(v,f'{prefix}/{i}',depth+1)

def table_audit(raw:bytes):
    if len(raw)>MAX_BYTES: return {'status':'TOO_LARGE'}
    text=raw.decode('utf-8-sig',errors='replace')
    first=text.splitlines()[:4]
    if not first:return {'status':'EMPTY'}
    try:
        dialect=csv.Sniffer().sniff('\n'.join(first),delimiters='\t,;')
        delimiter=dialect.delimiter
    except csv.Error:
        delimiter='\t'
    reader=csv.reader(io.StringIO(text),delimiter=delimiter)
    rows=[]
    for row in reader:
        rows.append(row)
        if len(rows)>=10000:break
    return {'status':'PARSED','n_preview_rows':max(0,len(rows)-1), 'n_columns':len(rows[0]) if rows else 0,
            'column_names':rows[0][:40] if rows else [],'delimiter':repr(delimiter), 'truncated_rows':len(rows)>=10000}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',default='output/migraine_mtbls13630')
    args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    result={'study':ACCESSION,'checked_at_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
        'primary_research_question':'Are public migraine-related LC-MS matrices reusable for pathway hypothesis generation?',
        'study_design_from_article':{'disease':'migraine with patent foramen ovale', 'cases':30,'healthy_controls':17,
            'validation_scope':'cross-sectional PFO subtype, NOT preventive-treatment response',
            'doi':'10.1016/j.isci.2026.115587'},
        'endpoints':{},'candidate_metadata_files':[], 'processed_data_validation':'NOT_TESTED'}
    payloads={}
    for name,url in ENDPOINTS.items():
        status,body=get(url)
        result['endpoints'][name]=status
        if body is not None: payloads[name]=body
        if name in ('study','files') and body is not None:
            doc=maybe_json(body)
            if doc is not None:
                result['endpoints'][name]['top_keys']=list(doc)[:35] if isinstance(doc,dict) else []
                result['candidate_metadata_files'] += [{'path':k,'value':v} for k,v in list(walk(doc))[:100]]
        if name=='sample' and body:
            result['sample_table_audit']=table_audit(body)
    result['candidate_metadata_files']=result['candidate_metadata_files'][:100]
    # A positive metadata fetch does NOT establish public access to an analyzable matrix.
    s=result['endpoints']['study']['status']
    if s=='OK' and result.get('sample_table_audit',{}).get('status')=='PARSED':
        result['decision']='HOLD_FOR_ASSAY_MATRIX_AND_GROUP_MAPPING'
    elif s=='OK':
        result['decision']='HOLD_SAMPLE_METADATA_ACCESS'
    else:
        result['decision']='HOLD_API_ACCESS'
    result['warnings']=[
        'Do not use PFO migraine metabolite associations to claim treatment effect or match them to Yonsei SCFA patients.',
        'Do not claim any individual metabolite replicated without quantitative assay matrix, IDs, measurement units, QC, and patient mapping.',
        'The 2026 UK Biobank migraine metabolomics/proteomics studies reduce novelty for general case-control discovery.',
        'If quantitative matrices become accessible, perform per-batch QC and covariate-adjusted, patient-level analysis; never mix arterial and venous paired samples as independent people.'
    ]
    (out/'audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    tab=result.get('sample_table_audit',{})
    report=['# MTBLS13630 공개 데이터 접근성 검증', '',f"- 상태: **{result['decision']}**",
        '- 연구: 편두통+PFO 30명 및 건강대조군 17명 (원논문 기준)',
        '- 연구용도: 편두통 대사경로 기전 가설, 예방치료 반응 예측의 직접 검증 아님.',
        f"- study API: {result['endpoints']['study']['status']}",
        f"- files API: {result['endpoints']['files']['status']}",
        f"- sample metadata: {result['endpoints']['sample']['status']}",
        f"- sample table preview: {tab.get('n_preview_rows','NA')} rows, {tab.get('n_columns','NA')} cols",
        '- 수치 대사체 행렬 확보·환자 매핑·분산검증: 아직 통과하지 않음',
        '','## 주요 공개 자료',
        '- 원논문: https://doi.org/10.1016/j.isci.2026.115587',
        '- MetaboLights: https://www.ebi.ac.uk/metabolights/MTBLS13630',
        '- 대비 연구: https://doi.org/10.1186/s10194-026-02268-4 (세브란스, 혈청 3종 SCFA)',
        '','## 다음 검증','1. metadata/assay 파일과 sample mapping 확보',
        '2. Quantitative processed LC-MS/GC-MS feature matrix 입수',
        '3. QC 및 PFO, sex, arterial/venous 측정 조건 분리',
        '4. Indoleacrylic acid / tryptophan / SCFA 기전 가설 후보 비교',
        '5. 최종 실험은 반복채혈과 약물·증상 기록으로 새로 설계', '']
    (out/'REPORT.md').write_text('\n'.join(report),encoding='utf-8')
    print(json.dumps({'decision':result['decision'],'endpoints':result['endpoints'],'sample_table':result.get('sample_table_audit')},ensure_ascii=False,indent=2))

if __name__=='__main__':main()