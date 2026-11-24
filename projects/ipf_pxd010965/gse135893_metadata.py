#!/usr/bin/env python3
"""Audit GEO GSE135893 public scRNA metadata for IPF plasma-protein gene mapping.

Fetches only GEO small supplemental files. Never assumes PF=IPF, never performs
cell-level tests using cells as independent donors, never claims protein origin.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import time
import urllib.request

ACC='GSE135893'
TARGETS=['HPR','GPX3','CPN1','FCN3','CFH','C3']
FILES=['GSE135893_IPF_metadata.csv.gz','GSE135893_genes.tsv.gz','GSE135893_barcodes.tsv.gz']

def fetch(name: str, max_bytes=12_000_000):
    urls=[
        f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={ACC}&format=file&file={name}',
        f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/{ACC}/suppl/{name}',
    ]
    errors=[]
    for url in urls:
        for retry in range(2):
            try:
                req=urllib.request.Request(url, headers={'User-Agent':'NCBI-GEO-research-reproduction/1.0','Accept':'application/octet-stream,*/*'})
                with urllib.request.urlopen(req,timeout=60) as response:
                    data=response.read(max_bytes+1)
                    if len(data)>max_bytes: raise ValueError('exceeds small-file size cap')
                    if data[:2]!=b'\x1f\x8b': raise ValueError('response not gzip')
                    return data,{'url':url,'compressed_bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
            except Exception as exc:
                errors.append(f'{type(exc).__name__}: {str(exc)[:120]}')
                time.sleep(0.7*(retry+1))
    raise RuntimeError(f'Could not fetch {name}: '+ '; '.join(errors[-4:]))

def read_gz_table(content,sep=',',has_header=True):
    raw=gzip.decompress(content).decode('utf-8-sig',errors='replace')
    reader=csv.reader(io.StringIO(raw),delimiter=sep)
    head=next(reader)
    rows=[row for row in reader if row]
    if not has_header:
        rows.insert(0,head)
        head=[f'column_{i+1}' for i in range(len(rows[0]))]
    return head,rows

def choose_cols(header):
    keywords=('patient','sample','donor','disease','condition','cell','cluster','annotation','type','group','barcode','origin','batch','subject','diagnosis','ident','status','tissue')
    return [i for i,v in enumerate(header) if any(k in v.lower() for k in keywords)]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='output/ipf_gse135893_metadata');a=ap.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    report={'accession':ACC,'checked_utc':dt.datetime.now(dt.timezone.utc).isoformat(), 'targets':TARGETS,
        'design_note':'PF includes pulmonary fibrosis subtypes; not every PF sample is IPF. Biological replication = donor.',
        'source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893','files':{},'decision':'PENDING'}
    for name in FILES:
        data,info=fetch(name)
        report['files'][name]=info
        sep=',' if name.endswith('.csv.gz') else '\t'
        hdr,rows=read_gz_table(data,sep,has_header=('metadata' in name))
        report['files'][name]['table']={'n_rows':len(rows),'n_cols':len(hdr),'header':hdr[:120],
                                      'first_row_sample':rows[0][:12] if rows else []}
        if 'metadata' in name:
            report['metadata_columns']=hdr
            report['metadata_count']=len(rows)
            interesting={}
            for i in choose_cols(hdr):
                vals=[r[i] if len(r)>i else '' for r in rows]
                count=Counter(vals)
                interesting[hdr[i]]={'unique':len(count),'top':count.most_common(15),
                                     'empty':count.get('',0)}
            report['metadata_column_profiles']=interesting
            report['first_column_n_unique']=len(set(r[0] for r in rows if r))
        elif 'genes' in name:
            genes=[]
            for j,row in enumerate(rows):
                for sym in TARGETS:
                    if any(x.upper()==sym for x in row):
                        genes.append({'row_index_zero_based':j,'gene_fields':row,'target':sym})
            report['gene_matches']=genes
            report['gene_number']=len(rows)
        elif 'barcodes' in name:
            report['barcode_number']=len(rows)
            report['barcode_unique']=len(set(r[0] for r in rows if r))
    ng=report.get('gene_number',0)
    nm=report.get('metadata_count',0)
    nb=report.get('barcode_number',0)
    matches=report.get('gene_matches',[])
    required=set(TARGETS)-{'C3'}
    have=set(g['target'] for g in matches)
    report['decision']=('GO_TARGET_ROW_STREAM' if (nm==nb and nm>0 and ng>10000 and required.issubset(have))
                       else 'HOLD_MAPPING_OR_METADATA_MISMATCH')
    (out/'metadata_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# GSE135893 IPF scRNA cell-origin audit — metadata stage', '',f"- status: **{report['decision']}**",
           f'- metadata rows: {nm}; barcodes: {nb}; gene rows: {ng}',
           f"- target gene matches: {', '.join(sorted(have))}",
           '- Do not use cell counts as independent donors; plasma protein may be of extra-pulmonary origin.',
           '- Next: gene-specific 1GB MTX stream, then donor-aware pseudo-bulk (if annotations permit).', '',
           '## Source SHA256']
    lines += [f"- {k}: `{v['sha256']}`" for k,v in report['files'].items()]
    lines += ['', '## Metadata columns reviewed']
    lines += [f"- {k}: {v['unique']} distinct values; examples {v['top'][:5]}" for k,v in report.get('metadata_column_profiles',{}).items()]
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['decision'],'n_meta':nm,'n_barcodes':nb,'n_genes':ng,'gene_matches':matches,'column_profile':report.get('metadata_column_profiles',{})},ensure_ascii=False,indent=2)[:22000])

if __name__=='__main__': main()