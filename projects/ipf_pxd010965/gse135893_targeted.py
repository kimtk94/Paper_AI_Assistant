#!/usr/bin/env python3
"""Stream GSE135893 raw MTX and extract six selected genes only.

Uses exact GEO barcode intersection with 114,396 annotated cells. Performs
sample-aware descriptive summaries comparing IPF vs Control, never cell-level
pseudo-replicated P values, and never claims plasma protein secretory origin.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
import time
import urllib.request

from gse135893_metadata import fetch,read_gz_table,ACC,TARGETS
MATRIX_NAME='GSE135893_matrix.mtx.gz'

def open_mtx():
    links=[f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={ACC}&format=file&file={MATRIX_NAME}',
           f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/{ACC}/suppl/{MATRIX_NAME}']
    errors=[]
    for url in links:
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'NCBI-GEO-targeted-gene-research/1.0','Accept':'application/octet-stream,*/*'})
            resp=urllib.request.urlopen(req,timeout=180)
            if resp.headers.get('Content-Type','').startswith('text/html'):
                raise ValueError('HTML response, not MTX gzip')
            return resp,url
        except Exception as exc:
            errors.append(f'{url}: {type(exc).__name__}: {str(exc)[:150]}')
            time.sleep(1)
    raise RuntimeError('MTX access failure: '+'; '.join(errors))

def norm_num(x, fallback=0.0):
    try:return float(x)
    except (ValueError,TypeError):return fallback

def median(x): return statistics.median(x) if x else None

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--out',default='output/ipf_gse135893_targeted');a=pa.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    result={'accession':ACC,'requested_targets':TARGETS,
        'checked_at_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
        'decision':'HOLD_UNVERIFIED','scientific_constraints':[
        '30 distinct donors across PF spectrum: only Diagnosis==IPF and Diagnosis==Control are compared',
        'Sample_Name collapses multiple orig.ident sublibraries to one donor',
        'Unit of inference is donor (not 114,396 individual cells)',
        'Plasma protein could originate in liver or immune system and tissue mRNA is not source proof',
        'Reported values are descriptive (no cell-level significance test)']}
    raw_meta,minfo=fetch('GSE135893_IPF_metadata.csv.gz')
    hdr,rows=read_gz_table(raw_meta,',',has_header=True)
    raw_barcodes,binfo=fetch('GSE135893_barcodes.tsv.gz')
    _,barcodes=read_gz_table(raw_barcodes,'\t',has_header=False)
    raw_genes,ginfo=fetch('GSE135893_genes.tsv.gz')
    _,genes=read_gz_table(raw_genes,'\t',has_header=False)
    idx={k:i for i,k in enumerate(hdr)}
    for k in ('Diagnosis','Sample_Name','celltype','nCount_RNA'):
        if k not in idx: raise RuntimeError('missing exact metadata field '+k)
    meta_by_barcode={r[0]:r for r in rows}
    if len(meta_by_barcode)!=len(rows):raise RuntimeError('duplicate annotated metadata barcode')
    barcode_to_matrix_index={r[0]:i for i,r in enumerate(barcodes) if r}
    n_match=sum(b in barcode_to_matrix_index for b in meta_by_barcode)
    result.update({'metadata_cells':len(rows),'barcodes_in_matrix':len(barcodes),
                   'barcode_intersection':n_match,'genes_n':len(genes),
                   'n_donors_IPF':len({r[idx['Sample_Name']] for r in rows if r[idx['Diagnosis']]=='IPF'}),
                   'n_donors_Control':len({r[idx['Sample_Name']] for r in rows if r[idx['Diagnosis']]=='Control'}),
                   'metadata_sha256':minfo['sha256'],'barcode_sha256':binfo['sha256'],'gene_sha256':ginfo['sha256']})
    if n_match != len(rows):
        result['decision']='HOLD_BARCODE_INTERSECTION';(out/'summary.json').write_text(json.dumps(result,indent=2));print(result);return
    gene_to_index={g[0].upper():j+1 for j,g in enumerate(genes) if len(g)>0}
    targets={gene_to_index[x.upper()]:x for x in TARGETS if x.upper() in gene_to_index}
    if len(targets)<5:raise RuntimeError(f'Only mapped {len(targets)} candidates')
    result['gene_row_indices_1_based']={v:k for k,v in targets.items()}
    # Restrict statistical scope to IPF and Control, exclude other interstitial lung disease diagnoses.
    meta_for_col={barcode_to_matrix_index[b]+1:{'group':r[idx['Diagnosis']],
        'donor':r[idx['Sample_Name']],'celltype':r[idx['celltype']],
        'nUMI':norm_num(r[idx['nCount_RNA']])} for b,r in meta_by_barcode.items()
        if r[idx['Diagnosis']] in ('IPF','Control')}
    groups=defaultdict(lambda:{'cells':0,'total_umi':0.0})
    for m in meta_for_col.values():
        key=(m['group'],m['donor'],m['celltype'])
        groups[key]['cells']+=1
        groups[key]['total_umi']+=m['nUMI']
    gstats=defaultdict(lambda:{'umi_count':0,'n_positive_cells':0,'sum_log1p_cpm10k':0.0})
    total_mtx_entries=0
    target_entries=0
    with_mtx=False
    resp=None
    try:
        resp,url=open_mtx();result['matrix_url']=url
        with gzip.GzipFile(fileobj=resp,mode='rb') as gz:
            # parse MatrixMarket header and dimensions
            first=gz.readline().decode('utf-8-sig').strip()
            if not first.startswith('%%MatrixMarket matrix coordinate'):
                raise RuntimeError('not coordinate MatrixMarket: '+first[:100])
            dims=None
            while True:
                line=gz.readline()
                if not line:break
                if line.startswith(b'%'):continue
                dims=[int(v) for v in line.split()]
                break
            if not dims or len(dims)!=3: raise RuntimeError('no MatrixMarket dimensions')
            result['mtx_dims']={'genes':dims[0],'barcodes':dims[1],'nonzero_entries':dims[2]}
            if dims[0]!=len(genes) or dims[1]!=len(barcodes):
                raise RuntimeError(f'matrix shape {dims[:2]} != gene/barcode lists {(len(genes),len(barcodes))}')
            # Stream ~1GB compressed MTX, never persist the full file or use dense expression matrix.
            for line in gz:
                if not line or line.startswith(b'%'):continue
                total_mtx_entries+=1
                # First field is one-based gene row.
                field=line.split(b' ',1)[0]
                try:row_index=int(field)
                except ValueError:
                    bits=line.split()
                    if len(bits)<3:continue
                    row_index=int(bits[0]); col_index=int(bits[1]);val=float(bits[2])
                else:
                    if row_index not in targets:continue
                    bits=line.split()
                    col_index=int(bits[1]);val=float(bits[2])
                if row_index not in targets:continue
                target_entries+=1
                m=meta_for_col.get(col_index)
                if m is None:continue
                key=(m['group'],m['donor'],m['celltype'],targets[row_index])
                t=gstats[key]
                t['umi_count']+=val
                if val>0:
                    t['n_positive_cells']+=1
                    t['sum_log1p_cpm10k']+=math.log1p(val*10000/(m['nUMI'] or 1.0))
    finally:
        if resp is not None:resp.close()
    result['mtx_parsed_entries']=total_mtx_entries
    result['mtx_target_entries']=target_entries
    if total_mtx_entries!=result['mtx_dims']['nonzero_entries']:
        result['decision']='HOLD_MTX_TRUNCATED'
    else:result['decision']='GO_DESCRIPTIVE_DONOR_ANALYSIS'
    # Donor x cell type table, one row per (donor, celltype, gene).
    donorrows=[]
    for (diag,donor,ct),s in sorted(groups.items()):
        if s['cells']<15:continue
        for gene in TARGETS:
            g=gstats.get((diag,donor,ct,gene),{})
            raw=g.get('umi_count',0)
            pos=g.get('n_positive_cells',0)
            donorrows.append({'diagnosis':diag,'donor':donor,'celltype':ct,'gene':gene,
                             'n_cells':s['cells'],'n_positive_cells':pos,
                             'detection_pct':round(pos/s['cells']*100,5),
                             'pseudobulk_cpm':round(raw/(s['total_umi'] or 1.0)*1e6,7),
                             'mean_log1p_norm':round(g.get('sum_log1p_cpm10k',0)/s['cells'],7)})
    with (out/'donor_celltype_gene.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(donorrows[0]) if donorrows else ['diagnosis','donor','celltype','gene']);writer.writeheader();writer.writerows(donorrows)
    per_ct=defaultdict(lambda:defaultdict(list))
    for r in donorrows:
        per_ct[(r['celltype'],r['gene'])][r['diagnosis']].append(r)
    results=[]
    for (ct,gene),by_diag in sorted(per_ct.items()):
        ipf=by_diag.get('IPF',[]);ctrl=by_diag.get('Control',[])
        if not ipf and not ctrl:continue
        def med(field,xx):return median([r[field] for r in xx])
        ipf_expr=med('pseudobulk_cpm',ipf);ctrl_expr=med('pseudobulk_cpm',ctrl)
        # Descriptive donor medians, not a false cell-level significance test.
        log2fc=(math.log2(1+ipf_expr)-math.log2(1+ctrl_expr)
                if ipf_expr is not None and ctrl_expr is not None else None)
        results.append({'celltype':ct,'gene':gene,'ipf_donors_n':len(ipf),'control_donors_n':len(ctrl),
                       'ipf_median_cpm':ipf_expr,'control_median_cpm':ctrl_expr,
                       'ipf_median_detect_pct':med('detection_pct',ipf),'control_median_detect_pct':med('detection_pct',ctrl),
                       'descriptive_log2fc_pseudobulk_cpm_plus_1':None if log2fc is None else round(log2fc,5),
                       'comparison_status':'DESCRIPTIVE' if len(ipf)>=3 and len(ctrl)>=3 else 'INSUFFICIENT_DONORS'})
    with (out/'celltype_gene_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]) if results else ['celltype','gene']);writer.writeheader();writer.writerows(results)
    result['n_donor_celltype_gene_rows']=len(donorrows);result['n_celltype_gene_summary_rows']=len(results)
    result['n_celltypes']=len({r['celltype'] for r in results});result['targets_mapped']=sorted(targets.values())
    (out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# IPF plasma proteins × GSE135893 lung single cell (targeted MTX stream)', '',
           f"- **Status:** {result['decision']}",
           f"- Matched annotated cells: {n_match:,}/{len(rows):,}; gene rows: {len(genes):,}",
           f"- IPF donors: {result['n_donors_IPF']}; Control donors: {result['n_donors_Control']}",
           f"- Matrix: {result.get('mtx_dims')}; parsed entries: {total_mtx_entries:,}",
           '- All other diagnoses excluded from IPF/control comparative summaries.',
           '- Data units: donor-celltype pseudo-bulk CPM (relative RNA), percent of expressing cells.',
           '- Plasma proteomics is a DIFFERENT cohort; no direct patient matched RNA-protein correlation.',
           '- HPR/CFH/FCN3/CPN1 may be made in liver; absence of lung expression does not invalidate plasma markers.',
           '- No donor-aware differential expression P-values claimed in this targeted descriptive analysis.',
           '', '## Key outputs', '- donor_celltype_gene.csv: group-level summaries, never independent cell P values',
           '- celltype_gene_summary.csv: diagnosis comparison, donor counts and descriptive directions',
           '- summary.json: exact matrix dimensions and mapping provenance', '',
           '## Next research design', '1. Candidate protein PRM in new, baseline Korean IPF plasma specimens.',
           '2. 12-month FVC decline as the preregistered outcome, adjust for age, sex, baseline FVC, CT extent, antifibrotics.',
           '3. Use scRNA as mechanism/context only; liver and circulating immunity cannot be resolved from lung transcript alone.',
           '']
    (out/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__': main()