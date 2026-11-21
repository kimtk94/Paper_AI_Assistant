#!/usr/bin/env python3
"""Read only public GSE136831 cell metadata/gene/barcode files, no 2GB MTX.

This is an access + matching audit, NOT a gene-expression or external validation.
"""
import collections, csv, datetime, gzip, hashlib, io, json, re, time, urllib.request
from pathlib import Path

GEO='GSE136831'
NAMES={
 'metadata':'GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz',
 'genes':'GSE136831_AllCells.GeneIDs.txt.gz',
 'barcodes':'GSE136831_AllCells.cellBarcodes.txt.gz',
}
TARGETS=['FCN3','GPX3','CFH','HPR','CPN1','C3']
URL='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136831&format=file&file='

def grab(name):
    url=URL+name
    last=''
    for retry in range(3):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'GEO-Research-Replicability/1.0','Accept':'*/*'})
            with urllib.request.urlopen(req,timeout=70) as f:
                content=f.read(25_000_001)
                if len(content)>25_000_000:raise RuntimeError('Max allowed download size exceeded')
            if content[:2]!=b'\x1f\x8b':raise ValueError('Content not gzipped - GEO may return an HTML redirect/error')
            return {'url':url,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()},gzip.decompress(content).decode('utf-8-sig','replace')
        except Exception as ex:
            last=str(ex)
            time.sleep(2*(retry+1))
    return {'url':url,'error':last},None

def audit_table(kind,txt):
    if txt is None:return {'status':'DOWNLOAD_HOLD'}
    head=txt.splitlines()[:4]
    result={'status':'OK','line_count':len(txt.splitlines()),'first_line_preview':head[0][:300] if head else ''}
    if kind=='barcodes':
        bar=[l.strip().split('\t')[0] for l in txt.splitlines() if l.strip()]
        result.update({'n_values':len(bar),'unique':len(set(bar)),'example_prefixes':list(dict.fromkeys(x[:5] for x in bar[:50]))[:8]})
        result['_barcodes']=set(bar)
    if kind=='genes':
        genes=txt.splitlines()
        exact={g:[] for g in TARGETS}
        for i,ln in enumerate(genes):
            vals=re.split(r'[\t;, ]+',ln.strip())
            for gene in TARGETS:
                if gene in vals:exact[gene].append(i)
        result.update({'n_values':len(genes),'exact_target_index':exact,'head_fields':[l[:150] for l in genes[:3]]})
    if kind=='metadata':
        reader=csv.DictReader(io.StringIO(txt),delimiter='\t')
        col=reader.fieldnames or []
        candidates=[c for c in col if any(term in c.lower() for term in ['donor','patient','diagnos','disease','condition','sample','cell.type','celltype','barcode','status','group','orig.ident'])]
        counts={c:collections.Counter() for c in candidates}
        n=0
        col_examples={}
        for row in reader:
            n+=1
            for c in candidates:
                v=(row.get(c) or '').strip()
                if v:counts[c][v]+=1
        result.update({'n_rows':n,'n_cols':len(col),'headers':col[:80], 'profiling':{c:{'distinct':len(v),'top':v.most_common(12)} for c,v in counts.items()}})
    return result

def main():
    out=Path('output/ipf_gse136831_metadata');out.mkdir(parents=True,exist_ok=True)
    report={'accession':GEO,'audit_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'research_goal':'Independent donor-aware cell type context for FCN3, GPX3, CFH, HPR, CPN1, C3',
            'evidence_type':'GEO public metadata only - NO expression matrix parsed, NO replication',
            'files':{},'audits':{}}
    for kind,name in NAMES.items():
        meta,txt=grab(name)
        report['files'][kind]=meta
        res=audit_table(kind,txt)
        if '_barcodes' in res:res.pop('_barcodes')
        report['audits'][kind]=res
    valid=all(x.get('status')=='OK' for x in report['audits'].values())
    report['decision']='GO_TO_TARGETED_MTX_IF_COHORT_MAPPING_CONFIRMED' if valid else 'HOLD_METADATA_OR_DOWNLOAD'
    report['warnings']=[
        'Metadata row count is cells/libraries, NOT donor count. Multiple libraries from one donor must be collapsed.',
        'Celltype names must be harmonized to GSE135893 and donor-specific IPF versus Control subsets validated.',
        'Raw 2GB MTX not downloaded; no FCN3/GPX3 disease effect can be asserted from this audit.',
        'FCN3 lung expression cannot establish secreted plasma protein cell-of-origin, especially liver contributions.',
        'The two scRNA studies may differ by source, chemistry, processing and severity.'
    ]
    (out/'audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    msg=['# GSE136831 independent scRNA reference feasibility','',f"Decision: {report['decision']}",
         '',f"Files: {[(k, v.get('bytes', 'ERROR')) for k,v in report['files'].items()]}",
         '',f"Metadata rows: {report['audits'].get('metadata',{}).get('n_rows','NA')}",
         f"Gene records: {report['audits'].get('genes',{}).get('n_values','NA')}",
         f"Barcode records: {report['audits'].get('barcodes',{}).get('n_values','NA')}",
         '', 'See audit.json for full column diagnostics.','',
         'Metadata access does NOT reproduce FCN3 or GPX3 lung signal. Next: subset selected gene rows from MTX, donor-aware pseudobulk, celltype harmonization.']
    (out/'REPORT.md').write_text('\n'.join(msg)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['decision'],'counts':{k:v.get('n_rows',v.get('n_values')) for k,v in report['audits'].items()}},ensure_ascii=False))

if __name__=='__main__':main()