#!/usr/bin/env python3
"""GSE136831 donor-aware six-gene RNA count audit. No dense 312k x 45k matrix.

Public raw counts and annotations; only IPF / control, donors independent.
This is mechanism triangulation, not paired plasma–lung evidence or biomarker proof.
"""
from __future__ import annotations
import argparse, csv, datetime as dt, gzip, hashlib, io, json, math, statistics, time, urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ACC='GSE136831'
BASE=f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/{ACC}/suppl/'
NCBI=f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={ACC}&format=file&file='
FILES={
  'meta':'GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz',
  'genes':'GSE136831_AllCells.GeneIDs.txt.gz',
  'barcodes':'GSE136831_AllCells.cellBarcodes.txt.gz',
  'matrix':'GSE136831_RawCounts_Sparse.mtx.gz',
}
TARGETS=('FCN3','GPX3','CFH','HPR','CPN1','C3')

def download_small(name):
    failures=[]
    for u in (NCBI+name, BASE+name):
        try:
            req=urllib.request.Request(u,headers={'User-Agent':'PublicIPFResearch-GSE136831/1.0'})
            with urllib.request.urlopen(req,timeout=90) as f:
                data=f.read(25_000_001)
            if len(data)>25_000_000 or data[:2]!=b'\x1f\x8b':raise ValueError('Unexpected file size or non-GZIP')
            return gzip.decompress(data).decode('utf-8-sig'),{'url':u,'compressed_bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
        except Exception as e:
            failures.append(str(e)[:150]);time.sleep(1)
    raise RuntimeError('Small file unavailable: '+name+'; '.join(failures))

def mtx_response():
    errors=[]
    for u in (BASE+FILES['matrix'],NCBI+FILES['matrix']):
        try:
            req=urllib.request.Request(u,headers={'User-Agent':'PublicIPFResearch-GSE136831/1.0','Accept':'application/octet-stream,*/*'})
            f=urllib.request.urlopen(req,timeout=180)
            if 'html' in f.headers.get('content-type',''):
                f.close();raise ValueError('Unexpected HTML response')
            return f,u
        except Exception as e:
            errors.append(type(e).__name__+': '+str(e)[:200]);time.sleep(2)
    raise RuntimeError('MTX download failed: '+' | '.join(errors))

class HashReader:
    def __init__(self,resp):self.r=resp;self.digest=hashlib.sha256();self.n=0
    def read(self,n=-1):
        b=self.r.read(n)
        if b:self.digest.update(b);self.n+=len(b)
        return b

def get_median(z):return statistics.median(z) if z else None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='output/ipf_gse136831_targeted');a=ap.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    audit={'accession':ACC,'date_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'targets':TARGETS,'state':'HOLD_NOT_RUN', 'scope':'Independent lung scRNA IPF/control, not paired plasma'}
    mtxt,mprov=download_small(FILES['meta']); gtxt,gprov=download_small(FILES['genes']);btxt,bprov=download_small(FILES['barcodes'])
    audit['input_provenance']={'meta':mprov,'genes':gprov,'barcodes':bprov}
    reader=csv.DictReader(io.StringIO(mtxt),delimiter='\t')
    keys=['CellBarcode_Identity','Subject_Identity','Disease_Identity','CellType_Category','Subclass_Cell_Identity','nUMI']
    if not set(keys).issubset(reader.fieldnames or []):raise RuntimeError('Missing metadata columns '+str(set(keys)-set(reader.fieldnames or [])))
    metadata=list(reader)
    bcodes=[x.strip() for x in btxt.splitlines() if x.strip()]
    if len(bcodes)!=312928 or len(set(bcodes))!=len(bcodes):raise RuntimeError('barcode count/uniqueness error')
    by_bar={r['CellBarcode_Identity'].strip():r for r in metadata}
    if len(by_bar)!=len(metadata) or len(by_bar)!=len(bcodes) or set(by_bar)!=set(bcodes):raise RuntimeError('Barcode join not exact')
    gene_reader=list(csv.reader(io.StringIO(gtxt),delimiter='\t'))
    if len(gene_reader)!=45948:raise RuntimeError('Unexpected gene table length')
    if gene_reader[0][:2]!=['Ensembl_GeneID','HGNC_EnsemblAlt_GeneID']:
        raise RuntimeError('Gene annotation header changed: '+str(gene_reader[0]))
    # Crucial: 1 header row in gene file; MTX row #1 maps to gene_reader[1].
    target_row_idx={i:row[1].strip() for i,row in enumerate(gene_reader[1:],start=1) if len(row)>=2 and row[1].strip() in TARGETS}
    if set(target_row_idx.values())!=set(TARGETS):raise RuntimeError('Incomplete or duplicated target gene mapping: '+str(target_row_idx))
    labels=Counter(r['Disease_Identity'].strip() for r in metadata)
    donor_groups=defaultdict(set)
    for r in metadata:donor_groups[r['Subject_Identity'].strip()].add(r['Disease_Identity'].strip())
    if any(len(gs)>1 for gs in donor_groups.values()):raise RuntimeError('Donor assigned multiple diseases')
    audit.update({'metadata_cells':len(metadata),'barcodes':len(bcodes),'genes':len(gene_reader)-1,'gene_row_indices':{gene:i for i,gene in target_row_idx.items()},
       'disease_cell_counts':dict(labels),'disease_donor_counts':dict(Counter(next(iter(v)) for v in donor_groups.values()))})
    # Group 1: categories; group 2: subtypes (names not guaranteed equivalent to GSE135893).
    idx_to_cell={};groups=defaultdict(lambda:{'cells':0,'umi':0.0})
    for col,bc in enumerate(bcodes,start=1):
        r=by_bar[bc]; diag=r['Disease_Identity'].strip()
        if diag not in ('IPF','Control'):continue
        donor=r['Subject_Identity'].strip(); major=r['CellType_Category'].strip();sub=r['Subclass_Cell_Identity'].strip()
        if major=='Multiplet' or not major:continue
        try:umi=float(r['nUMI'])
        except (ValueError,TypeError):continue
        if not math.isfinite(umi) or umi<=0:continue
        for typ,ct in [('category',major),('subclass',sub)]:
            if typ=='subclass' and not ct:continue
            key=(typ,ct,diag,donor)
            groups[key]['cells']+=1; groups[key]['umi']+=umi
        idx_to_cell[col]=(diag,donor,major,sub)
    exact_rows=set(target_row_idx)
    byte_rows={str(i).encode():target_row_idx[i] for i in exact_rows}
    counts=defaultdict(lambda:[0.0,0])
    response,url=mtx_response();hr=HashReader(response)
    nnz_seen=0;hits=0
    try:
        with gzip.GzipFile(fileobj=hr,mode='rb') as gz:
            if not gz.readline().startswith(b'%%MatrixMarket matrix coordinate'):
                raise RuntimeError('MTX not coordinate format')
            while True:
                ln=gz.readline()
                if not ln:raise RuntimeError('Truncated MTX before dimensions')
                if ln.startswith(b'%'):continue
                dims=list(map(int,ln.split()));break
            if len(dims)!=3 or dims[0]!=len(gene_reader)-1 or dims[1]!=len(bcodes):raise RuntimeError(f'MTX shape={dims} vs genes/barcodes {(len(gene_reader)-1,len(bcodes))}')
            for ln in gz:
                if not ln or ln.startswith(b'%'):continue
                nnz_seen+=1
                # MTX is column-oriented sparse input: inspect first field without parsing every row.
                sep=ln.find(b' ')
                if sep<0:sep=ln.find(b'\t')
                if sep<0:raise RuntimeError('Malformed MTX entry, no field separator')
                gene=byte_rows.get(ln[:sep])
                if gene is None:continue
                fields=ln.split()
                if len(fields)!=3:raise RuntimeError('Malformed target MTX line')
                col=int(fields[1]);value=float(fields[2])
                if value<0 or not math.isfinite(value):raise RuntimeError('Invalid count')
                if col not in idx_to_cell:continue
                diag,donor,major,sub=idx_to_cell[col]
                if value>0:
                    counts[('category',major,diag,donor,gene)][0]+=value;counts[('category',major,diag,donor,gene)][1]+=1
                    if sub:
                        counts[('subclass',sub,diag,donor,gene)][0]+=value;counts[('subclass',sub,diag,donor,gene)][1]+=1
                    hits+=1
            if nnz_seen!=dims[2]:raise RuntimeError(f'MTX truncated: {nnz_seen} != {dims[2]}')
            audit['mtx_dimensions']={'gene_rows':dims[0],'barcodes':dims[1],'nonzero':dims[2]}
            audit['parsed_nnz']=nnz_seen;audit['target_nonzero_entries']=hits
            audit['mtx_sha256']=hr.digest.hexdigest();audit['matrix_compressed_bytes']=hr.n;audit['mtx_source']=url
    finally:response.close()
    per_donor=[]
    for (typ,ct,diag,donor),s in sorted(groups.items()):
        if s['cells']<15:continue
        for gene in TARGETS:
            g=counts.get((typ,ct,diag,donor,gene),[0.0,0])
            per_donor.append({'type':typ,'celltype':ct,'group':diag,'donor':donor,'gene':gene,'n_cells':s['cells'],
                              'nonzero_cells':g[1],'detection_pct':round(100*g[1]/s['cells'],4),'pseudobulk_cpm':round(1e6*g[0]/s['umi'],6)})
    def write_csv(path,records,cols):
        with (out/path).open('w',newline='',encoding='utf-8') as h:
            w=csv.DictWriter(h,fieldnames=cols);w.writeheader();w.writerows(records)
    write_csv('donor_gene_celltype.csv',per_donor,['type','celltype','group','donor','gene','n_cells','nonzero_cells','detection_pct','pseudobulk_cpm'])
    grouped=defaultdict(lambda:defaultdict(list))
    for r in per_donor:grouped[(r['type'],r['celltype'],r['gene'])][r['group']].append(float(r['pseudobulk_cpm']))
    summary=[]
    for (typ,ct,gene),g in sorted(grouped.items()):
        a,b=g.get('IPF',[]),g.get('Control',[])
        if not a or not b:continue
        logfc=math.log2((1+get_median(a))/(1+get_median(b)))
        # Direction only; descriptive medians, no cell-level pseudo-replicated tests.
        x={'type':typ,'celltype':ct,'gene':gene,'n_ipf_donors':len(a),'n_control_donors':len(b),
          'ipf_median_cpm':round(get_median(a),4),'control_median_cpm':round(get_median(b),4),
          'log2_median_ratio_plus1':round(logfc,5),'status':'DESCRIPTIVE' if min(len(a),len(b))>=5 else 'SPARSE_DONOR_HOLD'}
        summary.append(x)
    write_csv('celltype_candidate_summary.csv',summary,['type','celltype','gene','n_ipf_donors','n_control_donors','ipf_median_cpm','control_median_cpm','log2_median_ratio_plus1','status'])
    audit['n_donor_gene_rows']=len(per_donor);audit['n_summary_rows']=len(summary)
    audit['decision']='GO_DESCRIPTIVE_EXTERNAL_COMPARISON'
    audit['warnings']=['Independent dataset, not matched PXD010965 plasma samples',
         'Celltype_Category and Subclass_Cell_Identity must be harmonized to Habermann GSE135893 before replication claims',
         'No treatment status, age/sex/collection batch covariate correction, or true patient prognosis',
         'Never treat individual cells as independent statistical replicates',
         'No plasma protein origin inference from lung RNA alone']
    (out/'audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    (out/'REPORT.md').write_text('# Independent GSE136831 IPF target transcript audit\n\n'+
       f"- Status: {audit['decision']}\n- Matrix {audit['mtx_dimensions']}\n- IPF/Control donors {audit['disease_donor_counts']}\n"+
       '- outputs: donor_gene_celltype.csv; celltype_candidate_summary.csv\n'+
       '- Counts are donor-level pseudobulk CPM medians with 15-cell minimum. No formal differential expression P values.\n',encoding='utf-8')
    print(json.dumps({'decision':audit['decision'],'parsed_nnz':nnz_seen,'target_nonzero_entries':hits,'n_summary_rows':len(summary),'selected':[r for r in summary if r['type']=='category' and r['gene'] in ('FCN3','GPX3') and r['celltype'] in ('Endothelial','Stromal','Epithelial')]},ensure_ascii=False))

if __name__=='__main__':main()