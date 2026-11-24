#!/usr/bin/env python3
"""Second-pass public dataset QC: explicit subject/ROI metadata and small-source quantification.
Not a clinical replication or drug response predictor; never infer counts/sample mapping by column order.
"""
import argparse,csv,gzip,io,json,re,collections,tarfile,hashlib
from pathlib import Path
import numpy as np
from geo_baseline import read_series,urlopen,series_path,GENES

def sample_profile(acc):
    samples,q=read_series(acc)
    return samples,q

def soft_names(acc):
    url=f'{series_path(acc)}/soft/{acc}_family.soft.gz'
    body,q=urlopen(url,maximum=18_000_000,timeout=70)
    out={'source':q,'gsm_file_pairs':[],'n_sample_headers':0}
    if body is None:return out
    try:
        txt=gzip.decompress(body).decode('utf8','replace')
        sample=None
        for line in txt.splitlines():
            if line.startswith('^SAMPLE = '): sample=line.split('=',1)[1].strip();out['n_sample_headers']+=1
            elif sample and line.startswith('!Sample_supplementary_file = '):
                val=line.split('=',1)[1].strip();out['gsm_file_pairs'].append({'accession':sample,'file':val})
    except Exception as e:out['parse_error']=str(e)[:240]
    return out

def therapy(out):
    acc='GSE213761'; samples,_=sample_profile(acc)
    src=soft_names(acc)
    target=f'{series_path(acc)}/suppl/GSE213761_RawCountsGEO.csv.gz'
    body,q=urlopen(target,maximum=18_000_000,timeout=70)
    audit={'accession':acc,'processed':q,'soft_mapping':{'fetch':src['source'],'n_samples':src['n_sample_headers'],'example_supplementary':src['gsm_file_pairs'][:8]},'status':'HOLD'}
    if body is None:return audit
    import pandas as pd
    mat=pd.read_csv(io.BytesIO(body),compression='gzip',index_col=0)
    audit.update({'gene_rows':mat.shape[0],'count_columns':mat.shape[1],
        'first_count_columns':list(mat.columns[:15]),
        'gene_ids_first':list(map(str,mat.index[:10])),
        'library_sizes':{'min':float(mat.sum(axis=0).min()),'median':float(mat.sum(axis=0).median()),'max':float(mat.sum(axis=0).max())},
        'missing_value_count':int(mat.isna().sum().sum()),
        'count_noninteger_rows_preview':int(((mat.iloc[:500].values-np.rint(mat.iloc[:500].values))!=0).sum())})
    # Direct GSM to matrix column mapping via SOFT file names; do NOT assume order.
    gsms={s['accession']:s for s in samples}
    mapping={}
    for link in src['gsm_file_pairs']:
        name=link['file'].split('/')[-1]
        for c in mat.columns:
            if c in name:
                mapping[c]=link['accession']
    # Individual RAW archive may preserve explicit GSM-to-internal-library links.
    raw_url=f'{series_path(acc)}/suppl/{acc}_RAW.tar'
    raw_body,raw_status=urlopen(raw_url,maximum=13_000_000,timeout=70)
    audit['raw_tar_probe']={'download':raw_status}
    if raw_body:
        try:
            with tarfile.open(fileobj=io.BytesIO(raw_body),mode='r:*') as tf:
                members=[m for m in tf.getmembers() if m.isfile()]
                member_names=[m.name for m in members]
                audit['raw_tar_probe']['member_count']=len(member_names)
                audit['raw_tar_probe']['first_names']=member_names[:15]
                audit['raw_tar_probe']['last_names']=member_names[-10:]
                canonical=lambda v: re.sub(r'[^a-z0-9]','',str(v).lower())
                canonical_cols={canonical(c):c for c in mat.columns}
                unmapped=[]
                for member in members:
                    name=member.name
                    gsm=re.search(r'GSM\d+',name)
                    if not gsm:continue
                    stem=re.sub(r'\.(?:txt|csv)(?:\.gz)?$', '', name.split('/')[-1],flags=re.I)
                    internal=re.sub(r'^GSM\d+_', '', stem)
                    col=canonical_cols.get(canonical(internal))
                    if col: mapping[col]=gsm.group(0)
                    else:unmapped.append(name)
                audit['raw_tar_probe']['unmatched_member_names_first']=unmapped[:15]
                audit['raw_tar_probe']['unmatched_matrix_columns_first']=[str(c) for c in mat.columns if c not in mapping][:15]
                # For every individual file, verify its gene ID/count rows against the matrix.
                # Missing identifiers are reconstructed only with *unique exact signatures*.
                variances=mat.var(axis=1).sort_values(ascending=False)
                topgenes=[str(v) for v in variances.index[:60]]
                signature_cols={c:tuple(float(mat.loc[int(g) if mat.index.dtype.kind in 'iu' else g,c]) for g in topgenes[:25]) for c in mat.columns}
                remaining={c for c in mat.columns if c not in mapping}
                for member in members:
                    gsm_match=re.search(r'GSM\d+',member.name)
                    if not gsm_match:continue
                    gsm=gsm_match.group(0)
                    if gsm in mapping.values():continue
                    binary=tf.extractfile(member).read()
                    if member.name.endswith('.gz'):binary=gzip.decompress(binary)
                    lines=binary.decode('utf8','replace').splitlines()
                    if 'sample_file_first_lines' not in audit['raw_tar_probe']:
                        audit['raw_tar_probe']['sample_file_first_lines']=lines[:6]
                    observed={}
                    for line in lines:
                        cells=re.split(r'[\t,; ]+',line.strip())
                        if len(cells)>=2 and cells[0] in topgenes:
                            try:observed[cells[0]]=float(cells[1])
                            except ValueError:continue
                    if len(observed)<10:continue
                    matched_cols=[]
                    genes=sorted(observed)[:30]
                    for col in remaining:
                        try:
                            if all(float(mat.loc[int(g) if mat.index.dtype.kind in 'iu' else g,col])==observed[g] for g in genes):matched_cols.append(col)
                        except (KeyError,ValueError):continue
                    if len(matched_cols)==1:
                        mapping[matched_cols[0]]=gsm;remaining.remove(matched_cols[0])
                audit['raw_tar_probe']['numeric_verified_map_n']=len(mapping)
        except Exception as e:audit['raw_tar_probe']['error']=str(e)[:350]
    audit['explicit_column_map_n']=len(mapping)
    audit['n_unique_gsm_from_map']=len(set(mapping.values()))
    audit['exemplar_map']=dict(list(mapping.items())[:8])
    if len(mapping)!=mat.shape[1] or len(set(mapping.values()))!=len(mapping):
        # Check series source names / supplementary links; never fabricate mapping.
        audit['status']='HOLD_MATRIX_SAMPLE_MAPPING'
        return audit
    matched=[]
    for col,gsm in mapping.items():
        s=gsms.get(gsm)
        if not s:continue
        m=re.match(r'^(\d+)_(Lesional|Nonlesional)_(V1|V12)$',s['title'],re.I)
        if m:matched.append({'column':col,'gsm':gsm,'patient':m.group(1),'tissue':m.group(2).lower(),'visit':m.group(3),'title':s['title'],
        'responder':next((v.split(':',1)[-1].strip() for k,v in s.items() if 'characteristics' in k and v.lower().startswith('responder:')),'UNKNOWN')})
    if len(matched)!=mat.shape[1]:audit['status']='HOLD_INCOMPLETE_PATIENT_MAPPING';return audit
    keys=collections.defaultdict(dict)
    for s in matched:keys[(s['patient'],s['tissue'])][s['visit']]=s
    pairs=[(v['V1'],v['V12']) for v in keys.values() if set(v)=={'V1','V12'}]
    les_pairs=[(a,b) for a,b in pairs if a['tissue']=='lesional']
    audit['patient_tissue_pairs']=len(pairs);audit['lesional_pairs']=len(les_pairs)
    audit['responder_pairs']=dict(collections.Counter(a['responder'] for a,b in les_pairs))
    # Save only per-library QC and exploratory globally normalized selected gene effects.
    import pandas as pd
    libs=mat.sum(axis=0).astype(float)
    logcpm=np.log2(mat.div(libs,axis=1)*1e6+1)
    # Gene index is NCBI Entrez GeneID; do not misinterpret as symbols.
    entrez={'IL17A':'3605','IL1B':'3553','TNF':'7124','CXCL8':'3576',
            'CXCL13':'10563','S100A8':'6279','S100A9':'6280','KRT17':'3872','CCL20':'6364'}
    candidate=[]
    for g in GENES['hs_therapy']:
        idx=[id for id in mat.index if str(id).strip()==entrez.get(g,'__MISSING__')]
        if not idx:continue
        row=logcpm.loc[idx[0]]
        diffs=np.array([float(row[b['column']]-row[a['column']]) for a,b in les_pairs])
        responder_deltas=[];nonresponder_deltas=[]
        for a,b in les_pairs:
            delta=float(row[b['column']]-row[a['column']])
            (responder_deltas if a['responder']=='YES' else nonresponder_deltas).append(delta)
        candidate.append({'gene':g,'entrez_id':entrez[g], 'n_matched_pairs':len(diffs),
                          'median_delta_log2cpm':round(float(np.median(diffs)),4),
                          'positive_pairs':int((diffs>0).sum()),'negative_pairs':int((diffs<0).sum()),
                          'n_responder_pairs':len(responder_deltas),'responder_median_delta':round(float(np.median(responder_deltas)),4) if responder_deltas else None,
                          'n_nonresponder_pairs':len(nonresponder_deltas),'nonresponder_median_delta':round(float(np.median(nonresponder_deltas)),4) if nonresponder_deltas else None})
    audit['candidate_gene_analysis']=candidate
    # Explicitly mapped 14 patient-paired lesional specimens for DESeq2.
    les_pairs.sort(key=lambda x:int(x[0]['patient']))
    selected=[s for a,b in les_pairs for s in (a,b)]
    columns=[s['column'] for s in selected]
    meta_rows=[{'column':s['column'],'patient':s['patient'],'visit':s['visit'], 'responder':s['responder'],'gsm':s['gsm']} for s in selected]
    if len(set(columns))==len(columns) and len(columns)==len(les_pairs)*2:
        mat.loc[:,columns].to_csv(out/'hs_therapy_paired_counts.csv',index=True,index_label='gene_entrez_id')
        with open(out/'hs_therapy_paired_metadata.csv','w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=list(meta_rows[0]));w.writeheader();w.writerows(meta_rows)
        audit['paired_model_input']='READY_COUNT_VERIFIED'
    if candidate:
        with open(out/'hs_therapy_paired_exploratory.csv','w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=list(candidate[0]));w.writeheader();w.writerows(candidate)
    audit['status']='QC_PAIRED_DESIGN_READY_EXPLORATORY' if len(les_pairs)>=5 else 'HOLD_TOO_FEW_PAIRS'
    with open(out/'hs_therapy_explicit_sample_map.csv','w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(matched[0]));w.writeheader();w.writerows(matched)
    return audit

def nmosd(out):
    acc='GSE242232';samples,_=sample_profile(acc)
    audit={'accession':acc,'cells':len(samples)}
    bytitle=collections.Counter(re.match(r'^(T-[A-Za-z0-9]+|Treg)',s['title']).group(1) if re.match(r'^(T-[A-Za-z0-9]+|Treg)',s['title']) else 'OTHER' for s in samples)
    audit['title_prefix_counts']=dict(bytitle)
    tracks={}
    for s in samples:
        m=re.search(r'\[(L\d+_Track-\d+)\]',s['title'])
        if m:tracks[m.group(1)]=s
    audit['track_ids_from_titles']=len(tracks)
    body,status=urlopen(f'{series_path(acc)}/suppl/{acc}_RawExpressionCountTable.csv.gz',maximum=18_000_000,timeout=70)
    audit['processed']=status
    if body is None:return audit
    # Stream gzipped CSV, no cells treated as patient replicates. Extract *only* target gene vectors.
    gz=gzip.GzipFile(fileobj=io.BytesIO(body)); reader=csv.reader(io.TextIOWrapper(gz,encoding='utf-8-sig',errors='replace'))
    header=next(reader)
    audit['columns']=len(header);audit['header_head']=header[:10]
    mapped={col:tracks[col] for col in header if col in tracks}
    audit['matched_track_columns']=len(mapped)
    audit['unique_matched_GSM']=len({s['accession'] for s in mapped.values()})
    gene_ids={'FOXP3':'ENSG00000049768','PDCD1':'ENSG00000188389',
              'CTLA4':'ENSG00000163599','LAG3':'ENSG00000089692',
              'CD40LG':'ENSG00000102245','IL21':'ENSG00000138684'}
    id_to_symbol={v:k for k,v in gene_ids.items()}
    candidate=[];first_genes=[]
    for i,row in enumerate(reader):
        if i<8: first_genes.append(row[:4])
        if len(row)!=len(header):continue
        # Gene_id is symbol or Ensembl depending on GEO export.
        names=[str(row[j]).upper() for j in range(min(4,len(row)))]
        hits=[id_to_symbol[row[1].upper()]] if len(row)>1 and row[1].upper() in id_to_symbol else []
        if not hits:continue
        vals=[]
        for j,c in enumerate(header):
            if c in mapped:
                try: vals.append((c,float(row[j])))
                except ValueError:pass
        # by *stimulation label*, not by disease status
        stats=collections.defaultdict(list)
        for c,v in vals:
            s=mapped[c];m=re.match(r'^(T-[A-Za-z0-9]+|Treg)',s['title'])
            stats[m.group(1) if m else 'Other'].append(v)
        for g in hits:
            for k,v in stats.items():
                candidate.append({'gene':g,'entrez_stable_id':gene_ids[g],'stimulation_title':k,'cells_with_read':len(v),
                  'fraction_detected':round(sum(vv>0 for vv in v)/len(v),4), 'mean_raw_count':round(sum(v)/len(v),3)})
    audit['first_gene_row_names']=first_genes
    audit['candidate_gene_by_title']=candidate
    audit['status']='CELL_LEVEL_MECHANISM_ONLY_NO_PATIENT_DIFFERENTIAL' if len(mapped)>=800 else 'HOLD_CELL_MAPPING'
    with open(out/'nmosd_gene_descriptive.csv','w',newline='',encoding='utf8') as f:
        if candidate:
            w=csv.DictWriter(f,fieldnames=list(candidate[0]));w.writeheader();w.writerows(candidate)
    return audit

def mg(out):
    acc='GSE227835';samples,_=sample_profile(acc)
    groups=collections.defaultdict(list)
    for s in samples:
        title=s.get('title','')
        if title.startswith('AChR-positive'):key='AChR-positive'
        elif title.startswith('Healthy control'):key='Healthy control'
        elif title.startswith('Seronegative MG pre'):key='Seronegative pre'
        elif title.startswith('Seronegative MG post'):key='Seronegative post'
        else:key='Other'
        groups[key].append(title)
    paired=[]
    before=[re.search(r'_(\d+)$',s) for s in groups['Seronegative pre']]
    after=[re.search(r'_(\d+)$',s) for s in groups['Seronegative post']]
    paired=sorted({m.group(1) for m in before if m}&{m.group(1) for m in after if m})
    return {'samples':len(samples),'groups':{k:len(v) for k,v in groups.items()},'seronegative_prepost_title_ids':len(paired),
            'GSE297547':'Series matrix fetch HTTP 403 in first audit, GEO HTML lists 25 samples; no MuSK status proven.',
            'status':'HOLD_NO_MUSK_SPECIFIC_GENOTYPED_OR_SEROTYPED_PATIENTS'}

def hs_spatial(out):
    acc='GSE294009';samples,_=sample_profile(acc)
    count=collections.Counter();subject=collections.defaultdict(set)
    for s in samples:
        disease=next((v.split(':',1)[-1].strip() for k,v in s.items() if 'characteristics' in k and v.lower().startswith('disease:')),'UNSPECIFIED')
        count[disease]+=1
        t=s.get('title','')
        prefix=t.split('|')[0].strip()
        subject[prefix].add(disease)
    diag=[]
    for d,n in count.most_common():
        stage=('Early' if d.lower().startswith('early') else 'Late' if d.lower().startswith('late') else 'Control' if d.lower().startswith('control') else 'Other')
        diag.append({'condition':d,'stage':stage,'n_rois':n})
    with open(out/'hs_spatial_roi_strata.csv','w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(diag[0]));w.writeheader();w.writerows(diag)
    # Use only known categories, do not equate internal slide prefix to an independent patient id.
    stages=collections.Counter()
    for d,n in count.items():
        stages['Early' if d.startswith('Early') else 'Late' if d.startswith('Late') else 'Control' if d.startswith('Control') else 'Other']+=n
    return {'total_geo_samples':len(samples),'stage_roi_counts':dict(stages),'roi_strata':diag,
            'unique_slide_title_prefixes':len(subject),'title_prefixes':dict(collections.Counter(s.get('title','').split('|')[0].strip() for s in samples).most_common(25)),
            'external_GSE249729':'2 HS lesions and one normal control; different platform and donors',
            'status':'ROI_STRATA_CONFIRMED_HOLD_PATIENT_AND_NUMERIC_MAPPING'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--role',required=True,choices=['hs_therapy','nmosd','musk_mg','hs_spatial']);p.add_argument('--out',required=True);a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    fn={'hs_therapy':therapy,'nmosd':nmosd,'musk_mg':mg,'hs_spatial':hs_spatial}[a.role]
    try:res=fn(out)
    except Exception as e:res={'status':'ERROR','error':str(e)[:300]}
    (out/'REFINED_QC.json').write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding='utf8')
    report=[f'# {a.role} second-pass evidence', '',f"**Status: {res.get('status')}**",'', 'Public GEO data only. No outreach or clinical result prediction.', '', '## Result', '', '```json',json.dumps(res,ensure_ascii=False,indent=2)[:12000], '```']
    (out/'REFINED_REPORT.md').write_text('\n'.join(report),encoding='utf8')
    print(json.dumps({'role':a.role,'status':res.get('status'),'key_counts':{k:res.get(k) for k in ('samples','groups','lesional_pairs','matched_track_columns','stage_roi_counts')}},ensure_ascii=False))
    # Scientific HOLD is an expected data-access outcome; technical ERROR is not.
    if res.get('status') == 'ERROR':
        raise SystemExit(1)

if __name__=='__main__':main()