#!/usr/bin/env python3
"""Recheck pinned public-study artifacts offline; no FVC or clinical inference."""
import argparse
import csv
import hashlib
import io
import json
import math
import statistics as st
import zipfile
from pathlib import Path

ARTIFACTS = {
    'plasma': (14909794711, '47061878edb85cd45e59bbf2be7d82a876ef5af5cff0619f222eea9a1633372f'),
    'source': (14910154655, 'e3c87b83537833a890a8b41f32d0850c3df872365ceafec3c0d3729a01501cbf'),
    'lung_a': (14911122272, '24daeb86b89d68f6077def5482145d324fdaa1581dff654ef7c46d52649e390a'),
    'lung_b': (14912096747, 'adaf366001dfe8b508d83c079291c97f1073d2a92bc325a6c025693eabcb7335'),
}
RAW_SHA = '1df37c4f3a2e130bd02ea87dc5da47ba1a6701efe24c657b358cc721bbd3f334'
TARGETS = {'FCN3': 'O75636', 'GPX3': 'P22352'}
MIN_CELLS = 20
MIN_DONORS = 5


def member(directory, kind, name):
    aid, sha = ARTIFACTS[kind]
    path = Path(directory)/f'github-actions-artifact-{aid}.zip'
    b = path.read_bytes()
    if hashlib.sha256(b).hexdigest() != sha:
        raise ValueError('STOP: artifact checksum mismatch: '+kind)
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        hits = [x for x in z.infolist() if x.filename == name]
        if len(hits) != 1 or hits[0].file_size > 5_000_000:
            raise ValueError('STOP: missing, ambiguous or oversized artifact member: '+name)
        return z.read(hits[0])


def table(b):
    return list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))))


def loo_effect(values, n_control):
    if not 2 <= n_control <= len(values)-2 or not all(math.isfinite(x) for x in values):
        raise ValueError('invalid groups for leave-one-out')
    full = st.mean(values[n_control:])-st.mean(values[:n_control])
    folds = [st.mean([v for j,v in enumerate(values) if j>=n_control and j!=i])-
             st.mean([v for j,v in enumerate(values) if j<n_control and j!=i])
             for i in range(len(values))]
    return {'log2fc': full, 'loo_min': min(folds), 'loo_max': max(folds),
            'same_nonzero_direction': sum(v*full>0 for v in folds), 'n_folds': len(folds)}


def plasma_review(directory):
    raw = member(directory,'plasma','ipf_pxd010965_quant/Proteins_public_source.csv')
    if hashlib.sha256(raw).hexdigest() != RAW_SHA:
        raise ValueError('STOP: unexpected raw abundance source SHA')
    rows = list(csv.reader(io.StringIO(raw.decode('utf-8-sig'))))
    if (len(rows)!=170 or any(len(r)!=49 for r in rows) or rows[0][12]!='Raw abundance' or
        [(i,x) for i,x in enumerate(rows[1]) if x in ('CONTROL','CASE')]!=[(12,'CONTROL'),(32,'CASE')] or
        len(set(rows[2][12:]))!=37 or not all(rows[2][12:]) or len({r[0] for r in rows[3:]})!=167):
        raise ValueError('STOP: unexpected plasma assay or protein identity')
    vals = [[float(x) for x in r[12:]] for r in rows[3:]]
    if not all(math.isfinite(x) and x>0 for r in vals for x in r):
        raise ValueError('STOP: invalid raw abundance')
    logs = [[math.log2(x) for x in r] for r in vals]
    med = [st.median(c) for c in zip(*logs)]
    archived = table(member(directory,'plasma','ipf_pxd010965_reanalysis/protein_reanalysis.csv'))
    by_acc = {r['accession']:r for r in archived}
    if len(archived)!=167 or len(by_acc)!=167 or set(by_acc)!={r[0] for r in rows[3:]}:
        raise ValueError('STOP: archived protein identities mismatch')
    results = {}
    for gene,acc in TARGETS.items():
        ix = next(i for i,r in enumerate(rows[3:]) if r[0]==acc)
        ref = by_acc[acc]
        if ref['gene']!=gene:
            raise ValueError('STOP: accession/gene identity mismatch')
        raw_loo = loo_effect(logs[ix],20)
        centered_loo = loo_effect([x-y for x,y in zip(logs[ix],med)],20)
        for v,k,q in ((raw_loo,'log2fc_case_minus_control','welch_fdr'),
                      (centered_loo,'sample_median_centered_log2fc','sample_median_centered_welch_fdr')):
            if not math.isclose(v['log2fc'],float(ref[k]),abs_tol=1e-9):
                raise ValueError('STOP: archived and replayed effects differ')
            v['archived_welch_fdr'] = float(ref[q])
            if not 0<=v['archived_welch_fdr']<=1:
                raise ValueError('STOP: invalid archived FDR')
        results[gene] = {'raw_log2_abundance':raw_loo,'sample_median_centered':centered_loo}
    return results


def source_map(rows):
    out={}
    for r in rows:
        key=r['donor']
        if not key or key in out or not r['diagnosis'] or not r['source']:
            raise ValueError('STOP: duplicate or incomplete donor source map')
        out[key]=(r['diagnosis'],r['source'])
    return out


def donor_summary(rows, gene, celltype, *, category=None, sources=None, site=None):
    groups={'IPF':[],'Control':[]}
    seen=set()
    for r in rows:
        if r['gene']!=gene or r['celltype']!=celltype or (category is not None and r['type']!=category):
            continue
        group=r.get('diagnosis',r.get('group')); donor=r['donor']
        if group not in groups or not donor or donor in seen:
            raise ValueError('STOP: ambiguous donor or group in selected celltype')
        seen.add(donor)
        if sources is not None:
            if donor not in sources or sources[donor][0]!=group:
                raise ValueError('STOP: donor/source group mismatch')
            if sources[donor][1]!=site:
                continue
        n=int(r['n_cells']); cpm=float(r['pseudobulk_cpm'])
        if n<1 or not math.isfinite(cpm) or cpm<0:
            raise ValueError('STOP: invalid donor-level CPM or cell count')
        if n>=MIN_CELLS: groups[group].append(cpm)
    a,b=groups['IPF'],groups['Control']
    ready=min(len(a),len(b))>=MIN_DONORS
    return {'gene':gene,'celltype':celltype,'category':category,'site':site,
            'ipf_donors':len(a),'control_donors':len(b),
            'log2_median_ratio_plus1':math.log2((1+st.median(a))/(1+st.median(b))) if ready else None,
            'status':'DESCRIPTIVE_ONLY' if ready else 'HOLD_SPARSE_DONORS'}


def evaluate(directory):
    a=json.loads(member(directory,'lung_a','summary.json'))
    b=json.loads(member(directory,'lung_b','audit.json'))
    src=json.loads(member(directory,'source','source_audit.json'))
    if (a['decision']!='GO_DESCRIPTIVE_DONOR_ANALYSIS' or a['metadata_sha256']!=src['source_sha256'] or
        a['mtx_parsed_entries']!=a['mtx_dims']['nonzero_entries'] or
        b['decision']!='GO_DESCRIPTIVE_EXTERNAL_COMPARISON' or b['parsed_nnz']!=b['mtx_dimensions']['nonzero']):
        raise ValueError('STOP: inconsistent archived source or matrix completeness gate')
    sources=source_map(table(member(directory,'source','donor_source_mapping.csv')))
    ar=table(member(directory,'lung_a','donor_celltype_gene.csv'))
    br=table(member(directory,'lung_b','donor_gene_celltype.csv'))
    if len(ar)!=a['n_donor_celltype_gene_rows'] or len(br)!=b['n_donor_gene_rows']:
        raise ValueError('STOP: archived donor table row count mismatch')
    matched=[donor_summary(ar,g,c,sources=sources,site='Vanderbilt') for g,c in
             [('FCN3','Endothelial Cells'),('GPX3','Endothelial Cells'),('GPX3','Fibroblasts'),('GPX3','Myofibroblasts')]]
    external=[donor_summary(br,g,c,category=k) for g,c,k in
              [('FCN3','Endothelial','category'),('GPX3','Endothelial','category'),('GPX3','Stromal','category'),
               ('GPX3','Fibroblast','subclass'),('FCN3','Epithelial','category'),('GPX3','Epithelial','category')]]
    return {'decision':'ARCHIVED_EVIDENCE_RECHECK__HOLD_CLINICAL_PROGNOSIS',
            'artifacts':ARTIFACTS,'min_cells_per_donor_celltype':MIN_CELLS,'min_donors_per_group':MIN_DONORS,
            'plasma':plasma_review(directory),'GSE135893_same_source':matched,'GSE136831':external,
            'clinical_prognosis_validated':False,'protein_tissue_origin_established':False,
            'limitations':['Same-cohort plasma sensitivity; FDR values copied from prior output, not recomputed.',
                           'Lung RNA cohorts are not paired with plasma, cell labels not fully harmonized.',
                           'Donor-level descriptive medians only; no FVC or clinical covariate validation.']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--artifact-dir',required=True);p.add_argument('--out',required=True,type=Path)
    a=p.parse_args(); result=evaluate(a.artifact_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
