#!/usr/bin/env python3
"""Exact 7-vs-7 group-label permutation of paired within-person log-normalized changes.
Strictly preselected genes, and this is an exploratory non-predictive comparison.
"""
from pathlib import Path
from itertools import combinations
from math import comb
import sys
import argparse
import numpy as np
import pandas as pd

def check_framing(df):
    assert df.shape[0]==28 and df.patient.nunique()==14
    assert set(df.visit)=={'V1','V12'} and set(df.responder)=={'NO','YES'}
    p=df.groupby('patient').agg(visits=('visit','nunique'),response=('responder','nunique'))
    assert (p.visits==2).all() and (p.response==1).all()
    assert df[df.visit=='V1'].responder.value_counts().to_dict()=={'YES':7,'NO':7}
    assert df.column.nunique()==28

def exact_perm_difference(deltas, is_resp):
    n=len(deltas);k=sum(is_resp)
    assert n==14 and k==7 and comb(n,k)==3432
    observed=float(np.mean(deltas[is_resp])-np.mean(deltas[~is_resp]))
    matches=0
    positions=np.arange(n)
    for a in combinations(range(n),k):
        mask=np.zeros(n,dtype=bool);mask[list(a)]=True
        diff=float(np.mean(deltas[mask])-np.mean(deltas[~mask]))
        if abs(diff)>=abs(observed)-1e-10:matches+=1
    return observed, matches/comb(n,k)

def bh(p):
    valid=[i for i,x in enumerate(p) if np.isfinite(x)]
    ranks=sorted(valid,key=lambda i:p[i]);ret=[float('nan')]*len(p)
    bound=1.
    for j in range(len(ranks)-1,-1,-1):
        i=ranks[j];r=j+1
        bound=min(bound,p[i]*len(ranks)/r);ret[i]=bound
    return ret

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    args=parser.parse_args();out=Path(args.out)
    m=pd.read_csv(out/'hs_therapy_paired_metadata.csv',dtype={'patient':str,'column':str})
    check_framing(m)
    expr=pd.read_csv(out/'hs_targets_log2_normalized_counts.csv',dtype={'entrez_gene_id':str})
    assert set(m.column).issubset(set(expr.columns))
    rows=[]
    for _,gene in expr.iterrows():
        vec=m[['column','patient','visit','responder']].copy()
        vec['expression']=vec.column.map({c:float(gene[c]) for c in m.column})
        piv=vec.pivot(index='patient',columns='visit',values='expression')
        resp=vec[vec.visit=='V1'].set_index('patient').responder.reindex(piv.index)
        assert not piv.isnull().any().any() and resp.notna().all()
        v=(piv.V12-piv.V1).to_numpy()
        is_yes=(resp=='YES').to_numpy()
        interaction,pchange=exact_perm_difference(v,is_yes)
        baseline,pbaseline=exact_perm_difference(piv.V1.to_numpy(),is_yes)
        rows.append({'gene':gene.symbol,'entrez_gene_id':gene.entrez_gene_id,
            'n_yes':int(is_yes.sum()),'n_no':int((~is_yes).sum()),
            'delta_median_yes':float(np.median(v[is_yes])),'delta_median_no':float(np.median(v[~is_yes])),
            'delta_mean_yes_minus_no':interaction,'exact_permutation_p_change':pchange,
            'baseline_mean_yes_minus_no_log2Norm':baseline,'exact_permutation_p_baseline':pbaseline})
    df=pd.DataFrame(rows)
    df['bh_over_preselected_genes_change']=bh(df.exact_permutation_p_change.tolist())
    df['bh_over_preselected_genes_baseline']=bh(df.exact_permutation_p_baseline.tolist())
    df=df.sort_values('exact_permutation_p_change')
    df.to_csv(out/'hs_targeted_exact_permutation_7v7.csv',index=False)
    md=['# Adalimumab response vs nonresponse: exact 7-vs-7 exploratory comparison','',
        '- 14 independent patients; 28 lesional biopsies; response status obtained from GEO metadata.',
        '- A full 3432 exact partition-label permutation compares mean within-patient log2-normalized changes.',
        '- BH adjusted ONLY across preselected observed genes; *not* a transcriptome-wide FDR.',
        '- Same patients were used to choose treatment response labels retrospectively; no baseline prognostic validation.',
        '- The response criterion is unverified against full original clinical records; do not assert HiSCR50 from GEO alone.',
        '- Adalimumab (anti-TNF) is NOT secukinumab (anti-IL17A). Spatial GEO data are not patient-matched.','',
        '|Gene|Interaction (mean log2 delta YES-NO)|Exact p|Prespecified BH|Baseline mean log2 YES-NO|Baseline exact p|',
        '|---|---:|---:|---:|---:|---:|']
    for _,v in df.iterrows():md.append(f"|{v.gene}|{v.delta_mean_yes_minus_no:+.2f}|{v.exact_permutation_p_change:.4g}|{v.bh_over_preselected_genes_change:.4g}|{v.baseline_mean_yes_minus_no_log2Norm:+.2f}|{v.exact_permutation_p_baseline:.4g}|")
    (out/'HS_EXACT_PERMUTATION_REPORT.md').write_text('\n'.join(md)+'\n',encoding='utf8')
    print(df[['gene','delta_mean_yes_minus_no','exact_permutation_p_change','bh_over_preselected_genes_change']].to_string(index=False))

if __name__=='__main__':main()