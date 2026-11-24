#!/usr/bin/env python3
"""Exploratory independent verification of IPF vs control plasma protein abundance.

Uses public PXD010965 Proteins.csv with two merged header rows and a real third
field-name row. Does NOT test progression/FVC or claim to rerun Progenesis QI.
"""
import argparse, csv, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

EXPECTED_SHA256='1df37c4f3a2e130bd02ea87dc5da47ba1a6701efe24c657b358cc721bbd3f334'

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True)
    parser.add_argument('--out',default='output/ipf_pxd010965_reanalysis')
    a=parser.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    raw=Path(a.input).read_bytes()
    checksum=hashlib.sha256(raw).hexdigest()
    if checksum!=EXPECTED_SHA256:
        raise ValueError('Source SHA256 changed; review original data before reanalysis')
    rows=list(csv.reader(raw.decode('utf-8-sig').splitlines()))
    if len(rows)<5:raise ValueError('Too few protein rows')
    width=len(rows[0]);assert width==49, f'Expected 49 columns, got {width}'
    if any(len(r)!=width for r in rows):raise ValueError('Ragged CSV')
    headers=rows[2]
    if headers[:2]!=['Accession','Peptide count'] or rows[0][12]!='Raw abundance':
        raise ValueError('Unexpected multi-line headers')
    starts=[(i,x) for i,x in enumerate(rows[1]) if x in ('CONTROL','CASE')]
    if starts!=[(12,'CONTROL'),(32,'CASE')]:
        raise ValueError('Mismatch in explicitly published group header positions')
    ids=headers[12:]
    if len(ids)!=37 or len(set(ids))!=37:raise ValueError('Not 37 unique assay labels')
    meta={"CONTROL":ids[:20],"CASE":ids[20:]}
    dat=pd.DataFrame(rows[3:],columns=[*headers[:12],*ids])
    features=dat['Accession']
    if features.isna().any() or features.duplicated().any(): raise ValueError('Nonunique protein accession rows')
    abundance=dat.loc[:,ids].apply(pd.to_numeric,errors='coerce')
    bad=(abundance<=0).sum().sum()
    missing=abundance.isna().sum().sum()
    if bad or missing:raise ValueError(f'Invalid raw abundances: nonpositive={bad}, missing={missing}')
    log2=np.log2(abundance.to_numpy(dtype=float))
    ctrl=log2[:,:20];case=log2[:,20:]
    # Sensitivity only: per-sample median centering is NOT the published Progenesis normalization.
    sample_med=np.median(log2,axis=0)
    centered=log2-sample_med[None,:]
    centered_effect=centered[:,20:].mean(axis=1)-centered[:,:20].mean(axis=1)
    centered_p=ttest_ind(centered[:,20:],centered[:,:20],axis=1,equal_var=False).pvalue
    centered_q=multipletests(centered_p,method='fdr_bh')[1]
    sample_median_group_p=float(ttest_ind(sample_med[:20],sample_med[20:],equal_var=False).pvalue)
    mean_diff=case.mean(axis=1)-ctrl.mean(axis=1)
    pv=ttest_ind(case,ctrl,axis=1,equal_var=False).pvalue
    adj=multipletests(pv,method='fdr_bh')[1]
    med_diff=np.median(case,axis=1)-np.median(ctrl,axis=1)
    mwu=np.array([mannwhitneyu(case[i,:],ctrl[i,:],alternative='two-sided').pvalue for i in range(len(dat))])
    mwu_adj=multipletests(mwu,method='fdr_bh')[1]
    q_original=pd.to_numeric(dat['q Value'],errors='coerce').to_numpy()
    p_original=pd.to_numeric(dat['Anova (p)'],errors='coerce').to_numpy()
    gene=dat['Description'].str.extract(r'GN=([^ ]+)')[0].fillna('')
    res=pd.DataFrame({'accession':features,'gene':gene,'description':dat['Description'],
        'original_anova_p':p_original,'original_q':q_original,
        'log2fc_case_minus_control':mean_diff,'median_log2fc_case_minus_control':med_diff,
        'welch_p':pv,'welch_fdr':adj,'mwu_p':mwu,'mwu_fdr':mwu_adj,
        'sample_median_centered_log2fc':centered_effect,
        'sample_median_centered_welch_p':centered_p,
        'sample_median_centered_welch_fdr':centered_q})
    res=res.sort_values('welch_fdr').reset_index(drop=True)
    res.to_csv(out/'protein_reanalysis.csv',index=False)
    qc=[]
    for i,name in enumerate(ids):
        col=log2[:,i]
        qc.append({'assay_id':name,'group':'CONTROL' if i<20 else 'CASE','median_log2_abundance':float(np.median(col)),
            'mean_log2_abundance':float(col.mean()),'n_proteins':len(col)})
    pd.DataFrame(qc).to_csv(out/'sample_qc.csv',index=False)
    # Is the reported study q-value rank broadly reproducible with a different test?
    ok=np.isfinite(q_original)&np.isfinite(adj)&(q_original>0)
    rho=float(spearmanr(-np.log10(q_original[ok]),-np.log10(np.maximum(adj[ok],1e-300))).statistic)
    # Leave-one-subject-out direction stability; inferential p-values are not recomputed.
    target=['HPR','CFH','CPN1','FCN3','C3','SERPINA3','GPX3']
    stability=[]
    for g in target:
        sub=res[res['gene']==g]
        if len(sub)!=1:continue
        accession=sub.iloc[0]['accession']
        j=int(np.where(features==accession)[0][0])
        vals=[case[j,:].mean()-np.delete(ctrl[j,:],k).mean() for k in range(20)] + [np.delete(case[j,:],k).mean()-ctrl[j,:].mean() for k in range(17)]
        sign=float(np.sign(mean_diff[j]));stable=sum(np.sign(v)==sign for v in vals)
        stability.append({'gene':g,'accession':accession,'mean_log2FC':float(mean_diff[j]),
            'welch_FDR':float(adj[j]),'published_q':float(q_original[j]),'same_direction_loo':int(stable),'total_loo':37})
    pd.DataFrame(stability).to_csv(out/'candidate_loo.csv',index=False)
    # Note: participant-level clinical covariates not supplied, raw abundance is not normalized value.
    summary={'dataset':'PXD010965','source_sha256':checksum,'statistical_status':'EXPLORATORY_RAW_ABUNDANCE',
        'n_proteins':len(dat),'controls':20,'cases':17,'n_sample_columns':len(ids),
        'complete_numeric_abundances':int(abundance.notna().sum().sum()),
        'welch_fdr_005':int((adj<.05).sum()),'welch_fdr_005_abs_log2fc_1':int(((adj<.05)&(np.abs(mean_diff)>1)).sum()),
        'mwu_fdr_005':int((mwu_adj<.05).sum()),'published_q_005':int((q_original<.05).sum()),
        'median_centered_welch_fdr_005':int((centered_q<.05).sum()),
        'median_centered_welch_fdr_005_abs_log2fc_1':int(((centered_q<.05)&(np.abs(centered_effect)>1)).sum()),
        'raw_vs_centered_fdr_overlap':int(((adj<.05)&(centered_q<.05)).sum()),
        'sample_median_group_welch_p':sample_median_group_p,
        'spearman_published_q_vs_new_FDR':rho,
        'group_sample_median_log2':{'CONTROL':float(np.median(np.median(ctrl,axis=0))),
                                    'CASE':float(np.median(np.median(case,axis=0)))},
        'top_candidates':stability,
        'warnings':['Raw abundance from Progenesis QI may include study-specific normalization and batch effects not reconstructed.',
         'Published ANOVA/q values are not the same statistical model as an independent Welch test on log2 raw abundance.',
         'Sample medians differ strongly between groups in the raw table; raw-abundance significance may reflect systematic intensity shifts.',
         'Per-sample median centering is only a sensitivity analysis and not equivalent to the authors published normalization.',
         'FVC progression labels, treatment covariates, smoking status, age/sex are not present in this file.',
         'Candidate lists cannot be used for clinical diagnosis or independent validation without new cohort.']}
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# IPF PXD010965 public plasma proteomics reanalysis', '',
      '- Source: 2020 FinnishIPF vs controls, no Korean subjects or progression labels',
      f"- {len(dat)} protein rows, {len(ids)} assays (control 20, case 17)",
      f"- Exploratory log2 Welch FDR<0.05: {summary['welch_fdr_005']}; FDR<0.05 and abs(log2FC)>1: {summary['welch_fdr_005_abs_log2fc_1']}",
      f"- Nonparametric MWU FDR<0.05: {summary['mwu_fdr_005']}; published q<0.05: {summary['published_q_005']}",
      f"- Median-centered sensitivity Welch FDR<0.05: {summary['median_centered_welch_fdr_005']} (|log2FC|>1: {summary['median_centered_welch_fdr_005_abs_log2fc_1']}); overlap raw={summary['raw_vs_centered_fdr_overlap']}",
      f"- Sample median raw log2 intensity difference Welch p={summary['sample_median_group_welch_p']:.3g}; caution batch/normalization.",
      f"- Spearman ranking concordance (published q vs new Welch FDR): {rho:.3f}",
      '','## Targeted protein checks']
    lines.extend([f"- {z['gene']}: log2FC={z['mean_log2FC']:.3f}; Welch FDR={z['welch_FDR']:.3g}; original q={z['published_q']:.3g}; LOO sign={z['same_direction_loo']}/37" for z in stability])
    lines.extend(['', '## Limitations',*['- '+z for z in summary['warnings']]])
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()