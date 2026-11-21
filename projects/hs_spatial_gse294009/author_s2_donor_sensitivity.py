#!/usr/bin/env python3
"""GSE294009 *descriptive* patient-aware check using official NIHMS2203138 Suppl Table S1/S2.
GSM->donor links explicit, never infer donor from slide. Exclude external 24 ROIs without DCC.
All stages from S1 with S2 mismatch flagged; comparisons are exploratory (no DEG p/FDR).
"""
import argparse,csv,collections,statistics,math,json
from pathlib import Path

TARGETS={'CXCL13','CCL20','S100A8','S100A9','IL17A','MS4A1','CD79A','JCHAIN','KRT17','CXCL8','CFH'}

def read(p):
    with open(p,newline='',encoding='utf8') as f:return list(csv.DictReader(f))

def save(path,rows):
    if not rows:return
    cols=list(dict.fromkeys(k for row in rows for k in row))
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)

def fval(v):
    try:
        q=float(v)
        return q if math.isfinite(q) else None
    except (TypeError,ValueError):return None

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',default='output/hs_spatial_gse294009');p.add_argument('--crosswalk',default='projects/hs_spatial_gse294009/author_s2_geo_crosswalk.csv');p.add_argument('--out',default='output/hs_spatial_gse294009');a=p.parse_args()
    folder=Path(a.out);folder.mkdir(parents=True,exist_ok=True)
    m=read(a.crosswalk);expr=read(Path(a.input)/'roi_gene_targeted_q3.csv')
    meta={x['gsm']:x for x in m}
    if len(meta)!=146 or len(m)!=146:raise RuntimeError('STOP original S2 GEO crosswalk not unique 146')
    if len({v['patient'] for v in m})!=16:raise RuntimeError('STOP S2 patient count mismatch')
    disc=[x for x in m if x['discordant_stage']=='1']
    if len(disc)!=4 or {x['patient'] for x in disc}!={'HS12'}:raise RuntimeError('STOP unexpected clinical S1/S2 stage conflict')
    gsm=set(x['gsm'] for x in expr)
    if not set(meta)<=gsm:raise RuntimeError('STOP S2 GSM missing in reconstructed DCC ROI')
    observed={x['gsm'] for x in expr if x['gene']==expr[0]['gene']}
    if len(observed)!=204:raise RuntimeError('STOP expression does not have all 204 ROI')
    count=collections.Counter();pairs=collections.defaultdict(list)
    for item in expr:
        subject=meta.get(item['gsm'])
        if not subject:continue
        val=fval(item.get('q3_scaled_ratio'))
        if val is None:continue
        gene=item['gene']; comp=item['compartment']
        rec={**subject,'gene':gene,'compartment':comp,'ratio':val}
        pairs[(subject['patient'],gene,comp)].append(rec)
        count['roi_gene_values']+=1
    donor=[]
    for (patient,gene,comp),xs in pairs.items():
        src=xs[0]
        for x in xs:
            if (x['institution'],x['clinical_stage_from_s1'],x['author_stage_from_roi'])!=(src['institution'],src['clinical_stage_from_s1'],src['author_stage_from_roi']):
                raise RuntimeError('Patient nonconstant stage or site: '+patient)
        donor.append({'patient_id':patient,'gene':gene,'compartment':comp,'n_roi':len(xs),'donor_median_q3_ratio':statistics.median(x['ratio'] for x in xs),
            's1_stage':src['clinical_stage_from_s1'],'s2_roi_stage':src['author_stage_from_roi'], 'site':src['institution'],'discordant':src['discordant_stage']})
    assert len({r['patient_id'] for r in donor})==16
    # Avoid sample-level GEO info in artifacts; save aggregate only.
    tests={'clinical_all':lambda d: True,'clinical_exclude_discordant':lambda d:d['discordant']!='1',
          'same_site_yonsei_exclude_discordant':lambda d:d['site']=='Yonsei University' and d['discordant']!='1',
          'geo_roi_label_all':lambda d:True}
    results=[]
    for t,fn in tests.items():
        stagekey='s2_roi_stage' if t=='geo_roi_label_all' else 's1_stage'
        for g,comp in sorted({(d['gene'],d['compartment']) for d in donor}):
            selected=[d for d in donor if d['gene']==g and d['compartment']==comp and fn(d) and d[stagekey] in ('Early','Late')]
            early=[float(d['donor_median_q3_ratio']) for d in selected if d[stagekey]=='Early']
            late=[float(d['donor_median_q3_ratio']) for d in selected if d[stagekey]=='Late']
            if not early or not late:continue
            m1=statistics.median(early);m2=statistics.median(late)
            results.append({'gene':g,'compartment':comp,'sensitivity':t,
                'early_donors':len(early),'late_donors':len(late),'early_donor_median_q3ratio':round(m1,6),
                'late_donor_median_q3ratio':round(m2,6), 'descriptive_log2_late_early_pseudo_0p1':round(math.log2((m2+.1)/(m1+.1)),4),
                'p_value':'NOT_COMPUTED','FDR':'NOT_COMPUTED'})
    # Key population stratum: source may influence both stage & gene expression.
    groupcounts=collections.Counter((x['clinical_stage_from_s1'],x['institution'],x['slide_preparation']) for x in m)
    printrows=[{'clinical_stage':k[0],'institution':k[1],'slide_prep':k[2],'n_roi':v,
       'n_donors':len({r['patient'] for r in m if (r['clinical_stage_from_s1'],r['institution'],r['slide_preparation'])==k})} for k,v in sorted(groupcounts.items())]
    result={'accession':'GSE294009','source':'Official NIHMS2203138 Supplementary Tables 1–3, Europe PMC / PMC13552446',
       'n_geo_s2_rows':len(m),'n_explicit_donor_ids':16,'n_HS_donors':15,
       'n_HS_geo_ROIs':sum(x['patient'].startswith('HS') for x in m),'n_AC_geo_ROIs':sum(x['patient'].startswith('AC') for x in m),
       'n_not_in_s2_of_204':204-len(m),'n_prior_study_HS_ROIs':24,
       'n_paper_HS_ROIs':sum(x['patient'].startswith('HS') for x in m)+24,
       'S1_S2_conflict':{'donor':'HS12','roi_ids':[x['gsm'] for x in disc],'S1':'Late Hurley III','S2':'Early ROI classification'},
       'donor_site_counts':printrows,
       'n_donor_gene_compartment_records':len(donor),'n_descriptive_contrasts':len(results),
       'status':'DESCRIPTIVE_DONOR_AWARE_WITH_BATCH_CONFOUNDING',
       'warnings':['No independent ROIs or ROI-level p/FDR.',
                   'S1/S2 HS12 stage discrepancy requires source correction; exclude conflicted donor in primary sensitivity.',
                   'GEO GSE294009 lacks 24 previously published HS ROIs included in S2; no exact full paper quant replication.',
                   'Multiple ROIs per donor; clinical data and study site confounded with stage.',
                   'Author normalization/LOQ and unknown exclusions cannot be equated with exploratory Q3 ratio.',
                   'Cross-study HS treatment response GSE213761 is different donor and Adalimumab; no joint prediction.']}
    save(folder/'AUTHOR_S2_SITE_DONOR_COUNTS.csv',printrows)
    save(folder/'AUTHOR_S2_DONOR_AGGREGATED_CONTRASTS.csv',results)
    (folder/'AUTHOR_S2_DONOR_AUDIT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    key=[r for r in results if r['gene'] in TARGETS and r['compartment'] in ('Epithelium_Epidermis','Dermis_Deep_Inflammation')]
    lines=['# GSE294009 author S2 patient-aware targeted GeoMx audit','','**Interpretation: descriptive only.** Explicit original patient IDs from author-published S2; see S1/S2 HS12 disagreement.','',f"- S2 listed GEO ROI = {len(m)}; HS=132, AC=14, plus historical external HS=24; manuscript HS ROI total=156.",
      '- The full GEO has 204 registrations; S2 lists only 146 of those plus 24 historic datasets.',
      '- Site/slide preparations and disease stage are confounded; no adjusted differential expression is estimated.',
      '', '## Same-compartment gene signals, donor medians','', '|Gene|Compartment|Sensitivity|Early donors / Late donors|Descriptive late-vs-early log₂ ratio|','|---|---|---|---:|---:|']
    for r in key:lines.append(f"|{r['gene']}|{r['compartment']}|{r['sensitivity']}|{r['early_donors']} / {r['late_donors']}|{r['descriptive_log2_late_early_pseudo_0p1']}|")
    lines+=['','## Required gates','- Clarify S1/S2 HS12 contradictory stage label with source study and re-verify external 24 prior ROIs.','- Inference must adjust site/slide/ROI QC and use independent donors, not the GEO 204 row count.','- No claim that donor RNA changes predict Secukinumab or Adalimumab response.']
    (folder/'AUTHOR_S2_DONOR_REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':result['status'],'n_donors':result['n_explicit_donor_ids'],'HS_S2':result['n_HS_geo_ROIs'],'external':result['n_prior_study_HS_ROIs'],'n_contrasts':len(results)},ensure_ascii=False))

if __name__=='__main__':main()