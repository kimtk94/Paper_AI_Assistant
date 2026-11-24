#!/usr/bin/env python3
"""Pre-registered descriptive sensitivity for GSE294009 GeoMx ROI expression.
No hypothesis tests: donor IDs and publication QC ROI membership not recovered.
"""
import argparse,csv,json,math,statistics,collections
from pathlib import Path

def read(path):
    with open(path,newline='',encoding='utf8') as f:return list(csv.DictReader(f))
def save(path,rows):
    if not rows:return
    with open(path,'w',newline='',encoding='utf8') as f:
        keys=list(dict.fromkeys(k for row in rows for k in row))
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def median_stage(group,stage):
    vals=[float(x['q3_scaled_ratio']) for x in group if x['stage']==stage and x['q3_scaled_ratio'] not in ('',None,'nan','None')]
    return (len(vals),statistics.median(vals)) if vals else (0,None)
def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='output/hs_spatial_gse294009');a=p.parse_args();out=Path(a.out)
    roi=read(out/'roi_gene_targeted_q3.csv'); qc=read(out/'roi_qc.csv')
    qc_map={x['gsm']:x for x in qc};assert len(qc_map)==204
    slide_stage=collections.defaultdict(lambda:collections.Counter())
    for x in qc:slide_stage[x['slide_label_not_donor']][x['stage']]+=1
    # Only explicitly present GEO slide label spanning both HS stages (not donor ID).
    mixed=[k for k,c in slide_stage.items() if c['Early'] and c['Late']]
    if not mixed:raise RuntimeError('No slide with both Early and Late HS; stop comparisons')
    selected={'all_q3_positive':lambda x:float(qc_map[x['gsm']]['Q3_endogenous_raw'])>0,
              'q3_at_least_five_sensitivity':lambda x:float(qc_map[x['gsm']]['Q3_endogenous_raw'])>=5,
              'within_mixed_slide_descriptive':lambda x:x['slide_label_not_donor'] in mixed and float(qc_map[x['gsm']]['Q3_endogenous_raw'])>0}
    comps=sorted({x['compartment'] for x in roi if x['stage'] in ('Early','Late')})
    targets=sorted({x['gene'] for x in roi})
    output=[]
    for comp in comps:
        for gene in targets:
            r=[x for x in roi if x['compartment']==comp and x['gene']==gene and x['stage'] in ('Early','Late')]
            if not r:continue
            for label,check in selected.items():
                rr=[x for x in r if check(x)]
                n1,med1=median_stage(rr,'Early');n2,med2=median_stage(rr,'Late')
                if not n1 or not n2:continue
                # Add a fixed 0.1 ratio pseudocount only for descriptive log fold summary, not count-based DE.
                log2ratio=math.log2((med2+0.1)/(med1+0.1))
                output.append({'gene':gene,'compartment':comp,'sensitivity':label,'n_early_roi':n1,'n_late_roi':n2,'early_median_gene_raw_over_q3':round(med1,6),'late_median_gene_raw_over_q3':round(med2,6),'descriptive_log2_late_early_with_point_one':round(log2ratio,4),'independent_donors':'UNKNOWN','statistical_p_value':'NOT_COMPUTED'})
    save(out/'slide_batch_sensitivity.csv',output)
    save(out/'slide_stage_distribution.csv',[{'slide_label_not_patient':k,**dict(v)} for k,v in sorted(slide_stage.items())])
    gene_list=['CXCL13','CCL20','S100A8','S100A9','MS4A1','CD79A','IL17A','KRT17']
    targets_important=[x for x in output if x['gene'] in gene_list and x['compartment'] in ['Epithelium_Epidermis','Dermis_Deep_Inflammation']]
    report=['# HS GSE294009 slide/batch/Q3 sensitivity audit','',
        f'- Stage-mixed slide labels: {mixed} (NOT patients)',
        '- ROI stage and slide are not balanced: several slides hold exclusively early or late ROI.',
        '- Original paper 156 QC ROIs remain unidentified among 204 GEO registrations.',
        '- No donor-level inference, no p-values, no false independent ROI sample sizes.',
        '- Q3>=5 is an exploratory normalization sensitivity threshold, not official GeoMx QC.',
        '','## Candidate descriptive sensitivity (selected same-compartment ROI)','',
        '|Gene|Compartment|Analysis|Early n / Late n|Early median|Late median|Descriptive log2 ratio|',
        '|---|---|---|---:|---:|---:|---:|']
    for x in targets_important:
        report.append(f"|{x['gene']}|{x['compartment']}|{x['sensitivity']}|{x['n_early_roi']} / {x['n_late_roi']}|{x['early_median_gene_raw_over_q3']}|{x['late_median_gene_raw_over_q3']}|{x['descriptive_log2_late_early_with_point_one']}|")
    (out/'SLIDE_BATCH_SENSITIVITY.md').write_text('\n'.join(report)+'\n',encoding='utf8')
    result={'status':'DESCRIPTIVE_ONLY','n_slides':len(slide_stage),'mixed_slides':mixed,'all_ROI':len(qc),'q3_positive_ROI':sum(float(x['Q3_endogenous_raw'])>0 for x in qc),'q3_at_least_five_ROI':sum(float(x['Q3_endogenous_raw'])>=5 for x in qc),'n_sensitivity_rows':len(output)}
    (out/'SLIDE_BATCH_AUDIT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()