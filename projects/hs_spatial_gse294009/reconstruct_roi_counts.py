#!/usr/bin/env python3
"""Public GEO GSE294009 GeoMx DCC->PKC gene-level targeted ROI audit.
No donor IDs supplied in GEO metadata; differences descriptive, NO p-values or DEG inference.
GEO 204 raw ROI != manuscript 156 QC-filtered ROI.
"""
import argparse, collections, csv, gzip, io, json, math, re, statistics, tarfile, hashlib
from pathlib import Path
from geo_file_probe import ACC,BASE,get,parse_series_matrix

TARGETS = ['CXCL13','MS4A1','CD79A','JCHAIN','IL17A','IL17F','CCL20','S100A8','S100A9','CXCL8','KRT17','CFH','COL4A5','KLK6','TNF','IL1B','MMP9','CD3D','CD68']

def quartile(vals,q=0.75):
    if not vals:return 0
    a=sorted(vals);p=(len(a)-1)*q;lo=int(p);hi=min(lo+1,len(a)-1)
    return a[lo]+(a[hi]-a[lo])*(p-lo)

def write_csv(path, rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='output/hs_spatial_gse294009');args=ap.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    audit={'accession':ACC,'analysis':'GEO GeoMx raw RTS_ID count mapping and per-ROI Q3 exploratory normalization', 'no_donor_inference':True}
    matrix,me=get(BASE+'/matrix/'+ACC+'_series_matrix.txt.gz',timeout=95)
    if matrix is None:raise RuntimeError('STOP: series matrix unavailable '+str(me))
    samples=parse_series_matrix(matrix)
    gsm={s['gsm']:s for s in samples}
    if len(gsm)!=204:raise RuntimeError('STOP: expected 204 distinct GEO GSM; observed '+str(len(gsm)))
    audit['series_sha256']=me.get('sha256')
    pkc,ps=get(BASE+'/suppl/GSE294009_Hs_R_NGS_WTA_v1.0.pkc.gz',limit=35_000_000,timeout=120)
    tar,ts=get(BASE+'/suppl/GSE294009_RAW.tar',limit=20_000_000,timeout=145)
    if pkc is None or tar is None:raise RuntimeError('STOP: missing PKC or DCC raw tar '+str((ps,ts)))
    audit['pkc_sha256']=ps['sha256'];audit['raw_tar_sha256']=ts['sha256']
    doc=json.loads(gzip.decompress(pkc))
    probes={}
    classes=collections.Counter()
    gene_probe_count=collections.Counter()
    for target in doc.get('Targets',[]):
        sym=str(target.get('DisplayName','')).strip();cls=str(target.get('CodeClass','')).strip();classes[cls]+=1
        for p in target.get('Probes',[]):
            code=str(p.get('RTS_ID','')).strip()
            if code:
                if code in probes and probes[code][0]!=sym:raise RuntimeError('Ambiguous RTS gene mapping: '+code)
                probes[code]=(sym,cls)
                gene_probe_count[sym]+=1
    if not probes:raise RuntimeError('STOP: PKC has no RTS mappings')
    audit.update({'pkc_target_count':len(doc.get('Targets',[])),'pkc_unique_RTS':len(probes),'pkc_CodeClass_counts':dict(classes), 'mapped_candidate_genes':[g for g in TARGETS if g in gene_probe_count],
                  'mapped_target_gene_probe_counts':{g:gene_probe_count[g] for g in TARGETS if g in gene_probe_count}})
    roi_rows=[];qc_rows=[];rt_status=collections.Counter();slides=collections.defaultdict(lambda:collections.Counter())
    all_rts=set();first_count_summary={};seen=set();sample_names=[];detected_targets=collections.Counter();bad_q3=[]
    with tarfile.open(fileobj=io.BytesIO(tar),mode='r:*') as tf:
        members=[m for m in tf.getmembers() if m.isfile() and '.dcc' in m.name.lower()]
        for i,mem in enumerate(members):
            match=re.search(r'GSM\d+',mem.name)
            if match is None: raise RuntimeError('DCC has no GSM: '+mem.name)
            accession=match.group()
            if accession not in gsm or accession in seen: raise RuntimeError('Missing/duplicate GSM: '+accession)
            seen.add(accession)
            b=tf.extractfile(mem).read()
            if mem.name.lower().endswith('.gz'):b=gzip.decompress(b)
            text=b.decode('utf8','replace')
            tag=re.search(r'<Code_Summary>\s*(.*?)\s*</Code_Summary>',text,re.S|re.I)
            if not tag:raise RuntimeError('Missing Code_Summary: '+accession)
            codes={}
            for line in tag.group(1).splitlines():
                items=line.strip().split(',')
                if len(items)!=2:continue
                key=items[0].strip();count=items[1].strip()
                if not re.fullmatch(r'RTS\d+',key) or not re.fullmatch(r'\d+',count):continue
                if key in codes:raise RuntimeError('DCC duplicate RTS: '+accession+' '+key)
                codes[key]=int(count)
            if not codes:raise RuntimeError('No numeric RTS counts in '+accession)
            unknown={c for c in codes if c not in probes};rt_status['unknown_probe_instances']+=len(unknown)
            rt_status['parsed_probe_instances']+=len(codes)
            if not first_count_summary:
                first_count_summary={'gsm':accession,'n_codes':len(codes),'n_matched':len(set(codes)&set(probes)),'first_keys':list(codes)[:12]}
            counts=collections.Counter()
            negative=[];endogenous=[]
            for code,v in codes.items():
                if code not in probes:continue
                sym,cls=probes[code]
                if cls.lower().startswith('endo'):
                    counts[sym]+=v
                elif 'negative' in cls.lower():negative.append(v)
            # Total genes in PKC not necessarily observed; include 0 for absent probes in geometric Q3.
            endogenes=[t['DisplayName'] for t in doc.get('Targets',[]) if str(t.get('CodeClass','')).lower().startswith('endo')]
            endocount=[counts.get(g,0) for g in endogenes]
            q3=quartile(endocount)
            q90=quartile(endocount,q=0.90)
            if q3<=0:bad_q3.append(accession)
            neg_med=statistics.median(negative) if negative else None
            # Exploratory WTA background estimate, NOT the paper's official LOQ.
            # Zeros are clamped to 1 before geometric statistics to keep the
            # transformation deterministic and conservative in low-count ROI.
            neg_log=[math.log(max(1.0,float(v))) for v in negative]
            neg_geom=math.exp(statistics.mean(neg_log)) if len(neg_log)>=2 else None
            neg_geosd=math.exp(statistics.stdev(neg_log)) if len(neg_log)>=2 else None
            neg_loq=max(2.0,neg_geom*neg_geosd**2) if neg_geom is not None else None
            q3_vs_neg=(q3/neg_geom) if neg_geom and q3>0 else None
            title=gsm[accession].get('title','')
            condition=gsm[accession].get('characteristics_1','').replace('disease:','').strip()
            stage='Early' if condition.startswith('Early_HS_') else 'Late' if condition.startswith('Late_HS_') else 'Control' if condition.startswith('Control_') else 'Other'
            tissue=(condition.replace('Early_HS_','').replace('Late_HS_','').replace('Control_','')) if stage!='Other' else condition
            slide=title.split('|')[0].strip()
            if stage!='Other' and stage!='Control':slides[slide][stage]+=1
            qc={'gsm':accession,'slide_label_not_donor':slide,'roi_title':title,'condition':condition,'stage':stage,'tissue_compartment':tissue,'n_codes':len(codes),'total_RTS_counts':sum(codes.values()),'n_endogenous_gene_detected':sum(x>0 for x in endocount),'Q3_endogenous_raw':round(q3,3),'q3_usable':int(q3>0),'negative_probe_median':neg_med,'n_negative_probes':len(negative),
                'negative_probe_geomean_clamp1':round(neg_geom,6) if neg_geom is not None else '',
                'negative_probe_geosd_clamp1':round(neg_geosd,6) if neg_geosd is not None else '',
                'exploratory_LOQ_clamp1_geomean_geosd2':round(neg_loq,6) if neg_loq is not None else '',
                'Q3_div_negative_geomean':round(q3_vs_neg,6) if q3_vs_neg is not None else '',
                'Q90_endogenous_raw':round(q90,6)}
            qc_rows.append(qc)
            for gene in TARGETS:
                if gene not in gene_probe_count:continue
                raw=counts.get(gene,0)
                # Scale each ROI to reference Q3 to make intensity comparisons descriptive only.
                # Raw Q3-normalized count as a ratio; no global mean need be guessed.
                roi_rows.append({'gsm':accession,'slide_label_not_donor':slide,'stage':stage,'condition':condition,'compartment':tissue,'gene':gene,'raw_count':raw,'gene_probe_count':gene_probe_count[gene],
                                 'q3_scaled_ratio':round(raw/q3,6) if q3>0 else None,
                                 'q90_scaled_ratio':round(raw/q90,6) if q90>0 else None,
                                 'log2_ratio_plus1':round(math.log2(1+raw/q3),6) if q3>0 else None,
                                 'detected':int(raw>0),
                                 'above_exploratory_LOQ':int(raw>=neg_loq) if neg_loq is not None and gene_probe_count[gene]==1 else ''})
    if seen!=set(gsm):raise RuntimeError('GSM correspondence not 204/204: '+str(len(seen)))
    all_slides=collections.Counter(x['slide_label_not_donor'] for x in qc_rows)
    # descriptive only: groups of ROI, no statistical replicates
    groups=collections.defaultdict(list)
    for row in roi_rows:
        groups[(row['stage'],row['compartment'],row['gene'])].append(row)
    summaries=[]
    for (stage,comp,g),rows in sorted(groups.items()):
        vals=[float(x['q3_scaled_ratio']) for x in rows if x['q3_scaled_ratio'] is not None]
        if not vals:continue
        summaries.append({'stage':stage,'compartment':comp,'gene':g,'n_roi':len(vals),'n_total_roi_in_group':len(rows),'n_excluded_zero_q3':len(rows)-len(vals),'median_q3_scaled':round(statistics.median(vals),6),'mean_q3_scaled':round(statistics.mean(vals),6),'detect_fraction':round(sum(x['detected'] for x in rows)/len(rows),4)})
    contrasts=[]
    lookup={(r['stage'],r['compartment'],r['gene']):r for r in summaries}
    compartments=sorted(set(x['compartment'] for x in summaries))
    for comp in compartments:
        for gene in TARGETS:
            early=lookup.get(('Early',comp,gene));late=lookup.get(('Late',comp,gene))
            if early is None or late is None:continue
            e=early['median_q3_scaled'];l=late['median_q3_scaled']
            contrasts.append({'gene':gene,'compartment':comp,'early_n_roi':early['n_roi'],'late_n_roi':late['n_roi'],'early_median_q3_scaled':e,'late_median_q3_scaled':l,'late_vs_early_log2_ratio_plus01':round(math.log2((l+.1)/(e+.1)),4),'no_pvalue_due_to_missing_donor_ID':True})
    write_csv(out/'roi_gene_targeted_q3.csv',roi_rows)
    write_csv(out/'roi_qc.csv',qc_rows)
    write_csv(out/'stage_compartment_gene_descriptive.csv',summaries)
    write_csv(out/'early_late_shared_compartment_descriptive.csv',contrasts)
    write_csv(out/'slide_stage_counts.csv',[{'slide_label_not_donor':slide,'early_ROI':v.get('Early',0),'late_ROI':v.get('Late',0),'is_mixed':bool(v.get('Early',0) and v.get('Late',0))} for slide,v in slides.items()])
    write_csv(out/'geo_roi_metadata.csv',samples)
    audit.update({'GSMs_matched':len(seen),'ROI_count':len(qc_rows),'slide_labels_n':len(all_slides),'slide_label_counts':dict(all_slides),'ROI_group_counts':dict(collections.Counter(x['stage'] for x in qc_rows)),'shared_compartment_contrasts_n':len(contrasts),'q3_unusable_ROI_count':len(bad_q3),'q3_unusable_GSMs':bad_q3[:50],'exploratory_loq_method':{'definition':'max(2, negative_probe_geoMean*geoSD^2) from raw DCC negative probes, clamp counts below 1 before log, sample SD; signal baseline only, not verified author LOQ; use only genes with exactly one RTS probe','missing_negative_probes_gsm':sum(1 for r in qc_rows if not r['exploratory_LOQ_clamp1_geomean_geosd2'])},
                  'quantification':{'definition':'raw DCC RTS counts summed by PKC DisplayName (endogenous targets); per-ROI ratio = gene raw counts / Q3 of all endogenous gene counts; no donor adjustment','unknown_probe_instances':rt_status['unknown_probe_instances'],'parsed_probe_instances':rt_status['parsed_probe_instances']},'first_dcc':first_count_summary})
    audit['status']='GO_ROI_TARGETED_DESCRIPTIVE_ONLY_HOLD_DONOR_LEVEL_INFERENCE'
    audit['clinical_limitations']=['No explicit 15-patient donor map in GEO series/sample metadata','Paper uses 156 ROI whereas GEO hosts 204; QC-filtered ROI list not reconstructed','ROI from same sample/slide not independent patients','No secukinumab treatment outcomes or patient-level pairing to GSE213761','Cannot estimate valid FDR for disease stage comparisons without donor group ID']
    (out/'RECONSTRUCTION_AUDIT.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf8')
    report=['# GSE294009 GeoMx DCC / PKC reconstruction','',f"**Decision:** {audit['status']}", '',f"- Original GEO DCC GSM matched: **{len(seen)}/204**",f"- PKC targets: {audit['pkc_target_count']}; RTS probes: {audit['pkc_unique_RTS']}",f"- Distinct slide labels (NOT donors): {len(all_slides)}",f"- ROI group strata: {audit['ROI_group_counts']}",f"- Q3 unusable ROIs: {len(bad_q3)} (raw counts retained, normalized intensity NA)",f"- Shared early/late compartment gene rows: {len(contrasts)}",'','## QC / statistical constraints']+['- '+s for s in audit['clinical_limitations']]+['','## Example early-vs-late descriptive candidates']
    for gene in ['CXCL13','CCL20','IL17A','S100A8','S100A9','KRT17','MS4A1','CD79A']:
        for x in contrasts:
            if x['gene']==gene and x['compartment'] in ['Dermis_Deep_Inflammation','Epithelium_Epidermis']:
                report.append(f"- {gene} / {x['compartment']}: Early {x['early_median_q3_scaled']}, Late {x['late_median_q3_scaled']} (n={x['early_n_roi']}/{x['late_n_roi']}); *no inferential p-value*")
    (out/'RECONSTRUCTION_REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf8')
    print(json.dumps({'status':audit['status'],'GSMs_matched':len(seen),'targets':audit['pkc_target_count'],'n_slides':len(all_slides),'n_contrasts':len(contrasts),'unknown_probe_instances':rt_status['unknown_probe_instances']},ensure_ascii=False))

if __name__=='__main__':main()