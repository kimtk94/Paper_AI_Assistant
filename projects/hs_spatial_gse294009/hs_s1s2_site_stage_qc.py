#!/usr/bin/env python3
"""S1/S2 xlsx versus GSE294009 crosswalk audit; no inferential statistics."""
import collections, csv, hashlib, json, re, xml.etree.ElementTree as ET, zipfile
from pathlib import Path

ROOT=Path('output/hs_spatial_supp')
CROSS=Path('projects/hs_spatial_gse294009/author_s2_geo_crosswalk.csv')
ORIGINAL=ROOT/'NIHMS2203138-supplement-Supplementary_Table_1-3.xlsx'
EXPECTED_SHA='fb65bdc9fd578eb4746c428faf9c7d4d4c37d6ef248205fcf36ae0df5317de9a'
NS={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
GSM=re.compile(r'GSM\d{7,9}')

def source_sheet_rows(path):
    with zipfile.ZipFile(path) as z:
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            root=ET.fromstring(z.read('xl/sharedStrings.xml'))
            shared=[''.join(q.text or '' for q in el.findall('.//x:t',NS))
                    for el in root.findall('x:si',NS)]
        sheets=sorted((p for p in z.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml',p)),
                      key=lambda s:int(re.search(r'sheet(\d+)',s).group(1)))
        for name in sheets:
            root=ET.fromstring(z.read(name))
            rows=[]
            for row in root.findall('.//x:sheetData/x:row',NS):
                values=[]
                for cell in row.findall('x:c',NS):
                    v=cell.find('x:v',NS)
                    if cell.get('t')=='s' and v is not None:
                        val=shared[int(v.text)]
                    elif cell.get('t')=='inlineStr':
                        val=''.join(k.text or '' for k in cell.findall('.//x:t',NS))
                    else:
                        val=v.text if v is not None and v.text else ''
                    if val:values.append(str(val))
                if values:rows.append(values)
            yield name,rows

def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    with CROSS.open(encoding='utf-8-sig',newline='') as f:
        records=list(csv.DictReader(f))
    gsm={r['gsm']:r for r in records}
    hs=[r for r in records if r['patient'].startswith('HS')]
    ac=[r for r in records if r['patient'].startswith('AC')]
    if (len(records),len(gsm),len(hs),len(ac))!=(146,146,132,14):
        raise ValueError('STOP: S2 crosswalk count or GSM uniqueness drift')
    donors=collections.defaultdict(list)
    for r in hs:donors[r['patient']].append(r)
    if len(donors)!=15 or len({x['patient'] for x in ac})!=1:
        raise ValueError('STOP: patient cardinality drift')
    for p,rs in donors.items():
        for key in ['institution','slide_preparation','clinical_stage_from_s1','author_stage_from_roi']:
            if len({r[key] for r in rs})!=1:raise ValueError('STOP donor inconsistent: '+p+'/'+key)
    disc=[r for r in hs if r['clinical_stage_from_s1']!=r['author_stage_from_roi']]
    if len(disc)!=4 or {r['patient'] for r in disc}!={'HS12'}:
        raise ValueError('STOP: unexpected S1/S2 discordance')
    scenarios={
      'S2_all':(hs,'author_stage_from_roi'),
      'S2_exclude_HS12':([r for r in hs if r['patient']!='HS12'],'author_stage_from_roi'),
      'Yonsei_FFPE_exclude_HS12':([r for r in hs if r['institution']=='Yonsei University' and r['patient']!='HS12'],'author_stage_from_roi'),
      'S1_all':(hs,'clinical_stage_from_s1')
    }
    counts={}
    for key,(rs,stage_key) in scenarios.items():
        counts[key]={s:{'roi':sum(r[stage_key]==s for r in rs),
                        'donors':len({r['patient'] for r in rs if r[stage_key]==s})}
                     for s in ('Early','Late')}
    site={}
    for r in hs:
        key=(r['institution'],r['slide_preparation'],r['author_stage_from_roi'])
        v=site.setdefault(key,{'roi':0,'donors':set()})
        v['roi']+=1;v['donors'].add(r['patient'])
    site_rows=[{'site':a,'preparation':b,'stage':c,'roi':v['roi'],'donors':len(v['donors'])}
               for (a,b,c),v in sorted(site.items())]
    comps=collections.defaultdict(lambda:collections.defaultdict(set))
    for r in scenarios['Yonsei_FFPE_exclude_HS12'][0]:
        comp=r['author_roi_classification'].split('_HS_',1)[-1]
        comps[comp][r['author_stage_from_roi']].add(r['patient'])
    coverage=[{'compartment':k,'early_donors':len(v.get('Early',set())),
               'late_donors':len(v.get('Late',set())),
               'status':'DESCRIPTIVE_ONLY' if v.get('Early') and v.get('Late') else 'NO_BOTH_STAGES'}
              for k,v in sorted(comps.items())]
    with (ROOT/'SITE_MATCHED_COMPARTMENT_COVERAGE.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['compartment','early_donors','late_donors','status'])
        w.writeheader();w.writerows(coverage)
    if not ORIGINAL.is_file():raise ValueError('STOP: official XLSX absent from public supplement run')
    original_sha=hashlib.sha256(ORIGINAL.read_bytes()).hexdigest()
    if original_sha!=EXPECTED_SHA:raise ValueError('STOP: published workbook digest changed: '+original_sha)
    found=collections.defaultdict(list);inventory=[];s2_no_gsm=[]
    for sheet,rows in source_sheet_rows(ORIGINAL):
        sheet_gsm=set()
        untagged=[]
        for index,row in enumerate(rows,1):
            joined=' | '.join(row)
            row_gsm=set(GSM.findall(joined))
            if not row_gsm:untagged.append({'row':index,'cells':row[:8]})
            for code in row_gsm:
                sheet_gsm.add(code)
                if code in gsm:
                    patient=gsm[code]['patient']
                    same_row=bool(re.search(r'(?<![A-Za-z0-9])'+re.escape(patient)+r'(?![A-Za-z0-9])',joined))
                    found[code].append({'sheet':sheet,'row':index,'same_row_patient':same_row})
        inventory.append({'sheet':sheet,'nonempty_rows':len(rows),'gsm_total':len(sheet_gsm),
                          'crosswalk_gsm_found':len(sheet_gsm & set(gsm))})
        if len(sheet_gsm & set(gsm))>=140:s2_no_gsm=untagged
    absent=sorted(set(gsm)-set(found))
    if absent:raise ValueError('STOP: '+str(len(absent))+' S2 GSM not independently found in official XLSX')
    exact=sum(any(hit['same_row_patient'] for hit in row) for row in found.values())
    report={
      'state':'PASS_S2_GSM_IN_SOURCE__HOLD_CAUSAL_INFERENCE',
      'source_sha256':original_sha,'source_GSM_found':len(found),
      'source_same_row_patient_explicit':exact,
      'patient_source_verdict':'EXPLICIT_SAME_ROW' if exact==146 else 'REVIEW_GROUPED_ROWS_NOT_PROVEN',
      'workbook_sheets':inventory,
      'S2_no_GSM_row_candidates':s2_no_gsm,
      'cohort':{'GEO_HS_ROI':132,'GEO_AC_ROI':14,'manuscript_HS_ROI':156,
                'not_in_GSE294009':24,'prior_24_original_source':'UNVERIFIED'},
      'HS12_discordant_GSM':[r['gsm'] for r in disc],
      'scenarios':counts,'site_stage_preparation':site_rows,
      'compartment_matched_donor_coverage':coverage,
      'limitations':['Donors not ROIs are independent','Partial stage/site/preparation confounding',
                     'No valid p/FDR or treatment-response inference','Historical 24 require source provenance',
                     'HS12 stage mismatch unresolved']
    }
    (ROOT/'HS_SOURCE_STAGE_SITE_QC.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'state':report['state'],'original_verified_GSM':len(found),
                      'explicit_patient_same_row':exact,'scenarios':counts,
                      'site_stage_preparation':site_rows,'compartments':coverage,
                      'workbook_sheets':inventory,'S2_no_GSM_row_candidates':s2_no_gsm},ensure_ascii=False))
if __name__=='__main__':main()
