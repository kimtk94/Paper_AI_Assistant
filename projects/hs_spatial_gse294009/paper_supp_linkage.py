#!/usr/bin/env python3
"""Find author-published GEO GSM -> ROI/donor linkage from public article supplements.
No donor guesses: only explicit table content should authorize mapping.
"""
import argparse, datetime as dt, json, re, urllib.request, urllib.error, xml.etree.ElementTree as ET, hashlib, zipfile, io
from pathlib import Path

PMCID='PMC13552446'; ROOT='https://pmc.ncbi.nlm.nih.gov/articles/'+PMCID
UA='Research-Reproducibility/1.0 (public metadata; provenance-first)'
URLS=[
 'https://www.ebi.ac.uk/europepmc/webservices/rest/'+PMCID+'/fullTextXML',
 'https://www.ebi.ac.uk/europepmc/webservices/rest/'+PMCID+'/supplementaryFiles',
 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/'+PMCID+'/unicode',
 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id='+PMCID,
 ROOT+'/',
]

def fetch(url,limit=12_000_000):
    try:
        req=urllib.request.Request(url,headers={'User-Agent':UA,'Accept':'application/xml,text/xml,application/zip,text/html,*/*'})
        with urllib.request.urlopen(req,timeout=40) as r:
            blob=r.read(limit+1);ret={'status':'OK','http':r.status,'bytes':len(blob),'content_type':r.headers.get('content-type',''),'url':url}
            if len(blob)>limit:ret['status']='TOO_LARGE';return None,ret
            ret['sha256']=hashlib.sha256(blob).hexdigest()
            return blob,ret
    except Exception as e:return None,{'status':'ERROR','url':url,'reason':str(e)[:300]}

def describe_xml(blob):
    root=ET.fromstring(blob)
    hits=[]; children=[]
    for node in root.iter():
        tag=node.tag.split('}')[-1]
        if tag in ('supplementary-material','inline-supplementary-material','supplement','ext-link','media','table-wrap','xref'):
            attrs={k.split('}')[-1]:v for k,v in node.attrib.items()}
            txt=' '.join(' '.join(node.itertext()).split())[:220]
            if tag in ('supplementary-material','inline-supplementary-material','media','ext-link') or ('supp' in txt.lower() and tag=='table-wrap'):
                hits.append({'tag':tag,'attrs':attrs,'text_preview':txt})
    return {'root':root.tag,'hits':hits[:220],'n_hits':len(hits)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='output/hs_spatial_supp');args=ap.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    report={'pmcid':PMCID,'date':dt.datetime.now(dt.timezone.utc).isoformat(),'source_urls':[], 'supplementary_candidates':[],'decision':'HOLD_NO_VERIFIED_DONOR_MAPPING'}
    blobs=[]
    for u in URLS:
        blob,status=fetch(u); report['source_urls'].append(status)
        if blob and blob.lstrip().startswith(b'<'):
            try:
                parsed=describe_xml(blob)
                report.setdefault('xml_documents',[]).append({'url':u,**parsed})
                if parsed['n_hits']>0:blobs.append((u,parsed))
            except Exception as e:report.setdefault('xml_errors',[]).append({'url':u,'error':str(e)[:240]})
        if blob and ('fullTextXML' in u):
            (out/'PMC13552446_article.xml').write_bytes(blob)
    candidates=set()
    for u,parsed in blobs:
        for z in parsed['hits']:
            for v in z['attrs'].values():
                if re.search(r'\.(?:xlsx|xls|csv|tsv|txt|zip|pdf|docx?|rds)(?:\?.*)?$',v,re.I):candidates.add(v)
    # The article may expose a specific supplementary table with an explicit GSM whitelist.
    for path in sorted(candidates):
        name=path.split('/')[-1].split('?')[0]
        if not name:continue
        urls=[path] if path.startswith('http') else [ROOT+'/bin/'+name,'https://pmc.ncbi.nlm.nih.gov/articles/instance/13552446/bin/'+name]
        cand={'name':name,'href':path,'sources':[]}
        for u in urls:
            blob,status=fetch(u,limit=12_000_000);cand['sources'].append(status)
            if blob and blob[:4] == b'PK\x03\x04' and re.search(r'\.xlsx$|\.zip$',name,re.I):
                # save supplementary workbook in raw untouched form for subsequent field/whitelist review
                (out/name).write_bytes(blob);cand['saved']=name;break
            elif blob and re.search(r'\.(csv|tsv|txt)$',name,re.I) and len(blob)<4_000_000:
                (out/name).write_bytes(blob);cand['saved']=name;break
            # Do not analyze PDF here: PDF needs a page-by-page source review.
        report['supplementary_candidates'].append(cand)
    # Europe PMC provides an authenticated/public supplementary archive in one ZIP.
    # Previous cap (12 MB) was too small; download once, extract only official tabular supplements.
    archive_url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+PMCID+'/supplementaryFiles'
    archive,stat=fetch(archive_url,limit=70_000_000)
    report['official_supplement_archive']=stat
    if archive and zipfile.is_zipfile(io.BytesIO(archive)):
        with zipfile.ZipFile(io.BytesIO(archive)) as z:
            report['official_archive_inventory']=[{'name':v.filename,'bytes':v.file_size,'compressed_bytes':v.compress_size} for v in z.infolist()][:100]
            for item in z.infolist():
                # Do not extract PDF here. Parsing PDF requires page-level review elsewhere.
                if item.filename.lower().endswith(('.xlsx','.csv','.tsv','.txt')) and item.file_size<=10_000_000:
                    name=item.filename.rsplit('/',1)[-1]
                    if not name:continue
                    raw=z.read(item)
                    (out/name).write_bytes(raw)
                    report.setdefault('official_saved_tabular',[]).append({'name':name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()})
    report['candidate_count']=len(candidates)
    report['saved_tabular_files']=[x['saved'] for x in report['supplementary_candidates'] if 'saved' in x]
    # Availability of a supplementary Table S2 does not mean it contains donor IDs; no automatic GO.
    (out/'SUPPLEMENT_DISCOVERY.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    print(json.dumps({'pmcid':PMCID,'source_statuses':[(x['url'].split('/')[-1],x['status']) for x in report['source_urls']], 'candidate_count':len(candidates),'saved':report['saved_tabular_files'],'official_tabular':report.get('official_saved_tabular',[]),'decision':report['decision']},ensure_ascii=False))

if __name__=='__main__':main()