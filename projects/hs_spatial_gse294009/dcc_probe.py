#!/usr/bin/env python3
"""Inspect public GeoMx DCC probes and PKC gene mapping, no inferential statistics."""
import argparse, collections, csv, gzip, hashlib, io, json, re, tarfile, xml.etree.ElementTree as ET
from pathlib import Path
from geo_file_probe import ACC, BASE, get

def xml_summary(b,limit=24):
    summary={'bytes_uncompressed':len(b),'preview':b[:650].decode('utf8','replace')}
    try:
        tree=ET.fromstring(b)
        ctr=collections.Counter();samples=collections.defaultdict(list);items=0
        for e in tree.iter():
            t=e.tag.split('}')[-1]
            ctr[t]+=1
            if len(samples[t])<2:
                samples[t].append({'attrib':dict(list(e.attrib.items())[:14]),'text':(e.text or '').strip()[:200], 'children':len(e)})
            items+=1
        summary.update({'kind':'XML','root_tag':tree.tag,'tag_counts':ctr.most_common(25),'tag_samples':dict(list(samples.items())[:limit]),'xml_elements':items})
    except Exception as e:summary.update({'kind':'NOT_XML','xml_error':str(e)[:200]})
    return summary

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='output/hs_spatial_gse294009');args=ap.parse_args();out=Path(args.out);out.mkdir(exist_ok=True,parents=True)
    audit={'accession':ACC,'stage':'DCC_FORMAT_AND_PKC_GENE_MAP_PROBE'}
    pkc_url=BASE+'/suppl/GSE294009_Hs_R_NGS_WTA_v1.0.pkc.gz'
    p,ps=get(pkc_url,limit=35_000_000,timeout=100);audit['pkc_download']=ps
    if p:
        try:
            decoded=gzip.decompress(p);audit['pkc']=xml_summary(decoded,limit=45)
            doc=json.loads(decoded)
            targets=doc.get('Targets',[])
            probes=[(t,pr) for t in targets for pr in t.get('Probes',[])]
            audit['pkc_json']={'top_level_keys':list(doc),'targets':len(targets),'probes':len(probes),'first_target_keys':list(targets[0]) if targets else [],'first_target':targets[0] if targets else {},'first_probe_keys':list(probes[0][1]) if probes else [], 'candidate_targets':[{ 'gene':t.get('DisplayName'),'probes':t.get('Probes',[])[:2]} for t in targets if t.get('DisplayName','').upper() in ('CXCL13','CCL20','MS4A1','CD79A','S100A8','S100A9','KRT17','IL17A','IL17F','CXCL8','JCHAIN','CFH','COL4A5','KLK6')][:15]}
            for c in ('RTS_ID','RTSID','ReporterID','ProbeID','ID','Name','Seq','RTS_Seq','RTS_Barcode'):
                audit['pkc_json'][f'{c}_examples']=[str(pr.get(c)) for _,pr in probes[:4]]
        except Exception as e:audit['pkc_error']=str(e)[:300]
    raw_url=BASE+'/suppl/'+ACC+'_RAW.tar';tar,t_status=get(raw_url,limit=20_000_000,timeout=140);audit['raw_download']=t_status
    if tar:
        try:
            with tarfile.open(fileobj=io.BytesIO(tar),mode='r:*') as archive:
                files=[x for x in archive.getmembers() if x.isfile()]
                audit['tar_members']=len(files)
                audit['tar_member_sample']=[x.name for x in files[:12]]
                cnt=collections.Counter(re.search(r'GSM\d+',v.name).group(0) if re.search(r'GSM\d+',v.name) else 'UNKNOWN' for v in files)
                audit['tar_unique_GSM']=len(cnt)-int('UNKNOWN' in cnt)
                audit['tar_dcc_count']=sum('.dcc' in x.name.lower() for x in files)
                picks=[]
                for m in files:
                    if len(picks)>=3:break
                    b=archive.extractfile(m).read()
                    if m.name.lower().endswith('.gz'):b=gzip.decompress(b)
                    picks.append({'file':m.name,'member_bytes':m.size,'summary':xml_summary(b,limit=35)})
                audit['dcc_samples']=picks
                if picks:
                    b=archive.extractfile(files[0]).read()
                    if files[0].name.endswith('.gz'):b=gzip.decompress(b)
                    dccs=dict((k,int(v)) for k,v in re.findall(rb'(?m)^(RTS\d+),(\d+)\r?$',b))
                    audit['dcc_example_code_count']=len(dccs)
                    audit['dcc_example_codes']=[k.decode() for k in list(dccs)[:15]]
                    if 'pkc_json' in audit:
                        for key in ['RTS_ID','RTSID','ReporterID','ProbeID','ID','RTS_Seq']:
                            vals=set(str(pr.get(key)) for t in doc.get('Targets',[]) for pr in t.get('Probes',[]))
                            audit['pkc_json'][key+'_matches_dcc_codes']=len({x.decode() for x in dccs}&vals)
        except Exception as e:audit['tar_error']=str(e)[:350]
    audit['decision']='PROBE_ONLY_NOT_YET_QUANTITATIVE';audit['cohort_warning']='ROI is nested in donors; cannot treat individual DCC as independent patient.'
    (out/'DCC_PKC_FORMAT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({'pkc_status':ps['status'],'pkc_root':audit.get('pkc',{}).get('root_tag'),'pkc_tags':audit.get('pkc',{}).get('tag_counts',[])[:10], 'tar_status':t_status['status'],'n_dcc':audit.get('tar_dcc_count'),'dcc_tags':audit.get('dcc_samples',[{}])[0].get('summary',{}).get('tag_counts',[])[:10]},ensure_ascii=False))

if __name__=='__main__':main()