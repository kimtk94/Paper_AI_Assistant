#!/usr/bin/env python3
"""Public provenance scan for the 2023 HS GeoMx study's 24 historical ROI.
A citation, GDS search hit, or S2 metadata row is not a raw count matrix.
"""
import argparse, datetime as dt, json, re, urllib.parse, urllib.request
import urllib.error, xml.etree.ElementTree as ET
from pathlib import Path

PMID="37590958"
DOI="10.1093/bjd/ljad274"
TITLE="Spatial transcriptomic profiling reveals the pathophysiology of early-stage hidradenitis suppurativa"
UA="HS-GeoMx-Source-Audit/1.0 academic source verification"
GSM=re.compile(r"\bGSM\d{6,9}\b",re.I)
GSE=re.compile(r"\bGSE\d{5,9}\b",re.I)
MAX=2500000

def fetch(name,url):
    item={"name":name,"url":url}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json,text/xml,text/html,*/*"})
        with urllib.request.urlopen(req,timeout=16) as f:
            b=f.read(MAX+1)
            item.update(http=f.status,bytes=len(b),type=f.headers.get("Content-Type",""))
            if len(b)>MAX:
                item["status"]="TOO_LARGE";return item,""
            item["status"]="OK"
            return item,b.decode("utf8",errors="replace")
    except (urllib.error.URLError,TimeoutError,ValueError,OSError) as e:
        item.update(status="ERROR",reason=str(e)[:180]);return item,""

def safe_json(s):
    try:return json.loads(s)
    except (ValueError,TypeError):return {}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="output/hs_spatial_historic24")
    a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    targets=[
      ("EuropePMC_record","https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="+urllib.parse.quote("EXT_ID:"+PMID+" AND SRC:MED")+"&format=json"),
      ("EuropePMC_fulltext","https://www.ebi.ac.uk/europepmc/webservices/rest/MED/"+PMID+"/fullTextXML"),
      ("Crossref","https://api.crossref.org/works/"+DOI),
      ("PubMed_GEO_elink","https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gds&id="+PMID),
      ("Original_BJD","https://academic.oup.com/bjd/article/189/5/643/7244647"),
      ("GEO_GSE294009","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294009"),
    ]
    for i,term in enumerate([f'"{TITLE}"[All Fields]', '"Joh CSY"[Author] AND "hidradenitis"[All Fields]', '"ljad274"[All Fields]'],1):
        targets.append(("GDS_search_"+str(i),"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json&retmax=30&term="+urllib.parse.quote(term)))
    report={
      "date":dt.datetime.now(dt.timezone.utc).isoformat(),
      "source_paper":{"doi":DOI,"pmid":PMID,"title":TITLE},
      "historical_roi_metadata":{"HS41":6,"HS42":6,"HS43":12},
      "historical_24_GSM_in_2026_S2":0,
      "endpoints":[],"GDS_candidate_ids":[],"PubMed_GEO_linked_ids":[],
      "crossref_links":[],"article_metadata":[],"accession_mentions":[],
      "verified_historical_24_raw_roi_files":0,
      "gate":"HOLD_NO_VERIFIED_PAPER_LINKED_2023_ROI_QUANT"
    }
    for name,url in targets:
        row,body=fetch(name,url);report["endpoints"].append(row)
        if not body:continue
        if name=="EuropePMC_record":
            for x in safe_json(body).get("resultList",{}).get("result",[])[:3]:
                report["article_metadata"].append({"pmid":x.get("id"),"pmcid":x.get("pmcid"),"doi":x.get("doi"),"title":x.get("title")})
        if name=="Crossref":
            ent=safe_json(body).get("message",{})
            report["crossref_links"]=[{"URL":q.get("URL"),"type":q.get("content-type")} for q in ent.get("link",[])][:12]
        if name=="PubMed_GEO_elink":
            try:
                root=ET.fromstring(body)
                report["PubMed_GEO_linked_ids"]=sorted(set(x.text for d in root.findall(".//LinkSetDb") if d.findtext("DbTo")=="gds" for x in d.findall(".//Link/Id") if x.text))
            except ET.ParseError:pass
        if name.startswith("GDS_search"):
            d=safe_json(body).get("esearchresult",{})
            report["GDS_candidate_ids"].append({"query":name,"ids":d.get("idlist",[]),"count":d.get("count","0"),"verified_same_paper":False})
        report["accession_mentions"].append({"endpoint":name,"gse":sorted(set(GSE.findall(body)))[:20],"gsm":sorted(set(GSM.findall(body)))[:20],
                                             "meaning":"unverified mentions are NOT 2023 expression evidence"})
    (out/"HISTORICAL_24_PUBLIC_SOURCE_AUDIT.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    md=["# BJD 2023 HS GeoMx historic 24 ROI — public data audit","",
        "**Status: HOLD — no verified original 24-ROI expression file.**","",
        "Paper DOI: https://doi.org/"+DOI+" (PMID "+PMID+")",
        "The 2026 author S2 links HS41=6, HS42=6 and HS43=12 ROI to the 2023 study without GSM accessions.",
        "This is a statement about this source audit, not proof the data do not exist anywhere.","",
        "| Endpoint | Status |","|---|---|"]
    for x in report["endpoints"]:md.append("| "+x["name"]+" | "+x["status"]+" |")
    md+=["","GDS candidates and repository mentions must be cross-verified by DOI/PMID plus exact ROI identifiers and expression files.",
         "Never mix 132 new GEO counts with 24 missing raw expression columns.",
         "Next gate: direct raw DCC or gene-by-ROI file + documented probe panel/LOQ/Q3 normalization + donor mapping.",
         "Until then permit 132 HS GEO ROI descriptive donor-level analyses only; no pooled 156-ROI inferential DE.",""]
    (out/"HISTORICAL_24_PUBLIC_SOURCE_AUDIT.md").write_text("\n".join(md),encoding="utf8")
    print(json.dumps({"gate":report["gate"],"endpoints":{x["name"]:x["status"] for x in report["endpoints"]},
        "GDS_candidates":report["GDS_candidate_ids"],"PubMed_GEO_linked_ids":report["PubMed_GEO_linked_ids"],
        "crossref_links":report["crossref_links"],"EuropePMC_article":report["article_metadata"]},ensure_ascii=False))

if __name__=="__main__":main()
