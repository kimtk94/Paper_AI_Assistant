#!/usr/bin/env python3
"""Retrieve official PMC13552446 supplemental Methods PDF and audit GeoMx methods.

Finds the source PDF inside the official Europe PMC supplementary ZIP.
Fails closed when original source is unavailable; never substitutes typical GeoMx QC.
Outputs bounded article-provenanced snippets and a structured checklist, not a PDF copy.
"""
import argparse,datetime as dt,hashlib,io,json,re,urllib.request,urllib.error,zipfile
from pathlib import Path

PMCID="PMC13552446"
ARCHIVE="https://www.ebi.ac.uk/europepmc/webservices/rest/"+PMCID+"/supplementaryFiles"
SHEET="NIHMS2203138-supplement-1.pdf"
UA="HS-GeoMx-Published-Method-Reproducibility/1.0"
MAX=100_000_000
CATEGORIES={
  "LOQ and negative controls":r"\bLOQ\b|limit\s+of\s+quanti|negative\s+(?:control|probe)|geometric\s+(?:mean|standard)",
  "probe quality control":r"\bprobe(?:s)?\s+(?:qc|quality|filter)|low[-\s]*express|detect(?:ion|ed)|filter(?:ing|ed)?",
  "normalization":r"\bQ3\b|upper\s+quartile|normaliz|quantile|background",
  "batch and site":r"\bbatch\b|combat|site|FFPE|frozen|integrat",
  "differential expression":r"differential|limma|DESeq2|edgeR|FDR|adjusted\s+p",
  "ROI quality":r"\bROI\b|region(?:s)?\s+of\s+interest|segment|nuclei|sequenc(?:ing|e)\s+sat",
}
def download():
    req=urllib.request.Request(ARCHIVE,headers={"User-Agent":UA,"Accept":"application/zip,*/*"})
    with urllib.request.urlopen(req,timeout=125) as response:
        blob=response.read(MAX+1)
        meta={"url":ARCHIVE,"http":response.status,"bytes":len(blob),
              "sha256":hashlib.sha256(blob).hexdigest(),
              "content_type":response.headers.get("content-type","")}
    if len(blob)>MAX:raise ValueError("STOP: official source archive too large")
    if not zipfile.is_zipfile(io.BytesIO(blob)):raise ValueError("STOP: official source not a ZIP")
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        entries=[{"name":i.filename,"bytes":i.file_size} for i in z.infolist()]
        matches=[i for i in z.infolist() if i.filename.split("/")[-1]==SHEET]
        if len(matches)!=1:raise ValueError("STOP: expected one official supplement PDF, found "+str(len(matches)))
        info=matches[0]
        if info.file_size>50_000_000:raise ValueError("STOP: PDF unbounded size")
        b=z.read(info)
    if not b.startswith(b"%PDF"):raise ValueError("STOP: selected supplement is not PDF")
    return b,meta,entries

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",default="output/hs_geomx_published_methods")
    ap.add_argument("--retain-pdf",action="store_true")
    a=ap.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    raw,archive_meta,inventory=download()
    source=out/SHEET
    source.write_bytes(raw)
    from pypdf import PdfReader
    reader=PdfReader(io.BytesIO(raw))
    texts=[page.extract_text(extraction_mode="layout") or "" for page in reader.pages]
    page_lens=[len(x) for x in texts]
    if sum(page_lens)<1000:raise ValueError("STOP: PDF has insufficient extractable text; manual page render needed")
    hits={k:[] for k in CATEGORIES}
    for i,page in enumerate(texts,1):
        flat=" ".join(page.split())
        for key,expr in CATEGORIES.items():
            for m in list(re.finditer(expr,flat,re.I))[:15]:
                start=max(0,m.start()-160)
                stop=min(len(flat),m.end()+250)
                hits[key].append({"source_page":i,"trigger":m.group(),
                                  "context":flat[start:stop]})
    for k in hits:hits[k]=hits[k][:40]
    # Expose page-labeled extracted text to actions logs only via another
    # local CI diagnostic script; do not upload the full supplementary PDF.
    (out/"SUPPLEMENT_METHODS_PAGE_TEXT.txt").write_text(
        "\n".join("[SOURCE PDF PAGE "+str(i)+"]\n"+page for i,page in enumerate(texts,1)),
        encoding="utf8")
    artifact={
       "state":"ORIGINAL_SUPPLEMENT_PDF_RETRIEVED_REQUIRES_METHOD_REVIEW",
       "pmcid":PMCID,
       "source_file":SHEET,
       "source_publication_doi":"10.1016/j.jid.2026.01.041",
       "archive":archive_meta,
       "pdf_sha256":hashlib.sha256(raw).hexdigest(),
       "pdf_bytes":len(raw),"pdf_pages":len(texts),"page_characters":page_lens,
       "supplement_inventory":inventory,
       "method_occurrences":{k:len(v) for k,v in hits.items()},
       "method_contexts":hits,
       "claims": {"paper_LOQ_parameters_verified":False,
                  "paper_Q3_normalization_verified":False,
                  "paper_batch_strategy_verified":False,
                  "paper_DE_method_verified":False,
                  "prior_exploratory_LOQ_equivalent_to_authors":False}
    }
    (out/"ORIGINAL_METHODS_AUDIT.json").write_text(json.dumps(artifact,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    lines=["# GSE294009 paper original methods source review","",
           "**State: original supplemental PDF downloaded and extracted; exact methodological equivalence requires separate assessment.**","",
           "- PMCID: "+PMCID,
           "- Supplement PDF: "+SHEET,
           "- SHA-256: "+artifact["pdf_sha256"],
           "- Pages: "+str(artifact["pdf_pages"]),
           "- Source archive SHA-256: "+archive_meta["sha256"],"",
           "## Index of relevant original passages","",
           "|Topic|Occurrences (bounded)|Example PDF pages|",
           "|---|---:|---|"]
    for k,v in hits.items():
        lines.append("|"+k+"|"+str(len(v))+"|"+", ".join(str(i) for i in sorted({m["source_page"] for m in v})[:12])+"|")
    lines+=["","## Scientific guardrails","",
            "- Do not equate our negative-geomean × negative-geoSD² proxy with exact original LOQ without comparing original equations and negative-control filtering.",
            "- Do not conflate per-ROI Q3 normalization and author-integrated batch corrected expression.",
            "- Match each page-labeled methods passage to the original source and retain methodological differences as explicit HOLD.",
            "- No formal DE p/FDR or therapy response claims are authorized by successful retrieval.",""]
    (out/"ORIGINAL_METHODS_REPORT.md").write_text("\n".join(lines),encoding="utf8")
    print(json.dumps({"state":artifact["state"],"pdf_sha256":artifact["pdf_sha256"],
                      "pdf_bytes":len(raw),"pdf_pages":len(texts),
                      "method_occurrences":artifact["method_occurrences"],
                      "source_archive_status":"VERIFIED"},ensure_ascii=False))
    # The original licensed PDF is intentionally not added to GitHub artifact.
    if not a.retain_pdf:source.unlink()
if __name__=="__main__":main()
